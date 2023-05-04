import pandas as pd
import geopandas as gpd
import numpy as np
import shapely
from shapely.geometry import Polygon
from shapely.ops import unary_union
import matplotlib.pyplot as plt
from addict import Dict

import mph
import jpype.types as jtypes
from shapely.geometry import Polygon

import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import matplotlib as mpl

from SQDMetal.Utilities.QiskitShapelyRenderer import QiskitShapelyRenderer
from SQDMetal.Utilities.PVD_Shadows import PVD_Shadows
from SQDMetal.Utilities.QUtilities import QUtilities

from qiskit_metal.renderers.renderer_mpl.mpl_renderer import QMplRenderer
import pandas as pd
import geopandas as gpd


class COMSOL_Model:
    _engine = None
    
    def __init__(self, model_name):
        self.model_name = model_name
        self._model = COMSOL_Model._engine.create(model_name)

    @staticmethod
    def init_engine(num_cores = 2):
        if COMSOL_Model._engine == None:
            COMSOL_Model._engine = mph.start(cores=num_cores)

    @staticmethod
    def close_all_models():
        COMSOL_Model._engine.clear()

    def initialize_model(self, design, simulations, **kwargs):
        self.design = design

        self.chip_len = QUtilities.parse_value_length(design.chips['main']['size']['size_x'])
        self.chip_wid = QUtilities.parse_value_length(design.chips['main']['size']['size_y'])
        self.chip_thickness = np.abs(QUtilities.parse_value_length(design.chips['main']['size']['size_z']))
        self.chip_centre = [QUtilities.parse_value_length(design.chips['main']['size'][x]) for x in ['center_x', 'center_y', 'center_z']]

        self.pad_x = kwargs.get('pad_x', 0.5e-3)
        self.pad_y = kwargs.get('pad_y', 0.5e-3)
        self.pad_z = kwargs.get('pad_z', 0.5e-3)

        self.bottom_grounded = kwargs.get('bottom_grounded', False)

        self._conds = []
        self._cond_polys = []   #Stores (polygon object, selection index)
        self._fine_mesh = []
        self._num_sel = 0

        self._num_polys = 0

        self._model.java.component().create("comp1", True)

        self._model.java.component("comp1").geom().create("geom1", 3)
        self._model.java.component("comp1").mesh().create("mesh1")        

        #N.B. The ordering of the datasets here is the defauly naming scheme: dset1, dset2 etc...
        self._cur_dset_index = 1
        #List of physics directives used amongst the simulations (e.g. Electrostatic parameters, RF Frequency parameters, Magnetic Field parameters etc...)
        self._cur_physics = []
        #List of study names already in use
        self._cur_studies = []
        #List of solution names already in use
        self._cur_solns = []

        self._cur_sims = simulations
        
        # #Prepare simulation parameters
        # for cur_sim in self._cur_sims:
        #     cur_sim._prepare_simulation()

        #Create the main bounding area
        self._model.java.component("comp1").geom("geom1").lengthUnit("m")
        self._create_block_centre('blk_chip', self.chip_len,self.chip_wid,self.chip_thickness, self.chip_centre[0], self.chip_centre[1],-self.chip_thickness*0.5)
        if self.bottom_grounded:
            self._create_block_centre('blk_boundary', self.chip_len+2*self.pad_x,self.chip_wid+2*self.pad_y,self.chip_thickness+self.pad_z,
                                            self.chip_centre[0], self.chip_centre[1], (self.pad_z-self.chip_thickness)*0.5)
        else:
            self._create_block_centre('blk_boundary', self.chip_len+2*self.pad_x,self.chip_wid+2*self.pad_y,self.chip_thickness+2*self.pad_z,
                                                    self.chip_centre[0], self.chip_centre[1], -self.chip_thickness*0.5)
        
        #Create the selections to get face IDs for the exterior boundary
        self._model.java.component("comp1").geom("geom1").feature("blk_boundary").set("selresultshow", "bnd")
        self._model.java.component("comp1").geom("geom1").selection().create("cBndOuter", "CumulativeSelection")
        self._model.java.component("comp1").geom("geom1").selection("cBndOuter").label("bound_outer")
        self._model.java.component("comp1").geom("geom1").feature("blk_boundary").set("contributeto", "cBndOuter")
        self._sel_ext_boundaries = "geom1_cBndOuter_bnd"

        #Create workplane and subsequent metallic geometry...
        self._model.java.component("comp1").geom("geom1").feature().create("wp1", "WorkPlane")
        self._model.java.component("comp1").geom("geom1").feature("wp1").set("unite", jtypes.JBoolean(True)) #Unite all objects...

        # self._model.java.component("comp1").geom("geom1").run()

    def add_metallic(self, layer_id, **kwargs):
        '''
        Adds metallic conductors from the Klayout design onto the surface layer of the chip simulation. The idea is to supply the final Klayout layout object
        and the metallic polygons are added to the COMSOL simulation.
        Inputs:
            - kLayoutObj - A Klayout object (i.e. the object used when calling the function to save to a GDS file)
            - layer_id - The index of the layer from which to take the metallic polygons (e.g. in the old tutorials, that would be layer_photo)
            - cell_num - (Optional) The cell index in the Klayout object in which the layer resides. Default value is taken to be 0.
        '''

        #Fresh update on PVD profiles...
        pvd_shadows = PVD_Shadows(self.design)

        thresh = kwargs.get('threshold',  -1)

        qmpl = QiskitShapelyRenderer(None, self.design, None)
        gsdf = qmpl.get_net_coordinates()

        filt = gsdf.loc[(gsdf['layer'] == layer_id) & (gsdf['subtract'] == False)]
        if filt.shape[0] == 0:
            return

        unit_conv = QUtilities.get_units(self.design)
        
        fuse_threshold = kwargs.get('fuse_threshold', 1e-12)

        #Merge the metallic elements
        metal_polys = shapely.unary_union(filt['geometry'])
        metal_polys = shapely.affinity.scale(metal_polys, xfact=unit_conv, yfact=unit_conv, origin=(0,0))
        metal_polys = self._fuse_polygons_threshold(metal_polys, fuse_threshold)
        #Calculate the individual evaporated elements if required
        evap_mode = kwargs.get('evap_mode', 'separate_delete_below')
        group_by_evaporations = kwargs.get('group_by_evaporations', False)
        if group_by_evaporations and evap_mode != 'merge':
            metal_evap_polys_separate = pvd_shadows.get_all_shadows(metal_polys, layer_id, 'separate')
            #Convert all MultiPolygons into individual polygons...
            if not isinstance(metal_evap_polys_separate, list):
                metal_evap_polys_separate = [metal_evap_polys_separate]
        #Calculate evaporated shadows
        evap_trim = kwargs.get('evap_trim', 20e-9)
        metal_evap_polys = pvd_shadows.get_all_shadows(metal_polys, layer_id, evap_mode, layer_trim_length=evap_trim)
        #Convert all MultiPolygons into individual polygons...
        if not isinstance(metal_evap_polys, list):
            metal_evap_polys = [metal_evap_polys]
        metal_polys_all = []
        metal_sel_ids = []
        for m, cur_poly in enumerate(metal_evap_polys):
            if isinstance(cur_poly, shapely.geometry.multipolygon.MultiPolygon):
                temp_cur_metals = [x for x in cur_poly.geoms]
                metal_polys_all += temp_cur_metals
                num_polys = len(temp_cur_metals)
            else:
                temp_cur_metals = [cur_poly] #i.e. it's just a lonely Polygon object...
                metal_polys_all += temp_cur_metals
                num_polys = 1

            if group_by_evaporations and evap_mode != 'merge':
                #Collect the separate polygons that live in the current evaporation layer
                cur_polys_separate = metal_evap_polys_separate[m]
                if isinstance(cur_polys_separate, shapely.geometry.multipolygon.MultiPolygon):
                    cur_polys_separate = [x for x in cur_polys_separate.geoms]
                #Find the separate polygon in which the given polygon fits... 
                cur_sel_inds = []
                for cur_metal in temp_cur_metals:
                    for sep_piece_ind, cur_metal_piece in enumerate(cur_polys_separate):
                        if shapely.intersection(cur_metal, cur_metal_piece).area >= 0.99*cur_metal.area:
                            cur_sel_inds += [len(metal_sel_ids) + sep_piece_ind]
                            break
                metal_sel_ids += cur_sel_inds
            else:
                #Just enumerate to all separate metallic pieces...
                metal_sel_ids += [x for x in range(len(metal_sel_ids), len(metal_sel_ids) + num_polys)]

        metal_sel_obj_names = {}
        unique_ids = np.unique(metal_sel_ids)
        for cur_id in unique_ids:
            metal_sel_obj_names[cur_id] = []
        for m in range(len(metal_polys_all)):
            pol_name = "pol"+str(self._num_polys+m)
            #Convert coordinates into metres...
            cur_poly = np.array(metal_polys_all[m].exterior.coords[:])
            cur_poly = self._simplify_geom(cur_poly, thresh)
            sel_x, sel_y, sel_r = self._create_poly(pol_name, cur_poly[:-1])
            poly_interiors = metal_polys_all[m].interiors
            if len(poly_interiors) > 0:
                pol_name_ints = []
                for ind,cur_int in enumerate(poly_interiors):
                    pol_name_int = f"pol{self._num_polys+m}S{ind}"
                    pol_name_ints.append(pol_name_int)
                    #Convert coordinates into metres...
                    cur_poly_int = np.array(cur_int.coords[:])
                    cur_poly_int = self._simplify_geom(cur_poly_int, thresh)
                    sel_x, sel_y, sel_r = self._create_poly(pol_name_int, cur_poly_int[:-1])
                #Subtract interiors from the polygon
                diff_name = f"difInt{self._num_polys+m}"
                self._model.java.component("comp1").geom("geom1").feature("wp1").geom().create(diff_name, "Difference")
                self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(diff_name).selection("input").set(pol_name)
                self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(diff_name).selection("input2").set(*pol_name_ints)
                select_obj_name = diff_name
            else:
                select_obj_name = pol_name

            metal_sel_obj_names[metal_sel_ids[m]] += [select_obj_name]
        self._num_polys += len(metal_polys_all)

        #Settle the conductor selections...
        for cur_id in metal_sel_obj_names:
            cur_sel_index = len(self._conds)
            select_3D_name = self._setup_selection_boundaries(cur_sel_index, metal_sel_obj_names[cur_id])            
            self._conds.append( ("geom1_", select_3D_name) )  #Prefix for accessing out of the geometry node...
            #Could optimise this line...
            self._cond_polys.append( ( shapely.unary_union([metal_polys_all[m] for m in range(len(metal_polys_all)) if metal_sel_ids[m] == cur_id]), cur_sel_index ) )

    @staticmethod
    def get_poly_cardinality(poly):
        if isinstance(poly, shapely.geometry.multipolygon.MultiPolygon):
            return len(poly.geoms)
        else:
            return 1

    def _fuse_polygons_threshold(self, polys, threshold=1e-12):
        if not isinstance(polys, list):
            polys = [polys]
        lePolys = []
        for cur_poly in polys:
            if isinstance(cur_poly, shapely.geometry.multipolygon.MultiPolygon):
                lePolys += list(cur_poly.geoms)
            else:
                lePolys += [cur_poly]
        lePolys = [p.buffer(threshold, join_style=2, cap_style=3) for p in lePolys]
        return shapely.unary_union(lePolys).buffer(-threshold, join_style=2, cap_style=3)

    def fuse_all_metals(self, fuse_threshold = 1e-12):
        new_cond_polys = []
        inds_left = [x for x in range(len(self._cond_polys))]
        new_ind_groups = {}
        while len(inds_left) > 0:
            inds_to_pop = [inds_left[0]]
            cur_blobotron = self._cond_polys[inds_left[0]][0]
            other_poly_inds = inds_left[1:]
            while True:
                found_another_ind = False
                for cur_ind in other_poly_inds:
                    cur_check_poly = self._cond_polys[cur_ind][0]
                    cur_check = self._fuse_polygons_threshold([cur_blobotron, cur_check_poly], fuse_threshold)
                    if COMSOL_Model.get_poly_cardinality(cur_check) < COMSOL_Model.get_poly_cardinality(cur_blobotron) + COMSOL_Model.get_poly_cardinality(cur_check_poly):
                        cur_blobotron = cur_check
                        inds_to_pop += [cur_ind]
                        other_poly_inds.remove(cur_ind)
                        found_another_ind = True
                if not found_another_ind:
                    break
            new_ind = self._cond_polys[inds_left[0]][1]
            new_cond_polys += [(cur_blobotron, new_ind)]
            new_ind_groups[new_ind] = inds_to_pop 
            for cur_ind_to_remove in inds_to_pop:
                inds_left.remove(cur_ind_to_remove)
        self._cond_polys = new_cond_polys

        self._model.java.component("comp1").geom("geom1").run()

        new_conds = []
        for cur_ind in new_ind_groups:
            select_3D_name = f"selFused{cur_ind}"
            self._model.java.component("comp1").geom("geom1").create(select_3D_name, "UnionSelection")
            self._model.java.component("comp1").geom("geom1").feature(select_3D_name).label(select_3D_name)
            self._model.java.component("comp1").geom("geom1").feature(select_3D_name).set("entitydim", jtypes.JInt(2))
            self._model.java.component("comp1").geom("geom1").feature(select_3D_name).set("input", jtypes.JArray(jtypes.JString)([self._conds[x][1] for x in new_ind_groups[cur_ind]]))
            new_conds += [('geom1_', select_3D_name)]
        self._conds = new_conds


    def _simplify_geom(self, cur_poly, threshold):
        if threshold > 0:
            fin_poly = [cur_poly[0]]
            for m in range(1,len(cur_poly)):
                if np.sqrt((cur_poly[m][0]-cur_poly[m-1][0])**2+(cur_poly[m][1]-cur_poly[m-1][1])**2) > threshold:
                    fin_poly += [cur_poly[m]]
            cur_poly = fin_poly
        return cur_poly

    def _setup_selection_boundaries(self, index, select_obj_name):
        if not isinstance(select_obj_name, list):
            select_obj_name = [select_obj_name] 
        select_name = f"sel{index}"
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().create(select_name, "ExplicitSelection")
        for cur_sel_name in select_obj_name:
            self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(select_name).selection("selection").set(cur_sel_name, 1)
        select_3D_name = f"cond{index}"
        self._model.java.component("comp1").geom("geom1").create(select_3D_name, "UnionSelection")
        self._model.java.component("comp1").geom("geom1").feature(select_3D_name).label(select_3D_name)
        self._model.java.component("comp1").geom("geom1").feature(select_3D_name).set("entitydim", jtypes.JInt(2))
        self._model.java.component("comp1").geom("geom1").feature(select_3D_name).set("input", f"wp1_{select_name}")
        return select_3D_name

    def add_ground_plane(self, **kwargs):
        #Assumes that all masks are simply closed without any interior holes etc... Why would you need one anyway?
        thresh = kwargs.get('threshold', -1)

        qmpl = QiskitShapelyRenderer(None, self.design, None)
        gsdf = qmpl.get_net_coordinates()

        filt = gsdf.loc[gsdf['subtract'] == True]
        #Now create a plane sheet covering the entire chip
        pol_name = "polgndBase"
        xL = self.chip_centre[0] - self.chip_len*0.5
        xR = self.chip_centre[0] + self.chip_len*0.5
        yL = self.chip_centre[1] - self.chip_wid*0.5
        yR = self.chip_centre[1] + self.chip_wid*0.5
        self._create_poly(pol_name, [[xL,yL], [xR,yL], [xR,yR], [xL,yR]])
        if filt.shape[0] > 0:
            unit_conv = QUtilities.get_units(self.design)
                
            space_polys = shapely.unary_union(filt['geometry'])
            if isinstance(space_polys, shapely.geometry.multipolygon.MultiPolygon):
                space_polys = [x for x in space_polys.geoms]
            else:
                space_polys = [space_polys] #i.e. it's just a lonely Polygon object...
            space_polys = [shapely.affinity.scale(x, xfact=unit_conv, yfact=unit_conv, origin=(0,0)) for x in space_polys]
            pol_name_ints = []
            for ind,cur_int in enumerate(space_polys):
                pol_name_int = f"polGND{ind}"
                pol_name_ints.append(pol_name_int)
                #Convert coordinates into metres...
                cur_poly_int = np.array(cur_int.exterior.coords[:])
                cur_poly_int = self._simplify_geom(cur_poly_int, thresh)
                sel_x, sel_y, sel_r = self._create_poly(pol_name_int, cur_poly_int[:-1])
            #Subtract interiors from the polygon
            diff_name = f"difGND"
            self._model.java.component("comp1").geom("geom1").feature("wp1").geom().create(diff_name, "Difference")
            self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(diff_name).selection("input").set(pol_name)
            self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(diff_name).selection("input2").set(*pol_name_ints)
            #Setup the resulting selection...
            self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(diff_name).set("selresult", jtypes.JBoolean(True))
            self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(diff_name).set("selresultshow", "bnd")
            select_2D_name = "cselGND"
            self._model.java.component("comp1").geom("geom1").feature("wp1").geom().selection().create(select_2D_name, "CumulativeSelection")
            self._model.java.component("comp1").geom("geom1").feature("wp1").geom().selection(select_2D_name).label("Cumulative Selection Ground")
            self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(diff_name).set("contributeto", select_2D_name)
            select_3D_name = "condGND"
            self._model.java.component("comp1").geom("geom1").create(select_3D_name, "UnionSelection")
            self._model.java.component("comp1").geom("geom1").runPre(select_3D_name)
            self._model.java.component("comp1").geom("geom1").feature(select_3D_name).set("entitydim", jtypes.JInt(2))
            self._model.java.component("comp1").geom("geom1").feature(select_3D_name).set("input", jtypes.JArray(jtypes.JString)(["wp1_"+select_2D_name]))
            cur_select_index = len(self._conds)
        else:
            select_obj_name = pol_name
            cur_select_index = len(self._conds)
            select_3D_name = self._setup_selection_boundaries(cur_select_index, select_obj_name)            
        self._conds.append(("geom1_", select_3D_name))  #Prefix for accessing out of the geometry node...
        #
        poly_sheet = Polygon([[xL,yL], [xR,yL], [xR,yR], [xL,yR]])
        space_whole = shapely.unary_union(space_polys)
        self._cond_polys += [(poly_sheet.difference(space_whole), cur_select_index)]

    def _extract_poly_coords(self, geom):
        #Inspired by: https://stackoverflow.com/questions/21824157/how-to-extract-interior-polygon-coordinates-using-shapely
        if isinstance(geom, Polygon):
            exterior_coords = geom.exterior.coords[:]
            interior_coords = []
            for interior in geom.interiors:
                interior_coords += [interior.coords[:]]
        elif isinstance(geom, shapely.geometry.multipolygon.MultiPolygon):
            exterior_coords = []
            interior_coords = []
            for part in list(geom.geoms):
                extC, intC = self._extract_poly_coords(part)  # Recursive call
                exterior_coords += extC
                interior_coords += intC
        else:
            raise ValueError('Unhandled geometry type: ' + repr(geom.type))
        return [exterior_coords], [interior_coords]

    def _create_block_corner(self, name, size_x,size_y,size_z, pos_x,pos_y,pos_z):
        '''
        Creates a block in the geometry.
        Inputs:
            - name - Unique name of the geometry object
            - size_x,size_y,size_z - Dimensions of the block in the prescribed units
            - pos_x,pos_y,pos_z    - Corner position of the block in the prescribed units
        '''
        self._model.java.component("comp1").geom("geom1").create(name, "Block")
        self._model.java.component("comp1").geom("geom1").feature(name).set("base", "corner")
        self._model.java.component("comp1").geom("geom1").feature(name).set("size", jtypes.JArray(jtypes.JDouble)([size_x,size_y,size_z]))
        self._model.java.component("comp1").geom("geom1").feature(name).set("pos", jtypes.JArray(jtypes.JDouble)([pos_x,pos_y,pos_z]))
        self._model.java.component("comp1").geom("geom1").feature(name).set('selresult', 'on') #To generate automatic selections...

    def _create_block_centre(self, name, size_x,size_y,size_z, pos_x,pos_y,pos_z):
        '''
        Creates a block in the geometry.
        Inputs:
            - name - Unique name of the geometry object
            - size_x,size_y,size_z - Dimensions of the block in the prescribed units
            - pos_x,pos_y,pos_z    - Corner position of the block in the prescribed units
        '''
        self._model.java.component("comp1").geom("geom1").create(name, "Block")
        self._model.java.component("comp1").geom("geom1").feature(name).set("base", "center")
        self._model.java.component("comp1").geom("geom1").feature(name).set("size", jtypes.JArray(jtypes.JDouble)([size_x,size_y,size_z]))
        self._model.java.component("comp1").geom("geom1").feature(name).set("pos", jtypes.JArray(jtypes.JDouble)([pos_x,pos_y,pos_z]))
        self._model.java.component("comp1").geom("geom1").feature(name).set('selresult', 'on') #To generate automatic selections...


    def _create_poly(self, name, poly_coords):
        '''
        Creates a polygon on workplane wp1. The return value is (x,y,r) where (x,y) is a point inside the polygon and r is the radius such that a sphere
        fits inside the polygon...
        Inputs:
            - name - Unique name of the geometry object
            - poly_coords - Coordinates (doesn't have to close) given as a list of lists: [[x1,y1], [x2,y2], ...]
        '''
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().create(name, "Polygon")
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(name).set('selresult', 'on') #To generate automatic selections...
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(name).set('selresultshow', 'bnd')
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(name).set("source", "table")
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(name).set('table', jtypes.JArray(jtypes.JDouble,2)(poly_coords))    #This is much faster...
        return self._find_point_in_polygon(poly_coords)
        
    def _find_point_in_polygon(self, poly_coords):
        min_ind = 0
        for row_no,cur_pt in enumerate(poly_coords):
            if cur_pt[1] <= poly_coords[min_ind][1]:
                min_ind = row_no
        #Find point inside polygon...
        min_ind2 = (min_ind + 1) % len(poly_coords)
        vec1 = np.array([ poly_coords[min_ind2][0]-poly_coords[min_ind][0], poly_coords[min_ind2][1]-poly_coords[min_ind][1] ])
        vec2 = np.array([ poly_coords[min_ind-1][0]-poly_coords[min_ind][0], poly_coords[min_ind-1][1]-poly_coords[min_ind][1] ])
        vec1 = vec1/np.linalg.norm(vec1)
        vec2 = vec2/np.linalg.norm(vec2)
        if np.dot(vec1,vec2) < -0.999847695: # 179 degrees
            vec3 = np.array([0.0,1.0])  #Assumption of being the lower-most point...
        else:
            vec3 = vec1 + vec2
            vec3 = vec3/np.linalg.norm(vec3)
        epsilon_mov = 10e-9
        rad = np.sqrt(min(1-np.dot(vec1,vec3)**2, 1-np.dot(vec2,vec3)**2)) * 0.9 * epsilon_mov
        vec3 *= epsilon_mov
        return (vec3[0]+poly_coords[min_ind][0], vec3[1]+poly_coords[min_ind][1], rad)

    def _create_boundary_selection_sphere(self, radius, pos_x,pos_y,pos_z=0.0):
        '''
        Creates a selection in which all boundaries within the sphere are selected (after the geometry has been fully built). Return value is the selection
        name/ID that can be used to query for domains later via the function _get_selection_boundaries.
        Inputs:
            - radius - Radius of sphere
            - pos_x,pos_y,pos_z - Position of sphere's centre
        '''
        sel_name = 'sel' + str(self._num_sel)
        self._num_sel += 1
        self._model.java.selection().create(sel_name, 'Ball')
        self._model.java.selection(sel_name).set('entitydim', '2')
        self._model.java.selection(sel_name).set('condition', 'intersects')
        self._model.java.selection(sel_name).set('posx', jtypes.JDouble(pos_x))
        self._model.java.selection(sel_name).set('posy', jtypes.JDouble(pos_y))
        self._model.java.selection(sel_name).set('posz', jtypes.JDouble(pos_z))
        self._model.java.selection(sel_name).set('r', jtypes.JDouble(radius))
        return sel_name

    def _get_selection_boundaries(self, sel_name):
        '''
        Returns a list of integers pertaining to the boundaries within a spherical-selection of ID given by sel_name.
        '''
        return [x for x in self._model.java.selection(sel_name).entities(2)]



    def _create_material(self, name, rel_permit, rel_permea, conductivity):
        '''
        Creates a polygon on workplane wp1.
        Inputs:
            - name - Unique name of the material
            - rel_permit - Relative permittivity of material
            - rel_permea - Relative permeability of material
            - conductivity - Electrical conductivity of material (sigma) given in S/m
            - selected_domain - Selected domain specified via commands like: self._model.java.selection('geom1_blk_chip_dom').entities(3) or a string
                                naming the selection.
        '''
        self._model.java.component("comp1").material().create(name, "Common")
        self._model.java.component("comp1").material(name).propertyGroup("def").set("relpermittivity", jtypes.JDouble(rel_permit))
        self._model.java.component("comp1").material(name).propertyGroup("def").set("relpermeability", jtypes.JDouble(rel_permea))
        self._model.java.component("comp1").material(name).propertyGroup("def").set("electricconductivity", jtypes.JDouble(conductivity))

    def _get_dset_name(self):
        retStr = f"dset{self._cur_dset_index}"
        self._cur_dset_index += 1
        return retStr
    
    def _add_physics(self, new_name):
        cntr = 0
        test_name = new_name
        while test_name in self._cur_physics:
            test_name = new_name + str(cntr)
        self._cur_physics += [test_name]
        return test_name

    def _add_study(self, new_name):
        cntr = 0
        test_name = new_name
        while test_name in self._cur_studies:
            test_name = new_name + str(cntr)
        self._cur_studies += [test_name]
        return test_name

    def _add_solution(self, new_name):
        cntr = 0
        test_name = new_name
        while test_name in self._cur_solns:
            test_name = new_name + str(cntr)
        self._cur_solns += [test_name]
        return test_name

    def _get_java_comp(self):
        return self._model.java

    def build_geom_mater_elec_mesh(self, mesh_structure='Normal', skip_meshing=False):    #mesh_structure can be 'Fine'
        '''
        Builds geometry, sets up materials, sets up electromagnetic parameters/ports and builds the mesh.
        '''
        #Build geometry
        self._model.java.component("comp1").geom("geom1").run()
        #
        select_3D_name = f"condAll"
        self._model.java.component("comp1").geom("geom1").create(select_3D_name, "UnionSelection")
        self._model.java.component("comp1").geom("geom1").feature(select_3D_name).label(select_3D_name)
        self._model.java.component("comp1").geom("geom1").feature(select_3D_name).set("entitydim", jtypes.JInt(2))
        self._model.java.component("comp1").geom("geom1").feature(select_3D_name).set("input", jtypes.JArray(jtypes.JString)([x[1] for x in self._conds]))

        #Create materials
        self._create_material('Vacuum', 1.0, 1.0, 0.0)
        self._model.java.component("comp1").material('Vacuum').selection().all()
        self._create_material('Si', 11.7, 1.0, 0.0)     #Don't use 4.35e-4 for conductivity as it can cause issues with convergence in inductance simulations...
        self._model.java.component("comp1").material('Si').selection().set(self._model.java.selection('geom1_blk_chip_dom').entities(3))
        #Set the conductive boundaries to Aluminium (ignored by s-parameter and capacitance-matrix simulations due to PEC, but inductor simulations require conductivity)
        # cond_bounds = [self._get_selection_boundaries(x)[0] for x in self._conds]
        self._create_material('Metal', 0.0, 1.0, 3.77e7)
        self._model.java.component("comp1").material("Metal").selection().geom("geom1", 2)
        self._model.java.component("comp1").material("Metal").selection().named("geom1_condAll")

        for cur_sim in self._cur_sims:
            cur_sim._run_premesh()

        #Create mesh
        for mesh_ind, cur_fine_struct in enumerate(self._fine_mesh):
            cur_polys = self._get_selection_boundaries(cur_fine_struct[0])
            if len(cur_polys) == 0:
                continue
            mesh_name = "ftri" + str(mesh_ind)
            self._model.java.component("comp1").mesh("mesh1").create(mesh_name, "FreeTri")
            self._model.java.component("comp1").mesh("mesh1").feature(mesh_name).create("size1", "Size")
            self._model.java.component("comp1").mesh("mesh1").feature(mesh_name).selection().set(jtypes.JArray(jtypes.JInt)(cur_polys))
            self._model.java.component("comp1").mesh("mesh1").feature(mesh_name).feature("size1").set('custom', jtypes.JBoolean(True))
            self._model.java.component("comp1").mesh("mesh1").feature(mesh_name).feature("size1").set('hmaxactive', jtypes.JBoolean(True))
            self._model.java.component("comp1").mesh("mesh1").feature(mesh_name).feature("size1").set('hminactive', jtypes.JBoolean(True))
            self._model.java.component("comp1").mesh("mesh1").feature(mesh_name).feature("size1").set('hmin', jtypes.JDouble(cur_fine_struct[1]))
            self._model.java.component("comp1").mesh("mesh1").feature(mesh_name).feature("size1").set('hmax', jtypes.JDouble(cur_fine_struct[2]))
        # self._model.java.component("comp1").mesh("mesh1").feature("size").set("custom", jtypes.JBoolean(True))
        # self._model.java.component("comp1").mesh("mesh1").feature("size").set("hmax", jtypes.JDouble(8e-4))
        # self._model.java.component("comp1").mesh("mesh1").feature("size").set("hmin", jtypes.JDouble(10e-6))
        self._model.java.component("comp1").mesh("mesh1").create("ftet10", "FreeTet")
        self._model.java.component("comp1").mesh("mesh1").feature("ftet10").create("size1", "Size")
        if mesh_structure == 'Fine':
            self._model.java.component("comp1").mesh("mesh1").feature("ftet10").feature("size1").set("hauto", jtypes.JInt(4))
        if not skip_meshing:
            self._model.java.component("comp1").mesh("mesh1").run()

        #Prepare simulation parameters
        for cur_sim in self._cur_sims:
            cur_sim._prepare_simulation()

        #Activate the appropriate physics (to be solved)
        for cur_sim in self._cur_sims:
            cur_req_phys = cur_sim._get_required_physics()
            for cur_physic in self._cur_physics:
                cur_active = (cur_physic in cur_req_phys)
                if isinstance(cur_sim._sub_study, list):
                    for cur_substudy in cur_sim._sub_study:
                        self._model.java.study(cur_sim._study).feature(cur_substudy).activate(cur_physic, cur_active)
                else:                    
                    self._model.java.study(cur_sim._study).feature(cur_sim._sub_study).activate(cur_physic, cur_active)


    def save(self, file_name):
        self._model.save(file_name)


class COMSOL_Simulation_Base:
    def _prepare_simulation(self):
        pass
    def _run_premesh(self):
        pass
    def _get_required_physics(self):
        return []
    def run(self):
        pass
