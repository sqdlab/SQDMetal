# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0

import geopandas as gpd
import numpy as np
import shapely
from shapely.geometry import Polygon
import matplotlib.pyplot as plt

import mph
import jpype.types as jtypes


from SQDMetal.Utilities.QiskitShapelyRenderer import QiskitShapelyRenderer
from SQDMetal.Utilities.PVD_Shadows import PVD_Shadows
from SQDMetal.Utilities.QUtilities import QUtilities


from SQDMetal.Utilities.Materials import Material
from SQDMetal.Utilities.ShapelyEx import ShapelyEx

class COMSOL_Model:
    _engine = None
    
    def __init__(self, model_name):
        self.model_name = model_name
        self._model = COMSOL_Model._engine.create(model_name)

    @staticmethod
    def init_engine(num_cores = 2):
        if COMSOL_Model._engine is None:
            COMSOL_Model._engine = mph.start(cores=num_cores)

    @staticmethod
    def close_all_models():
        COMSOL_Model._engine.clear()

    def initialize_model(self, design, simulations, **kwargs):
        self.design = design

        restrict_rect = kwargs.get('segment_rectangle', None)   #Given as [x1,y1,x2,y2]
        if isinstance(restrict_rect, list) or isinstance(restrict_rect, np.ndarray) or isinstance(restrict_rect, tuple):
            self.restrict_rect = [min([restrict_rect[0], restrict_rect[2]]), min([restrict_rect[1], restrict_rect[3]]),
                                  max([restrict_rect[0], restrict_rect[2]]), max([restrict_rect[1], restrict_rect[3]])]
            self.chip_len = self.restrict_rect[2]-self.restrict_rect[0]
            self.chip_wid = self.restrict_rect[3]-self.restrict_rect[1]
            self.chip_centre = [(self.restrict_rect[0]+self.restrict_rect[2])*0.5,
                                (self.restrict_rect[1]+self.restrict_rect[3])*0.5,
                                QUtilities.parse_value_length(design.chips['main']['size']['center_z'])]
        else:
            self.chip_len = QUtilities.parse_value_length(design.chips['main']['size']['size_x'])
            self.chip_wid = QUtilities.parse_value_length(design.chips['main']['size']['size_y'])
            self.chip_centre = [QUtilities.parse_value_length(design.chips['main']['size'][x]) for x in ['center_x', 'center_y', 'center_z']]
            self.restrict_rect = [self.chip_centre[0]-0.5*self.chip_len, self.chip_centre[1]-0.5*self.chip_wid,
                                  self.chip_centre[0]+0.5*self.chip_len, self.chip_centre[1]+0.5*self.chip_wid]

        self.chip_thickness = np.abs(QUtilities.parse_value_length(design.chips['main']['size']['size_z']))

        self.pad_x = kwargs.get('pad_x', 0.5e-3)
        self.pad_y = kwargs.get('pad_y', 0.5e-3)
        self.pad_z = kwargs.get('pad_z', 0.5e-3)

        self._include_loss_tangents = kwargs.get('include_loss_tangents', False)

        self.bottom_grounded = kwargs.get('bottom_grounded', False)

        self._conds = []
        self._cond_polys = []   #Stores (polygon object, selection index)
        self._fine_mesh = []
        self._num_sel = 0

        self._num_polys = 0

        self._model.java.component().create("comp1", True)

        self._model.java.component("comp1").geom().create("geom1", 3)
        self._model.java.component("comp1").mesh().create("mesh1")        

        #N.B. The ordering of the datasets here is the default naming scheme: dset1, dset2 etc...
        self._cur_dset_index = 1
        #List of physics directives used amongst the simulations (e.g. Electrostatic parameters, RF Frequency parameters, Magnetic Field parameters etc...)
        self._cur_physics = []
        #List of study names already in use
        self._cur_studies = []
        #List of solution names already in use
        self._cur_solns = []

        self._cur_sims = simulations

        self._resolution = kwargs.get('resolution', 4)
        
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
        Adds metallic conductors from the Qiskit-Metal design object onto the surface layer of the chip simulation. If the particular layer has
        fancy PVD evaporation steps, the added metallic layer will account for said steps and merge the final result. In addition, all metallic
        elements that are contiguous are merged into single blobs.

        Inputs:
            - layer_id - The index of the layer from which to take the metallic polygons. If this is given as a LIST, then the metals in the specified
                         layer (0 being ground plane) will be fused and then added into COMSOL.
            - threshold - (Optional) Defaults to -1. This is the threshold in metres, below which consecutive vertices along a given polygon are
                          combined into a single vertex. This simplification helps with meshing as COMSOL will not overdo the meshing. If this
                          argument is negative, the argument is ignored.
            - fuse_threshold - (Optional) Defaults to 1e-12. This is the minimum distance between metallic elements, below which they are considered
                               to be a single polygon and thus, the polygons are merged with the gap filled. This accounts for floating-point errors
                               that make adjacent elements fail to merge as a single element, due to infinitesimal gaps between them.
            - evap_mode - (Optional) Defaults to 'separate_delete_below'. These are the methods upon which to separate or merge overlapping elements
                          across multiple evaporation steps. See documentation on PVD_Shadows for more details on the available options.
            - group_by_evaporations - (Optional) Defaults to False. If set to True, if elements on a particular evaporation step are separated due
                                      to the given evap_mode, they will still be selected as a part of the same conductor (useful for example, in
                                      capacitance matrix simulations).
            - evap_trim - (Optional) Defaults to 20e-9. This is the trimming distance used in certain evap_mode profiles. See documentation on
                          PVD_Shadows for more details on its definition.
            - multilayer_fuse - (Optional) Defaults to False. Flattens everything into a single layer (careful when using this with evap_mode).
            - smooth_radius - (Optional) Defaults to 0. If above 0, then the corners of the metallic surface will be smoothed via this radius. Only works
                              if multilayer_fuse is set to True.
            - ground_cutout - (Optional) MUST BE SPECIFIED if layer_id is 0. It is a tuple (x1, x2, y1, y2) for the x and y bounds
        '''
        thresh = kwargs.get('threshold', -1)
        kwargs['restrict_rect'] = self.restrict_rect
        xL = self.chip_centre[0] - self.chip_len*0.5
        xR = self.chip_centre[0] + self.chip_len*0.5
        yL = self.chip_centre[1] - self.chip_wid*0.5
        yR = self.chip_centre[1] + self.chip_wid*0.5
        kwargs['ground_cutout'] = (xL, xR, yL, yR)
        kwargs['resolution'] = self._resolution
        metal_polys_all, metal_sel_ids = QUtilities.get_metals_in_layer(self.design, layer_id, **kwargs)

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

    def _add_cond(self, poly_coords):
        cur_sel_index = self._num_polys

        selObjName = f"Uclip{cur_sel_index}"
        self._create_poly(selObjName, poly_coords)
        select_3D_name = self._setup_selection_boundaries(len(self._conds), selObjName)
        self._conds.append( ("geom1_", select_3D_name) )  #Prefix for accessing out of the geometry node...
        self._cond_polys.append((shapely.Polygon(poly_coords), cur_sel_index))

        self._num_polys += 1

    @staticmethod
    def get_poly_cardinality(poly):
        if isinstance(poly, shapely.geometry.multipolygon.MultiPolygon):
            return len(poly.geoms)
        else:
            return 1

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
                    cur_check = ShapelyEx.fuse_polygons_threshold([cur_blobotron, cur_check_poly], fuse_threshold)
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

    def reorder_conds_by_comps(self, comp_list):
        qmpl = QiskitShapelyRenderer(None, self.design, None)
        gsdf = qmpl.get_net_coordinates(self._resolution)
        pvdSh = PVD_Shadows(self.design)
        unit_conv = QUtilities.get_units(self.design)
        
        temp_conds_and_polys = []
        for cur_comp in comp_list:
            #Construct the net PVD shadowed polygon
            polys = gsdf.loc[(gsdf['component'] == self.design.components[cur_comp].id) & (~gsdf['subtract'])]['geometry'].to_list()
            if len(polys) > 1:
                poly = shapely.unary_union(polys)
                assert isinstance(poly, shapely.Polygon), f"Component \'{cur_comp}\' is not a contiguous polygonal object."
            else:
                poly = polys[0]
            poly = shapely.affinity.scale(poly, xfact=unit_conv, yfact=unit_conv, origin=(0,0))
            poly = pvdSh.get_all_shadows(poly, self.design.components[cur_comp].options.layer, 'merge')
            for m in range(len(self._conds)):
                if ShapelyEx.chk_within(poly, self._cond_polys[m][0]):
                    cond = self._conds.pop(m)
                    cond_poly = self._cond_polys.pop(m)
                    temp_conds_and_polys += [(cond, cond_poly[0])]
                    break
        #Just add unlabelled polygons...
        residual_conds = [(self._conds[m], self._cond_polys[m][0]) for m in range(len(self._conds))]
        if len(residual_conds) > 0:
            print('SQDMetal Warning: Residual conductors found that were unlabelled.')
        temp_conds_and_polys += residual_conds

        self._conds = []
        self._cond_polys = []
        for m in range(len(temp_conds_and_polys)):
            self._conds += [temp_conds_and_polys[m][0]]
            self._cond_polys += [(temp_conds_and_polys[m][1], m+1)] 

    def _simplify_geom(self, cur_poly, threshold):
        if threshold > 0:
            fin_poly = [cur_poly[0]]
            for m in range(1,len(cur_poly)):
                if np.sqrt((cur_poly[m][0]-cur_poly[m-1][0])**2+(cur_poly[m][1]-cur_poly[m-1][1])**2) > threshold:
                    fin_poly += [cur_poly[m]]
            cur_poly = fin_poly
        return cur_poly

    def _setup_selection_boundaries(self, index, select_obj_name, sel_3D_prefix='cond'):
        if not isinstance(select_obj_name, list):
            select_obj_name = [select_obj_name] 
        select_name = f"sel{sel_3D_prefix}{index}"
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().create(select_name, "ExplicitSelection")
        for cur_sel_name in select_obj_name:
            self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(select_name).selection("selection").set(cur_sel_name, 1)
        select_3D_name = f"{sel_3D_prefix}{index}"
        self._model.java.component("comp1").geom("geom1").create(select_3D_name, "UnionSelection")
        self._model.java.component("comp1").geom("geom1").feature(select_3D_name).label(select_3D_name)
        self._model.java.component("comp1").geom("geom1").feature(select_3D_name).set("entitydim", jtypes.JInt(2))
        self._model.java.component("comp1").geom("geom1").feature(select_3D_name).set("input", f"wp1_{select_name}")
        return select_3D_name

    def _get_gnd_plane_gaps(self, fuse_threshold):
        qmpl = QiskitShapelyRenderer(None, self.design, None)
        gsdf = qmpl.get_net_coordinates(self._resolution)

        filt = gsdf.loc[gsdf['subtract']]

        unit_conv = QUtilities.get_units(self.design)

        if filt.shape[0] > 0:
            space_polys = shapely.unary_union(filt['geometry'].buffer(0))
            space_polys = shapely.affinity.scale(space_polys, xfact=unit_conv, yfact=unit_conv, origin=(0,0))
            space_polys = ShapelyEx.fuse_polygons_threshold(space_polys, fuse_threshold)
            if isinstance(space_polys, shapely.geometry.multipolygon.MultiPolygon):
                space_polys = [x for x in space_polys.geoms]
            else:
                space_polys = [space_polys] #i.e. it's just a lonely Polygon object...
            return space_polys
        else:
            return []

    def add_ground_plane(self, **kwargs):
        '''
        Adds metallic ground-plane from the Qiskit-Metal design object onto the surface layer of the chip simulation.

        Inputs:
            - threshold - (Optional) Defaults to -1. This is the threshold in metres, below which consecutive vertices along a given polygon are
                          combined into a single vertex. This simplification helps with meshing as COMSOL will not overdo the meshing. If this
                          argument is negative, the argument is ignored.
            - fuse_threshold - (Optional) Defaults to 1e-12. This is the minimum distance between metallic elements, below which they are considered
                               to be a single polygon and thus, the polygons are merged with the gap filled. This accounts for floating-point errors
                               that make adjacent elements fail to merge as a single element, due to infinitesimal gaps between them.
        '''

        #Assumes that all masks are simply closed without any interior holes etc... Why would you need one anyway?
        thresh = kwargs.get('threshold', -1)

        space_polys = self._get_gnd_plane_gaps(kwargs.get('fuse_threshold',1e-12))
        
        #Now create a plane sheet covering the entire chip
        pol_name = "polgndBase"
        xL = self.chip_centre[0] - self.chip_len*0.5
        xR = self.chip_centre[0] + self.chip_len*0.5
        yL = self.chip_centre[1] - self.chip_wid*0.5
        yR = self.chip_centre[1] + self.chip_wid*0.5
        self._create_poly(pol_name, [[xL,yL], [xR,yL], [xR,yR], [xL,yR]])
        if len(space_polys) > 0:
            pol_name_cuts = []
            for ind,cur_cut in enumerate(space_polys):
                pol_name_cut = f"polGND{ind}"
                pol_name_cuts.append(pol_name_cut)
                cur_poly_cut = np.array(cur_cut.exterior.coords[:])
                cur_poly_cut = self._simplify_geom(cur_poly_cut, thresh)
                #
                poly_interiors = cur_cut.interiors
                if len(poly_interiors) > 0:
                    pol_name_ints = []
                    for ind,cur_int in enumerate(poly_interiors):
                        pol_name_int = f"{pol_name_cut}S{ind}"
                        pol_name_ints.append(pol_name_int)
                        #Convert coordinates into metres...
                        cur_poly_int = np.array(cur_int.coords[:])
                        cur_poly_int = self._simplify_geom(cur_poly_int, thresh)
                        sel_x, sel_y, sel_r = self._create_poly(pol_name_int, cur_poly_int[:-1])
                    #Subtract interiors from the polygon
                    pol_ext_name = f"{pol_name_cut}E"
                    sel_x, sel_y, sel_r = self._create_poly(pol_ext_name, cur_poly_cut[:-1])
                    diff_name = pol_name_cut
                    self._model.java.component("comp1").geom("geom1").feature("wp1").geom().create(diff_name, "Difference")
                    self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(diff_name).selection("input").set(pol_ext_name)
                    self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(diff_name).selection("input2").set(*pol_name_ints)
                else:
                    sel_x, sel_y, sel_r = self._create_poly(pol_name_cut, cur_poly_cut[:-1])
            #Subtract cuts from the polygon
            diff_name = "difGND"
            self._model.java.component("comp1").geom("geom1").feature("wp1").geom().create(diff_name, "Difference")
            self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(diff_name).selection("input").set(pol_name)
            self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(diff_name).selection("input2").set(*pol_name_cuts)
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

    def fine_mesh_in_rectangle(self, x1, y1, x2, y2, minElementSize=1e-7, maxElementSize=5e-6):
        assert x2>x1 and y2>y1, "Ensure (x1,y1) is the bottom-left corner while (x2,y2) is the top-right corner of the rectangle."
        ind = len(self._fine_mesh)
        rect_name = f"mesh_rect{ind}"
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().create(rect_name, "Rectangle")
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(rect_name).set("size", jtypes.JArray(jtypes.JDouble)([ x2-x1, y2-y1 ]) )
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(rect_name).set("pos", jtypes.JArray(jtypes.JDouble)([ x1, y1 ]))
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().run(rect_name) #Makes selection easier...
        #
        sel_mesh_name = f"selMesh{ind}"
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().create(sel_mesh_name, "ExplicitSelection")
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(sel_mesh_name).selection("selection").set(rect_name, 1)
        #
        leLine = shapely.LineString([[x1, y1], [x2, y1], [x2, y2], [x1, y2], [x1, y1]])
        self._fine_mesh += [{'type':'all', 'poly':leLine, 'sel_rect':"wp1_"+sel_mesh_name, 'minElem':minElementSize, 'maxElem':maxElementSize}]

    def fine_mesh_conductors_in_rectangle(self, x1, y1, x2, y2, minElementSize=1e-7, maxElementSize=5e-6):
        assert x2>x1 and y2>y1, "Ensure (x1,y1) is the bottom-left corner while (x2,y2) is the top-right corner of the rectangle."
        ind = len(self._fine_mesh)
        rect_name = f"mesh_rect{ind}"
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().create(rect_name, "Rectangle")
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(rect_name).set("size", jtypes.JArray(jtypes.JDouble)([ x2-x1, y2-y1 ]) )
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(rect_name).set("pos", jtypes.JArray(jtypes.JDouble)([ x1, y1 ]))
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().run(rect_name) #Makes selection easier...
        #
        sel_mesh_name = f"selMesh{ind}"
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().create(sel_mesh_name, "ExplicitSelection")
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(sel_mesh_name).selection("selection").set(rect_name, 1)
        #
        leLine = shapely.LineString([[x1, y1], [x2, y1], [x2, y2], [x1, y2], [x1, y1]])
        self._fine_mesh += [{'type':'conds', 'poly':leLine, 'sel_rect':"wp1_"+sel_mesh_name, 'minElem':minElementSize, 'maxElem':maxElementSize}]

    def fine_mesh_in_poly(self, poly, minElementSize=1e-7, maxElementSize=5e-6, qmUnits=True, buffer_size=0):
        if buffer_size > 0 and qmUnits:
            buffer_size *= QUtilities.get_units(self.design)
        leCoords = np.array(poly.buffer(buffer_size, join_style=2, cap_style=3).exterior.coords[:])
        if qmUnits:
            leCoords *= QUtilities.get_units(self.design)
        ind = len(self._fine_mesh)
        pol_name = f"mesh_poly{ind}"
        sel_x, sel_y, sel_r = self._create_poly(pol_name, leCoords)
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().run(pol_name) #Makes selection easier...
        #
        sel_mesh_name = f"selMesh{ind}"
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().create(sel_mesh_name, "ExplicitSelection")
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(sel_mesh_name).selection("selection").set(pol_name, 1)
        #
        self._fine_mesh += [{'type':'all', 'poly':shapely.LineString(leCoords), 'sel_rect':"wp1_"+sel_mesh_name, 'minElem':minElementSize, 'maxElem':maxElementSize}]

    def fine_mesh_around_comp_boundaries(self, list_comp_names, minElementSize=1e-7, maxElementSize=5e-6, **kwargs):
        kwargs['restrict_rect'] = self.restrict_rect
        kwargs['resolution'] = self._resolution
        list_polys = QUtilities.get_perimetric_polygons(self.design, list_comp_names, **kwargs)
        self.fine_mesh_in_polys(list_polys, minElementSize, maxElementSize, qmUnits=False, buffer_size=kwargs.get('buffer_size',0))

    def fine_mesh_in_polys(self, list_polys, minElementSize=1e-7, maxElementSize=5e-6, qmUnits=True, buffer_size=0):
        ind = len(self._fine_mesh)
        pol_names = []
        leLines = []
        for m, poly in enumerate(list_polys):
            if buffer_size > 0 and qmUnits:
                buffer_size *= QUtilities.get_units(self.design)
            leCoords = np.array(poly.buffer(buffer_size, join_style=2, cap_style=3).exterior.coords[:])
            if leCoords.size == 0:
                continue
            if qmUnits:
                leCoords *= QUtilities.get_units(self.design)
            leLines += [shapely.LineString(leCoords)]
            pol_name = f"mesh_poly{ind}_{m}"
            sel_x, sel_y, sel_r = self._create_poly(pol_name, leCoords)
            self._model.java.component("comp1").geom("geom1").feature("wp1").geom().run(pol_name) #Makes selection easier...
            pol_names += [pol_name]
        sel_mesh_name = f"selMesh{ind}"
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().create(sel_mesh_name, "ExplicitSelection")
        for pol_name in pol_names:
            self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(sel_mesh_name).selection("selection").set(pol_name, 1)
        #
        self._fine_mesh += [{'type':'all', 'poly':shapely.MultiLineString(leLines), 'sel_rect':"wp1_"+sel_mesh_name, 'minElem':minElementSize, 'maxElem':maxElementSize}]

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

    def _create_boundary_selection_sphere(self, radius, pos_x,pos_y,pos_z=0.0, is_edge = False):
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
        self._model.java.selection(sel_name).set('entitydim', '2' if not is_edge else '1')
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



    def _create_material(self, name, rel_permit, rel_permea, conductivity, loss_tangent=0):
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
        if loss_tangent > 0:
            self._model.java.component("comp1").material(name).propertyGroup("def").set("relpermittivity", f"{rel_permit}*(1-{loss_tangent}i)")
        else:
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

    def _dset_exists(self, dset_name):
        return dset_name in [x.tag() for x in self._model.java.result().dataset()]

    def _numeval_exists(self, numeval_name):
        return numeval_name in [x.tag() for x in self._model.java.result().numerical()]

    def build_geom_mater_elec_mesh(self, mesh_structure='Normal', skip_meshing=False, substrate_material=None, substrate_permittivity=11.45, **kwargs):    #mesh_structure can be 'Fine'
        '''
        Builds geometry, sets up materials, sets up electromagnetic parameters/ports and builds the mesh.
        '''

        #Prepare the dielectric selection...
        pol_name = "polSurface"
        x1 = self.chip_centre[0] - self.chip_len*0.5
        x2 = self.chip_centre[0] + self.chip_len*0.5
        y1 = self.chip_centre[1] - self.chip_wid*0.5
        y2 = self.chip_centre[1] + self.chip_wid*0.5
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().create(pol_name, "Rectangle")
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(pol_name).set("size", jtypes.JArray(jtypes.JDouble)([ x2-x1, y2-y1 ]) )
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(pol_name).set("pos", jtypes.JArray(jtypes.JDouble)([ x1, y1 ]))
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().run(pol_name) #Makes selection easier...
        #
        sel_surface_name = "selSurface"
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().create(sel_surface_name, "ExplicitSelection")
        self._model.java.component("comp1").geom("geom1").feature("wp1").geom().feature(sel_surface_name).selection("selection").set(pol_name, 1)

        #Build geometry
        self._model.java.component("comp1").geom("geom1").run()
        #
        select_3D_name = "condAll"
        self._model.java.component("comp1").geom("geom1").create(select_3D_name, "UnionSelection")
        self._model.java.component("comp1").geom("geom1").feature(select_3D_name).label(select_3D_name)
        self._model.java.component("comp1").geom("geom1").feature(select_3D_name).set("entitydim", jtypes.JInt(2))
        self._model.java.component("comp1").geom("geom1").feature(select_3D_name).set("input", jtypes.JArray(jtypes.JString)([x[1] for x in self._conds]))
        self._model.java.component("comp1").geom("geom1").run("condAll")
        #
        diff_sel_name = "dielectricSurface"
        self._model.java.component("comp1").geom("geom1").create(diff_sel_name, "DifferenceSelection")
        self._model.java.component("comp1").geom("geom1").feature(diff_sel_name).set("entitydim", jtypes.JInt(2))
        self._model.java.component("comp1").geom("geom1").feature(diff_sel_name).set("add", jtypes.JArray(jtypes.JString)([f"wp1_{sel_surface_name}"]))
        self._model.java.component("comp1").geom("geom1").feature(diff_sel_name).set("subtract", jtypes.JArray(jtypes.JString)(["condAll"]))

        #Create materials
        if isinstance(substrate_material, Material):
            Er = substrate_material.permittivity
            loss_tangent = substrate_material.loss_tangent
        else:
            Er = np.real(substrate_permittivity)
            loss_tangent = -np.imag(substrate_permittivity)/Er
        #
        self._create_material('Vacuum', 1.0, 1.0, 0.0)
        self._model.java.component("comp1").material('Vacuum').selection().all()
        self._create_material('Substrate', Er, 1.0, 0.0, loss_tangent)     #Don't use 4.35e-4 for conductivity as it can cause issues with convergence in inductance simulations...
        self._model.java.component("comp1").material('Substrate').selection().set(self._model.java.selection('geom1_blk_chip_dom').entities(3))
        #Set the conductive boundaries to Aluminium (ignored by s-parameter and capacitance-matrix simulations due to PEC, but inductor simulations require conductivity)
        # cond_bounds = [self._get_selection_boundaries(x)[0] for x in self._conds]
        self._create_material('Metal', 0.0, 1.0, 3.77e7)
        self._model.java.component("comp1").material("Metal").selection().geom("geom1", 2)
        self._model.java.component("comp1").material("Metal").selection().named("geom1_condAll")

        for cur_sim in self._cur_sims:
            cur_sim._run_premesh()

        #Create mesh
        for m, fmData in enumerate(self._fine_mesh):
            sel_mesh_name = f'selFineMesh{m}'
            self._model.java.component("comp1").geom("geom1").create(sel_mesh_name, "IntersectionSelection")
            self._model.java.component("comp1").geom("geom1").feature(sel_mesh_name).set("entitydim", jtypes.JInt(2))
            sel_list = [fmData['sel_rect']]
            if fmData['type'] == 'conds':
                sel_list += ["condAll"]
            self._model.java.component("comp1").geom("geom1").feature(sel_mesh_name).set("input", jtypes.JArray(jtypes.JString)(sel_list))
            self._model.java.component("comp1").geom("geom1").feature(sel_mesh_name).label(f"selFineMesh{m}")
            self._model.java.component("comp1").geom("geom1").run(sel_mesh_name)

            mesh_name = "ftri" + str(m)
            self._model.java.component("comp1").mesh("mesh1").create(mesh_name, "FreeTri")
            self._model.java.component("comp1").mesh("mesh1").feature(mesh_name).selection().named("geom1_"+sel_mesh_name)
            self._model.java.component("comp1").mesh("mesh1").feature(mesh_name).create("size1", "Size")
            self._model.java.component("comp1").mesh("mesh1").feature(mesh_name).feature("size1").set("hauto", jtypes.JInt(1)) #Setting to extremely fine and then changing it - does this even do anything?
            self._model.java.component("comp1").mesh("mesh1").feature(mesh_name).feature("size1").set('custom', jtypes.JBoolean(True))
            self._model.java.component("comp1").mesh("mesh1").feature(mesh_name).feature("size1").set('hmaxactive', jtypes.JBoolean(True))
            self._model.java.component("comp1").mesh("mesh1").feature(mesh_name).feature("size1").set('hminactive', jtypes.JBoolean(True))
            self._model.java.component("comp1").mesh("mesh1").feature(mesh_name).feature("size1").set('hmin', jtypes.JDouble(fmData['minElem']))
            self._model.java.component("comp1").mesh("mesh1").feature(mesh_name).feature("size1").set('hmax', jtypes.JDouble(fmData['maxElem']))
        # self._model.java.component("comp1").mesh("mesh1").feature("size").set("custom", jtypes.JBoolean(True))
        # self._model.java.component("comp1").mesh("mesh1").feature("size").set("hmax", jtypes.JDouble(8e-4))
        # self._model.java.component("comp1").mesh("mesh1").feature("size").set("hmin", jtypes.JDouble(10e-6))
        self._model.java.component("comp1").mesh("mesh1").create("ftet10", "FreeTet")
        self._model.java.component("comp1").mesh("mesh1").feature("ftet10").create("size1", "Size")
        mesh_auto = {'Extremely fine' : 1, 'Extra fine' : 2, 'Finer' : 3, 'Fine' : 4, 'Normal' : 5, 'Coarse' : 6, 'Coarser' : 7, 'Extra coarse' : 8, 'Extremely coarse' : 9, 'Custom' : 1}
        assert mesh_structure in mesh_auto, f"Predefined mesh type \'{mesh_structure}\' is not supported/recognised."
        self._model.java.component("comp1").mesh("mesh1").feature("ftet10").feature("size1").set("hauto", jtypes.JInt(mesh_auto[mesh_structure]))
        if mesh_structure == 'Custom':
            self._model.java.component("comp1").mesh("mesh1").feature("ftet10").feature("size1").set('custom', jtypes.JBoolean(True))
            if 'minElementSize' in kwargs:
                self._model.java.component("comp1").mesh("mesh1").feature("ftet10").feature("size1").set('hminactive', jtypes.JBoolean(True))
                self._model.java.component("comp1").mesh("mesh1").feature("ftet10").feature("size1").set('hmin', jtypes.JDouble(kwargs['minElementSize']))
            if 'maxElementSize' in kwargs:
                self._model.java.component("comp1").mesh("mesh1").feature("ftet10").feature("size1").set('hmaxactive', jtypes.JBoolean(True))
                self._model.java.component("comp1").mesh("mesh1").feature("ftet10").feature("size1").set('hmax', jtypes.JDouble(kwargs['maxElementSize']))
            kwargs
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


    def plot(self, plot_chip_boundaries = True, plot_fine_mesh_bounds = True):
        names = [x[1] for x in self._cond_polys]
        geoms = [x[0] for x in self._cond_polys]

        if plot_fine_mesh_bounds:
            for m, cur_fine_mesh in enumerate(self._fine_mesh):
                if cur_fine_mesh['type'] == 'conds':
                    cur_name = f'FineMesh{m}(Conds)'
                elif cur_fine_mesh['type'] == 'all':
                    cur_name = f'FineMesh{m}(Region)'
                else:
                    continue

                geoms += [cur_fine_mesh['poly']]
                names += [cur_name]

        gdf = gpd.GeoDataFrame({'names':names}, geometry=geoms)
        fig, ax = plt.subplots(1)

        if plot_chip_boundaries:
            minX = (self.chip_centre[0] - self.chip_len*0.5)
            maxX = (self.chip_centre[0] + self.chip_len*0.5)
            minY = (self.chip_centre[1] - self.chip_wid*0.5)
            maxY = (self.chip_centre[1] + self.chip_wid*0.5)
            data = np.array([(minX,minY),(maxX,minY),(maxX,maxY),(minX,maxY),(minX,minY)])
            ax.plot(data[:,0], data[:,1], 'k')

        gdf.plot(ax = ax, column='names', cmap='jet', alpha=0.2, categorical=True, legend=True, aspect='equal')

    def save(self, file_name):
        self._model.java.component("comp1").view("view1").set("transparency", jtypes.JBoolean(True))    #Just a nice helpful feature to stop one extra button click...
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
