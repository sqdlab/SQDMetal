from qiskit_metal.qgeometries.qgeometries_handler import QGeometryTables
from qiskit_metal.renderers.renderer_mpl.mpl_renderer import QMplRenderer
from SQDMetal.Utilities.QiskitShapelyRenderer import QiskitShapelyRenderer
from SQDMetal.Utilities.ShapelyEx import ShapelyEx
from qiskit_metal.qlibrary.terminations.launchpad_wb import LaunchpadWirebond
from SQDMetal.Utilities.QUtilities import QUtilities
import gmsh
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import colormaps
import shapely
import qiskit_metal 


class GMSH_Geometry_Builder:

    def __init__(self, design, fillet_resolution, gmsh_verbosity=1):
        #gmsh_verbosity: 0: silent except for fatal errors, 1: +errors, 2: +warnings, 3: +direct, 4: +information, 5: +status, 99: +debug
        
        self.design = design
        self.fillet_resolution = fillet_resolution

        unit_conv = 1e-3

        #Get dimensions of chip base and convert to design units in 'mm'
        self.center_x = QUtilities.parse_value_length(design.chips['main']['size']['center_x']) / unit_conv
        self.center_y = QUtilities.parse_value_length(design.chips['main']['size']['center_y']) / unit_conv
        self.center_z = QUtilities.parse_value_length(design.chips['main']['size']['center_z']) / unit_conv
        self.size_x = QUtilities.parse_value_length(design.chips['main']['size']['size_x']) / unit_conv
        self.size_y = QUtilities.parse_value_length(design.chips['main']['size']['size_y']) / unit_conv
        self.size_z = QUtilities.parse_value_length(design.chips['main']['size']['size_z']) / unit_conv

        #Initialize the GMSH API and name the model
        gmsh.model.remove()
        gmsh.finalize()
        gmsh.initialize()
        gmsh.model.add('qiskit_to_gmsh')
        self._verbosity = gmsh_verbosity
        gmsh.option.setNumber('General.Verbosity', self._verbosity)
        gmsh.option.setNumber("General.Terminal", self._verbosity)

      
    def construct_geometry_in_GMSH(self, metallic_layers, ground_plane, sim_constructs, fine_meshes, fuse_threshold, **kwargs):
        '''This function takes the existing geometry in the Qiskit Metal design and constructs them in GMSH.
        
        Args:
            None.

        Returns:
            
        '''

        #Do pre-processing in shapely to get metallic elements and dielectric cutouts ready to build in GMSH.
        #Note: metals list contains ground plane. Dielectric gaps are the difference between the dielectric cutout
        #and the metals
        metals, dielectric_gaps = self._process_qiskit_geometries(metallic_layers, ground_plane, fuse_threshold, **kwargs)

        #Plot the shapely metals for user to see device
        # geoms = metals# + ground_plane
        # metal_names = [str(i) for i,_ in enumerate(geoms)]
        # gdf = gpd.GeoDataFrame({'names':metal_names}, geometry=geoms)
        # fig, ax = plt.subplots()
        # gdf.plot(ax = ax, column='names', cmap='tab10', alpha=0.75, categorical=True, legend=True)
        # plt.show()

        #Draw shapely metal and dielectric gap polygons into GMSH. If the polygon has an interior, the interior sections are subtracted from
        #the exterior boundary
        metal_list = self._create_gmsh_geometry_from_shapely_polygons(metals)                           #create all metal traces
        dielectric_gap_list = self._create_gmsh_geometry_from_shapely_polygons(dielectric_gaps)         #create all the gaps between the metal traces and the ground plane 

        #Create the chip base (dielectric substrate) in Gmsh
        chip_base = self._create_chip_base()
        
        #Create the airbox which houses the chip base
        bottom_grounded = True
        air_box = self._draw_air_box(bottom_grounded)

        fragment_list = []


        # #create JJ's as lumped ports, if they exist
        # junctions = []
        # if not self.design.qgeometry.tables['junction'].empty:
        #     junctions, jj_inductance = self._create_JJ_lumped_ports()
        
        # #create junction list for boolean operations in gmsh
        # junction_list = [(2,x) for x in junctions]

        # #if there are junctions in the design, create the physical group
        # jj_dict = {}
        # if junctions:
        #     for i,junction in enumerate(junctions):
        #         jj_name = 'jj_' + str(i)
        #         jj_physical_group = gmsh.model.addPhysicalGroup(2, [junction], name = jj_name)
        #         jj_dict[jj_name] = (jj_physical_group, jj_inductance)

        #Process simulation constructs - e.g. ports
        lePortPolys = {}
        aux_list = []
        metal_cut_list_for_ports = [(x[0],x[1]) for x in metal_list]
        for cur_sim_poly in sim_constructs:
            sim_poly = self._draw_polygon_in_GMSH_from_coords(np.array(cur_sim_poly[1])[:-1]*1e3) #i.e. ignore closed loop part as function closes it automatically and use mm
            sim_poly,leMap = gmsh.model.occ.cut([(2,sim_poly)], metal_cut_list_for_ports, removeObject=True, removeTool=False)
            gmsh.model.occ.synchronize()
            gmsh.model.geo.synchronize()
            #
            sim_poly_list = [(2,x[1]) for x in sim_poly]
            fragment_list += sim_poly_list
            aux_list += sim_poly_list
            lePortPolys[cur_sim_poly[0]] = [x[1] for x in sim_poly]

        #list for metals for capacitance simulation
        metal_cap_physical_group = []
        metal_cap_names = []
        
        #Dielectric gap list needs to be updated to account for the presence of ports/junctions etc...
        if len(aux_list) > 0:
            cut_list = aux_list #ports_list + junction_list
            dielectric_gap_list_for_cut = [(x[0],x[1]) for x in dielectric_gap_list]
            dielectric_gap_list,lemap = gmsh.model.occ.cut(dielectric_gap_list_for_cut, cut_list, removeObject=True, removeTool=False)
        fragment_list = metal_list + dielectric_gap_list + fragment_list

        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()

        fragment_list = [(x[0],x[1]) for x in fragment_list]
        #See here: https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/tutorials/python/t16.py for details on the mapping...
        chip, chip_map = gmsh.model.occ.fragment([(3,chip_base)], fragment_list, removeObject=True, removeTool=True)
        surface_remappings = {}
        for e in zip([(3,chip_base)] + fragment_list, chip_map):
            # print("parent " + str(e[0]) + " -> child " + str(e[1]))
            surface_remappings[e[0]] = e[1]

        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()

        #Create a physical group for each metal in the design. This is important for capacitance simulations
        for i,metal in enumerate(metal_list):
            metal_name = 'metal_' + str(i)
            cur_metal_surfaces = surface_remappings[metal]
            metal_physical_group = gmsh.model.addPhysicalGroup(2, [x[1] for x in cur_metal_surfaces], name = metal_name)
            metal_cap_physical_group.append(metal_physical_group)
            metal_cap_names.append(metal_name)

        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()

        #Create a physical group for each port/junction simulation construct...
        cur_port_names = [x for x in lePortPolys]
        for cur_port in cur_port_names:
            cur_port_surfaces = [(2,x) for x in lePortPolys[cur_port]]
            cur_port_surfaces_final = []
            for cur_port_surface in cur_port_surfaces:
                cur_port_surfaces_final += surface_remappings[cur_port_surface]
            lePortPolys[cur_port] = gmsh.model.addPhysicalGroup(2, [x[1] for x in cur_port_surfaces_final], name = cur_port)

        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()

        #Add physical group for dielectric gap list
        leDielectricGaps = gmsh.model.addPhysicalGroup(2, [x[1] for x in dielectric_gap_list], name = 'dielectric_gaps')

        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()

        #@DavidSommers - basically the mesh elements must be conformal. Thus, as the fragments have been done with
        #surface metals onto the substrate, this needs to be done onto the air-box as well...
        #
        #Add in physical groups for each part of the device
        leDielectric = gmsh.model.addPhysicalGroup(3, [x[1] for x in chip if x[0] == 3], name = 'dielectric_substrate') #Only taking 3D elements from chip...
        #
        #Fragment the dielectric volume with the airbox
        chip_and_air_box, chip_and_air_box_map = gmsh.model.occ.fragment(chip, [(3, air_box)], removeObject=True, removeTool=True)
        
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()
        
        leAirBox = gmsh.model.addPhysicalGroup(3, [air_box], name = 'air_box')
        
        #Get all gmsh elements which comprise the airbox and chip base.
        #These elements can be used to determine the far-field conditions. 
        _,air_box_surfaces = gmsh.model.getAdjacencies(3,air_box)
        _,dielectric_box_surfaces = gmsh.model.getAdjacencies(3,chip_base)

        #Define far-field boundary condition based on whether the dielectric substrate has a back-side ground plane.
        #For the case where the substrate is grounded there is an extra element to add.
        if bottom_grounded == True:
            far_field_surfaces = air_box_surfaces[:6].tolist() + [air_box_surfaces[0]-1]
        else:
            far_field_surfaces = air_box_surfaces[:6]

        leFarField = gmsh.model.addPhysicalGroup(2, far_field_surfaces, name = 'far_field')
        
        #synchronise
        gmsh.model.geo.synchronize()
        gmsh.model.occ.synchronize()

        fine_mesh_elems = []
        for cur_fine_mesh in fine_meshes:
            if cur_fine_mesh['type'] == 'box':
                cur_mesh_attrb = {}
                x1, x2 = cur_fine_mesh['x_bnds']
                y1, y2 = cur_fine_mesh['y_bnds']
                #Get it in mm...
                cur_mesh_attrb['region'] = self._create_gmsh_geometry_from_shapely_polygons([ShapelyEx.rectangle(x1*1e3, y1*1e3, x2*1e3, y2*1e3)])
                cur_mesh_attrb['mesh_min'] = cur_fine_mesh['min_size'] * 1e3
                cur_mesh_attrb['mesh_max'] = cur_fine_mesh['max_size'] * 1e3
                cur_mesh_attrb['taper_dist_min'] = cur_fine_mesh['taper_dist_min'] * 1e3
                cur_mesh_attrb['taper_dist_max'] = cur_fine_mesh['taper_dist_max'] * 1e3
                #
                fine_mesh_elems.append(cur_mesh_attrb)
            elif cur_fine_mesh['type'] == 'comp_bounds':
                cur_mesh_attrb = {}
                comp_outlines = QUtilities.get_perimetric_polygons(self.design, cur_fine_mesh['list_comp_names'], fuse_threshold=fuse_threshold, resolution=self.fillet_resolution, unit_conv=1, metals_only=cur_fine_mesh['metals_only'])    #Get it in mm...
                cur_mesh_attrb['region'] = self._create_gmsh_geometry_from_shapely_polygons(comp_outlines)
                cur_mesh_attrb['mesh_min'] = cur_fine_mesh['min_size'] * 1e3
                cur_mesh_attrb['mesh_max'] = cur_fine_mesh['max_size'] * 1e3
                cur_mesh_attrb['taper_dist_min'] = cur_fine_mesh['taper_dist_min'] * 1e3
                cur_mesh_attrb['taper_dist_max'] = cur_fine_mesh['taper_dist_max'] * 1e3
                #
                fine_mesh_elems.append(cur_mesh_attrb)
            elif cur_fine_mesh['type'] == 'path':
                cur_mesh_attrb = {}
                #Get it in mm...
                lePoints = [gmsh.model.geo.addPoint(x[0]*1e3,x[1]*1e3, 0) for x in cur_fine_mesh['path']]
                leLines = [gmsh.model.geo.addLine(lePoints[m-1], lePoints[m]) for m in range(1, len(lePoints))]
                gmsh.model.occ.synchronize()
                gmsh.model.geo.synchronize()
                #
                cur_mesh_attrb['path'] = leLines
                cur_mesh_attrb['mesh_min'] = cur_fine_mesh['min_size'] * 1e3
                cur_mesh_attrb['mesh_max'] = cur_fine_mesh['max_size'] * 1e3
                cur_mesh_attrb['taper_dist_min'] = cur_fine_mesh['taper_dist_min'] * 1e3
                cur_mesh_attrb['taper_dist_max'] = cur_fine_mesh['taper_dist_max'] * 1e3
                #
                fine_mesh_elems.append(cur_mesh_attrb)


        ret_dict = {
            'air_box'   : self._conv_to_lists(leAirBox),
            'metals'    : self._conv_to_lists(metal_cap_physical_group),
            'metalsShapely'    : metals,
            'far_field' : self._conv_to_lists(leFarField),
            'ports'     : lePortPolys,
            'dielectric': self._conv_to_lists(leDielectric),
            'dielectric_gaps': self._conv_to_lists(leDielectricGaps),
            'fine_mesh_elems': fine_mesh_elems
        }
        return ret_dict

    def _conv_to_lists(self, data):
        if isinstance(data, list):
            return data
        else:
            return [data]

    def _process_qiskit_geometries(self, metallic_layers, ground_plane, fuse_threshold, **kwargs):
        fuse_threshold /= 1e-3  #Convert to mm...
        threshold = kwargs.get('threshold', 1e-9) / 1e-3  #Convert to mm...

        #Remove any previous model and add new gmsh model
        gmsh.model.remove()
        gmsh.finalize()
        gmsh.initialize()
        gmsh.model.add('qiskit_to_gmsh')
        gmsh.option.setNumber('General.Verbosity', self._verbosity)
        gmsh.option.setNumber("General.Terminal", self._verbosity)

        unit_conv = 1   #TODO: Because it is stuck in mm?!

        #Handle the ground plane...
        qmpl = QiskitShapelyRenderer(None, self.design, None)
        gsdf = qmpl.get_net_coordinates(self.fillet_resolution)
        filt = gsdf.loc[gsdf['subtract'] == True]
        #TODO: It should only take stuff from the metallic layers specified when calculating spaces???
        space_polys = ShapelyEx.fuse_polygons_threshold(filt['geometry'].buffer(0), fuse_threshold)
        space_polys = ShapelyEx.shapely_to_list(space_polys)
        space_polys = [shapely.affinity.scale(x, xfact=unit_conv, yfact=unit_conv, origin=(0,0)) for x in space_polys]
        #
        metal_surface = shapely.geometry.box(self.center_x - 0.5*self.size_x, self.center_y - 0.5*self.size_y,
                                             self.center_x + 0.5*self.size_x, self.center_y + 0.5*self.size_y)
        if not ground_plane['omit']:
            ground_plane_poly = shapely.difference(metal_surface, shapely.geometry.multipolygon.MultiPolygon(space_polys))
        #
        #Gather all polygons into contiguous groups
        metal_polys = []
        if not ground_plane['omit'] and not ground_plane_poly.is_empty:
            metal_polys.append(ground_plane_poly)
        for cur_layer in metallic_layers:
            if cur_layer['type'] == 'design_layer':
                cur_layer['unit_conv'] = unit_conv
                cur_layer['resolution'] = self.fillet_resolution
                cur_layer['threshold'] = threshold
                metal_polys_all, metal_sel_ids = QUtilities.get_metals_in_layer(self.design, **cur_layer)
                metal_polys += metal_polys_all
            elif cur_layer['type'] == 'Uclip':
                if cur_layer['clip_type'] == 'inplaneLauncher':
                    Uclip = QUtilities.get_RFport_CPW_groundU_Launcher_inplane(self.design, cur_layer['qObjName'], cur_layer['thickness_side'], cur_layer['thickness_back'], cur_layer['separation_gap'], cur_layer['unit_conv_extra']*1e3) #TODO: Stuck in mm?!
                elif cur_layer['clip_type'] == 'inplaneRoute':
                    Uclip = QUtilities.get_RFport_CPW_groundU_Route_inplane(self.design, cur_layer['route_name'], cur_layer['pin_name'], cur_layer['thickness_side'], cur_layer['thickness_back'], cur_layer['separation_gap'], cur_layer['unit_conv_extra']*1e3) #TODO: Stuck in mm?!
                metal_polys.append(shapely.Polygon(Uclip))
        #Try to fuse any contiguous polygons...
        metal_polys = ShapelyEx.fuse_polygons_threshold(metal_polys)
        #If there are any MultiPolygons, convert them into normal polygons...
        new_polys = ShapelyEx.shapely_to_list(metal_polys)

        substrate = ShapelyEx.rectangle(self.center_x-self.size_x/2, self.center_y-self.size_y/2, self.center_x+self.size_x/2, self.center_y+self.size_y/2)
        dielectric_gaps = shapely.difference(substrate, metal_polys)
        dielectric_gaps = ShapelyEx.shapely_to_list(dielectric_gaps)

        return new_polys, dielectric_gaps

    def _create_gmsh_geometry_from_shapely_polygons(self, polygons):
        polygons_list = []
        for m,poly in enumerate(polygons):
            poly_simplified = poly.simplify(1e-9) #this removes points that are spaced too closely together - TODO: investigate why this is still required?
            cur_coords = poly_simplified.exterior.coords[:-1]
            if len(cur_coords) == 0 or poly_simplified.area < 1e-18:
                continue
            gmsh_exterior = self._draw_polygon_in_GMSH_from_coords(cur_coords) #remove last coord from list because it's repeated
            if len(poly_simplified.interiors) >= 1:
                interiors = []
                for _,interior in enumerate(poly_simplified.interiors):
                    gmsh_interior = self._draw_polygon_in_GMSH_from_coords(interior.coords[:-1]) #remove last coord from list because it's repeated
                    interiors.append((2,gmsh_interior))
                gmsh_surface, gmsh_surface_map = gmsh.model.occ.cut([(2,gmsh_exterior)], 
                                                                    interiors,
                                                                    removeObject=True, removeTool=True)
                if len(gmsh_surface) > 0:
                    polygons_list.append((2,gmsh_surface[0][1]))
            else:
                polygons_list.append((2,gmsh_exterior))
        
        return polygons_list


    def _draw_polygon_in_GMSH_from_coords(self, coords):
        '''This function takes the coordinates of shapely polygons and then draws them into GMSH.
        
        Args:
            coords - coordinates of polygons as taken from the processed shapely objects.

        Returns:
            surface - GMSH surface ID of newly drawn polygon.
        '''
        #gets polygon coordinates depending on the type of argument fed to the function
        x_coords = []
        y_coords = []

        for _,coordinates in enumerate(coords):
            x_coords.append(coordinates[0])
            y_coords.append(coordinates[1])

        #define lists to store points and lines
        points = []
        lines = []
    
        #create 2D points in gmsh
        for i,_ in enumerate(x_coords):
            point = gmsh.model.occ.addPoint(x_coords[i], y_coords[i], 
                                            self.design.parse_value(self.design.chips['main'].size.center_z))
            points.append(point)

        #update model with the created geometry items
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()

        #draw lines from points
        for j,value in enumerate(points):
            if(j<len(points)-1):
                line = gmsh.model.occ.add_line(points[j], points[j+1])
                lines.append(line)

        line = gmsh.model.occ.add_line(points[len(points)-1], points[0])
        lines.append(line)      
        
        #create curved loop
        curve_loop = gmsh.model.occ.add_curve_loop(lines)
        
        #create_surface
        surface = gmsh.model.occ.add_plane_surface([curve_loop])
        
        #update model with the created geometry items
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()

        return surface


    def _create_chip_base(self):
        '''Creates the dielectric chip volume in GMSH from the dimensions specified in Qiskit Metal.

        Args:
            None.

        Returns:
            Tuple of integers with the first value specifying the dimension and the second value the ID of the object in GMSH.
        
        '''
        
        #define region of device for simulation


        #surface of chip base rests at z = 0
        chip_base = gmsh.model.occ.addBox(self.center_x-self.size_x/2, self.center_y-self.size_y/2, self.center_z, self.size_x, self.size_y, -np.abs(self.size_z))
        gmsh.model.occ.synchronize()

        return chip_base
    

    def _draw_air_box(self, bottom_grounded):
        '''Creates the airbox which the chip will be placed in.

        Args:
            None.

        Returns:
            Tuple of integers with the first value specifying the dimension and the second value the ID of the object in GMSH. 
        '''

        #for air box choose increase in dimensions as a percentage of the chip dimensions
        increase = 0.20 
        air_box_delta_x = increase * self.size_x
        air_box_delta_y = increase * self.size_y
        
        #bottom left corner of air box
        x_point = self.center_x-self.size_x/2-air_box_delta_x/2 
        y_point = self.center_y-self.size_y/2-air_box_delta_y/2

        #Check to place ground plane on back side by choosing statring z point
        if bottom_grounded == True: 
            z_point = self.center_z - np.abs(self.size_z)
            air_box_delta_z = 2 * np.abs(self.size_z)
        else:
            z_point = self.center_z - 2 * np.abs(self.size_z)
            air_box_delta_z = 3 * np.abs(self.size_z)

        air_box = gmsh.model.occ.addBox(x_point, y_point, z_point, self.size_x + air_box_delta_x, self.size_y + air_box_delta_y, air_box_delta_z)
        gmsh.model.occ.synchronize()

        return air_box

    def _check_port_orientation(self, vec_perp):
        thresh = 0.9999

        xDot = np.dot(vec_perp, [1,0])
        if xDot > thresh:
            return '+X', '-X'
        elif xDot < -thresh:
            return '-X', '+X'
        
        yDot = np.dot(vec_perp, [0,1])
        if yDot > thresh:
            return '+Y', '-Y'
        elif yDot < -thresh:
            return '-Y', '+Y'
        
        assert False, f"AWS Palace requires RF Lumped Ports to be aligned with the x/y axes. Here the port is pointing: {vec_perp}."


    def _create_JJ_lumped_ports(self):

        #get the dataframe with the JJ parameters
        junctions_df = self.design.qgeometry.tables['junction']

        #get number of rows/columns in data frame
        rows, cols = junctions_df.shape

        junctions = []
        jj_inductance = []
        for junction_no in range(rows):
            
            inductance = junctions_df.loc[junction_no].hfss_inductance
            inductance = float(inductance[:-2]) * 1e-9

            #coords has two points which represent the JJ as a linestring
            jj_coords = junctions_df.loc[junction_no].geometry.coords[:]
            jj_width = junctions_df.loc[junction_no].width

            #find left-most coord if JJ orientation is horizontal
            if jj_coords[0][1] == jj_coords[1][1]:
                if jj_coords[0][0] < jj_coords[1][0]:
                    left_most_coord = jj_coords[0]
                    right_most_coord = jj_coords[1]
                else:
                    left_most_coord = jj_coords[1]
                    right_most_coord = jj_coords[0]

                left_x = left_most_coord[0]
                right_x = right_most_coord[0]
                left_y = left_most_coord[1] - jj_width/2
                delta_x = right_x - left_x
            
            else:
                if jj_coords[0][1] < jj_coords[1][1]:
                    bottom_most_coord = jj_coords[0]
                    up_most_coord = jj_coords[1]
                else:
                    bottom_most_coord = jj_coords[1]
                    up_most_coord = jj_coords[0]
                
                left_x = bottom_most_coord[0] - jj_width/2
                left_y = bottom_most_coord[1]
                delta_y = up_most_coord[1] - bottom_most_coord[1]

            junction = gmsh.model.occ.addRectangle(left_x, left_y, 0, jj_width, delta_y)
            gmsh.model.occ.synchronize()

            junctions.append(junction)
            jj_inductance.append(inductance)

        return junctions, jj_inductance
    
