# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0

from SQDMetal.Utilities.ShapelyEx import ShapelyEx
from SQDMetal.Utilities.QUtilities import QUtilities
import gmsh
import numpy as np

from SQDMetal.Utilities.GeometryProcessors.GeomBase import GeomBase
from SQDMetal.Utilities.GeometryProcessors.GeomQiskitMetal import GeomQiskitMetal

class GMSH_Geometry_Builder:

    def __init__(self, geom_processor:GeomBase, fillet_resolution, gmsh_verbosity=1):
        #gmsh_verbosity: 0: silent except for fatal errors, 1: +errors, 2: +warnings, 3: +direct, 4: +information, 5: +status, 99: +debug
        
        self._geom_processor = geom_processor
        self.fillet_resolution = fillet_resolution

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
        metals, dielectric_gaps = self._geom_processor.process_layers(metallic_layers, ground_plane, fuse_threshold=fuse_threshold, fillet_resolution=self.fillet_resolution, unit_conv=1e-3, **kwargs)    #It's in mm...

        unit_conv = 1e-3

        #Get dimensions of chip base and convert to design units in 'mm'
        chip_centre = self._geom_processor.chip_centre
        self.center_x = chip_centre[0] / unit_conv
        self.center_y = chip_centre[1] / unit_conv
        self.center_z = chip_centre[2] / unit_conv
        self.size_x = self._geom_processor.chip_size_x / unit_conv
        self.size_y = self._geom_processor.chip_size_y / unit_conv
        self.size_z = self._geom_processor.chip_size_z / unit_conv

        bottom_grounded = True

        #Plot the shapely metals for user to see device
        # geoms = metals# + ground_plane
        # metal_names = [str(i) for i,_ in enumerate(geoms)]
        # gdf = gpd.GeoDataFrame({'names':metal_names}, geometry=geoms)
        # fig, ax = plt.subplots()
        # gdf.plot(ax = ax, column='names', cmap='tab10', alpha=0.75, categorical=True, legend=True)
        # plt.show()


        ####CONSTRUCT AIR-BOX AND CARVE OUT THE SUBSTRATE FROM IT####
        #
        #Create substrate and air box
        tag3D_chip_base = self._create_chip_base()
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()
       
        ####CONSTRUCT THE SURFACE ELEMENTS ON TOP OF THE SUBSTRATE####
        #
        #Draw shapely metal and dielectric gap polygons into GMSH. If the polygon has an interior, the interior sections are subtracted from
        #the exterior boundary
        metal_list = self._create_gmsh_geometry_from_shapely_polygons(metals)                           #create all metal traces
        dielectric_gap_list = self._create_gmsh_geometry_from_shapely_polygons(dielectric_gaps)         #create all the gaps between the metal traces and the ground plane 
        #
        fragment_list = []
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
        #List for metals for capacitance simulation
        metal_cap_physical_group = []
        metal_cap_names = []
        #Dielectric gap list needs to be updated to account for the presence of ports/junctions etc...
        if len(aux_list) > 0:
            cut_list = aux_list #ports_list + junction_list
            dielectric_gap_list_for_cut = [(x[0],x[1]) for x in dielectric_gap_list]
            dielectric_gap_list,lemap = gmsh.model.occ.cut(dielectric_gap_list_for_cut, cut_list, removeObject=True, removeTool=False)
        fragment_list = metal_list + dielectric_gap_list + fragment_list
        #
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()
        #
        fragment_list = [(x[0],x[1]) for x in fragment_list]
        #See here: https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/tutorials/python/t16.py for details on the mapping...
        chip, chip_map = gmsh.model.occ.fragment([(3,tag3D_chip_base)], fragment_list, removeObject=True, removeTool=True)
        surface_remappings = {}
        for e in zip([(3,tag3D_chip_base)] + fragment_list, chip_map):
            # print("parent " + str(e[0]) + " -> child " + str(e[1]))
            surface_remappings[e[0]] = e[1]
        #
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()
        #Create a physical group for each metal in the design. This is important for capacitance simulations
        for m,metal in enumerate(metal_list):
            cur_metal_surfaces = surface_remappings[metal]
            metal_cap_physical_group.append([x[1] for x in cur_metal_surfaces])
        #
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
        #
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()
        #Add physical group for dielectric gap list
        leDielectricGaps = [x[1] for x in dielectric_gap_list]
        #
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()

        ####ADD IN AIR-BOX####
        #Note that the substrate volume already encompasses all the changes made to its surface etc...
        air_box = self._draw_air_box(bottom_grounded)
        #The two pieces intersect in 3D volume and may have intersecting 2D surfaces - they must be reconciled into one (i.e. conformal meshing;
        #otherwise, Gmsh will create the same mesh twice in that intersecting 3D volume etc.)
        chip_and_air_box, chip_and_air_box_map = gmsh.model.occ.fragment([(3,tag3D_chip_base)], [(3, air_box)], removeObject=True, removeTool=True)
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()
        #NOTE: THIS ASSUMES THAT THERE ARE ONLY 2 3D VOLUMES IN PLAY UP TO THIS POINT
        #The air-box will have the substrate carved out of it. Find the 3D tags for the substrate (it'll be solitary) and the air-box...
        tag3D_chip_base = chip_and_air_box_map[0][0][1]
        tag3D_air_box = [x for x in chip_and_air_box_map[1] if x[1] != tag3D_chip_base][0][1]
        #Get the outer-region for the far-field. Just find this by checking which planes have a centre-of-mass on the extents...
        all_surfaces = gmsh.model.occ.getEntities(2)
        far_field_surfaces = []
        for dim, tag in all_surfaces:
            com = gmsh.model.occ.getCenterOfMass(dim, tag)
            if np.min(np.abs(self._extents.T - np.array(com))) < 1e-6:
                far_field_surfaces.append(tag)
        leFarField = gmsh.model.addPhysicalGroup(2, far_field_surfaces, name = 'far_field')
        #        
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
                assert isinstance(self._geom_processor, GeomQiskitMetal), "The geometry must be in Qiskit-Metal format to use comp_bounds..."
                cur_mesh_attrb = {}
                comp_outlines = QUtilities.get_perimetric_polygons(self._geom_processor.design, cur_fine_mesh['list_comp_names'], fuse_threshold=fuse_threshold, resolution=self.fillet_resolution, unit_conv=1, metals_only=cur_fine_mesh['metals_only'])    #Get it in mm...
                cur_mesh_attrb['region'] = self._create_gmsh_geometry_from_shapely_polygons(comp_outlines)
                cur_mesh_attrb['mesh_min'] = cur_fine_mesh['min_size'] * 1e3
                cur_mesh_attrb['mesh_max'] = cur_fine_mesh['max_size'] * 1e3
                cur_mesh_attrb['taper_dist_min'] = cur_fine_mesh['taper_dist_min'] * 1e3
                cur_mesh_attrb['taper_dist_max'] = cur_fine_mesh['taper_dist_max'] * 1e3
                #
                fine_mesh_elems.append(cur_mesh_attrb)
            elif cur_fine_mesh['type'] == 'arb_polys':
                cur_mesh_attrb = {}
                polys = cur_fine_mesh['polys']
                polys = shapely.affinity.scale(polys, xfact=1e3, yfact=1e3, origin=(0,0))  #Get it into mm
                polys = ShapelyEx.shapely_to_list(polys)
                cur_mesh_attrb['region'] = self._create_gmsh_geometry_from_shapely_polygons(polys)
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

        lemetals_physgrps = []
        air_box_tags = [tag3D_air_box]
        #Process the metals into groups (extrude if required)
        full_3D_params = kwargs.get('full_3D_params', None)
        if full_3D_params is None or full_3D_params['metal_thickness'] == 0:
            for m,cur_metal_entities in enumerate(metal_cap_physical_group):
                metal_name = 'metal_' + str(m)
                metal_physical_group = gmsh.model.addPhysicalGroup(2, cur_metal_entities, name = metal_name)
                lemetals_physgrps.append(metal_physical_group)
                metal_cap_names.append(metal_name)  #TODO: Look into why is this even stored?
        else:
            metal_thickness = full_3D_params['metal_thickness'] * 1e3

            new_extrusion_metals = []   #Holds (dim,tag) pairs for 3D volumes
            for cur_metal_entities in metal_cap_physical_group:
                ex_metals = gmsh.model.occ.extrude([(2,x) for x in cur_metal_entities], 0, 0, metal_thickness)
                gmsh.model.occ.synchronize()
                gmsh.model.geo.synchronize()
                new_extrusion_metals += [x for x in ex_metals if x[0] == 3]

            gmsh.model.occ.synchronize()
            gmsh.model.geo.synchronize()
            chip_and_air_box, chip_and_air_box_map = gmsh.model.occ.fragment([(3,tag3D_air_box)], new_extrusion_metals, removeObject=True, removeTool=True)
            gmsh.model.occ.synchronize()
            gmsh.model.geo.synchronize()

            for m in range(len(new_extrusion_metals)):
                #Collate all 3D parts inside this mapping...
                tags_3D = [cur_tag[1] for cur_tag in chip_and_air_box_map[m+1] if cur_tag[0]==3]
                cur_metal_surfaces = []
                for cur_tag_3D in tags_3D:
                    #Note that here getAdjacencies returns 4D and 2D entities (empty list for former)
                    cur_metal_surfaces += gmsh.model.getAdjacencies(3, cur_tag_3D)[1].tolist()

                metal_name = 'metal_' + str(m)
                metal_physical_group = gmsh.model.addPhysicalGroup(2, cur_metal_surfaces, name = metal_name)
                lemetals_physgrps.append(metal_physical_group)
                metal_cap_names.append(metal_name)  #TODO: Look into why is this even stored?

            #Figure out the 3D tag for the air-box...
            new_extrude_tags = []
            for x in chip_and_air_box_map[1:]:
                new_extrude_tags += [y[1] for y in x if y[0]==3]
            air_box_tags = [x[1] for x in chip_and_air_box_map[0] if x[1] not in new_extrude_tags and x[0]==3]
        #Process trenching into dielectric
        if full_3D_params is not None and full_3D_params['substrate_trenching'] > 0:
            trench_depth = full_3D_params['substrate_trenching'] * 1e3

            trenches = gmsh.model.occ.extrude(dielectric_gap_list, 0, 0, -trench_depth)
            gmsh.model.occ.synchronize()
            gmsh.model.geo.synchronize()
            new_trench_extrusions = [x for x in trenches if x[0] == 3]

            gmsh.model.occ.synchronize()
            gmsh.model.geo.synchronize()
            chip_and_trenches, chip_and_trenches_map = gmsh.model.occ.fragment([(3,tag3D_chip_base)], new_trench_extrusions, removeObject=True, removeTool=True)
            gmsh.model.occ.synchronize()
            gmsh.model.geo.synchronize()

            trenches_3D = []
            for cur_trench in chip_and_trenches_map[1:]:
                trenches_3D += [y[1] for y in cur_trench if y[0]==3]
            chip_tags = [x[1] for x in chip_and_trenches_map[0] if x[1] not in trenches_3D and x[0]==3]
            leDielectric = gmsh.model.addPhysicalGroup(3, chip_tags, name = 'dielectric_substrate')
            
            all_surfaces = gmsh.model.occ.getEntities(2)
            leDielectricGaps = []
            for dim, tag in all_surfaces:
                com = gmsh.model.occ.getCenterOfMass(dim, tag)
                if com[2] < -0.25*trench_depth and com[2] > -1.25*trench_depth:
                    leDielectricGaps.append(tag)
            air_box_tags += trenches_3D
        else:
            leDielectric = gmsh.model.addPhysicalGroup(3, [tag3D_chip_base], name = 'dielectric_substrate')
        
        leDielectricGaps = gmsh.model.addPhysicalGroup(2, leDielectricGaps, name = 'dielectric_gaps')
        leAirBox = gmsh.model.addPhysicalGroup(3, air_box_tags, name = 'air_box')
        

        ret_dict = {
            'air_box'   : self._conv_to_lists(leAirBox),
            'metals'    : lemetals_physgrps,
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
        #define lists to store points and lines
        points = []
        lines = []
    
        #create 2D points in gmsh
        for coord in coords:
            point = gmsh.model.occ.addPoint(coord[0], coord[1], self._geom_processor.chip_centre[2])
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
        if bottom_grounded: 
            z_point = self.center_z - np.abs(self.size_z)
            air_box_delta_z = 2 * np.abs(self.size_z)
        else:
            z_point = self.center_z - 2 * np.abs(self.size_z)
            air_box_delta_z = 3 * np.abs(self.size_z)

        air_box = gmsh.model.occ.addBox(x_point, y_point, z_point, self.size_x + air_box_delta_x, self.size_y + air_box_delta_y, air_box_delta_z)
        gmsh.model.occ.synchronize()
        self._extents = np.array([(x_point, x_point + self.size_x + air_box_delta_x),
                                  (y_point, y_point + self.size_y + air_box_delta_y),
                                  (z_point, z_point + air_box_delta_z)])

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

