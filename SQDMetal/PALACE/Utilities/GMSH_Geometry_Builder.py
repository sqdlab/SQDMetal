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

        #Handy tid-bits to use in debugging if using assert False and Notebooks...
        # import gmsh
        # gmsh.fltk.run()

        #Do pre-processing in shapely to get metallic elements and dielectric cutouts ready to build in GMSH.
        #Note: metals list contains ground plane. Dielectric gaps are the difference between the dielectric cutout
        #and the metals
        metals, dielectric_gaps = self._geom_processor.process_layers(metallic_layers, ground_plane, fuse_threshold=fuse_threshold, fillet_resolution=self.fillet_resolution, unit_conv=1e-3, **kwargs)    #It's in mm...

        self._unit_conv = 1e-3

        #Get dimensions of chip base and convert to design units in 'mm'
        chip_centre = self._geom_processor.chip_centre
        self.center_x = chip_centre[0] / self._unit_conv
        self.center_y = chip_centre[1] / self._unit_conv
        self.center_z = chip_centre[2] / self._unit_conv
        self.size_x = self._geom_processor.chip_size_x / self._unit_conv
        self.size_y = self._geom_processor.chip_size_y / self._unit_conv
        self.size_z = self._geom_processor.chip_size_z / self._unit_conv

        #Plot the shapely metals for user to see device
        # geoms = metals# + ground_plane
        # metal_names = [str(i) for i,_ in enumerate(geoms)]
        # gdf = gpd.GeoDataFrame({'names':metal_names}, geometry=geoms)
        # fig, ax = plt.subplots()
        # gdf.plot(ax = ax, column='names', cmap='tab10', alpha=0.75, categorical=True, legend=True)
        # plt.show()

        dict_all_entities = {}

        #Create substrate and air box
        dict_all_entities['chip'] = [(3,self._create_chip_base())]
        dict_all_entities['airbox'] = [(3,self._draw_air_box(kwargs.get('boundary_distances', {})))]
        #Draw shapely metal and dielectric gap polygons into GMSH. If the polygon has an interior, the interior sections are subtracted from
        #the exterior boundary
        metal_names = []
        for m,cur_metal in enumerate(self._create_gmsh_geometry_from_shapely_polygons(metals)):
            metal_name = 'metal_' + str(m)
            metal_names.append(metal_name)
            dict_all_entities[metal_name] = [cur_metal]
        dict_all_entities['dielectric_gaps'] = self._create_gmsh_geometry_from_shapely_polygons(dielectric_gaps)         #create all the gaps between the metal traces and the ground plane
        #Process simulation constructs - e.g. ports
        wave_ports = []
        sim_construct_names = []
        for cur_sim_construct in sim_constructs:
            if cur_sim_construct['type'] == 'wave':
                wave_ports.append(cur_sim_construct['metadata'])
                continue
            cur_key, cur_coords = cur_sim_construct['name'], cur_sim_construct['coords']
            sim_poly = self._draw_polygon_in_GMSH_from_coords(np.array(cur_coords)[:-1]*1e3) #i.e. ignore closed loop part as function closes it automatically and use mm
            sim_construct_names.append(cur_key)
            dict_all_entities[cur_key] = [(2,sim_poly)]
        #
        self._gmsh_sync()

        #Now cut out overlapping sections between simulation constructs (e.g. ports) and metals
        if len(sim_construct_names) > 0:
            cleanup_sim_construct_metal_overlap =  kwargs.get('cleanup_sim_construct_metal_overlap', 'preserve_sim_constructs')
            if cleanup_sim_construct_metal_overlap == 'preserve_sim_constructs':
                self._gmsh_cut_overlapping_parts(dict_all_entities, metal_names, sim_construct_names)
            elif cleanup_sim_construct_metal_overlap == 'preserve_metals':
                self._gmsh_cut_overlapping_parts(dict_all_entities, sim_construct_names, metal_names)
            else:
                assert cleanup_sim_construct_metal_overlap == '', "The argument cleanup_sim_construct_metal_overlap must be either '', 'preserve_sim_constructs' or 'preserve_metals'"
            self._gmsh_cut_overlapping_parts(dict_all_entities, ['dielectric_gaps'], sim_construct_names)

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

        #Process the metals into groups (extrude if required)
        full_3D_params = kwargs.get('full_3D_params', None)
        metal_tags_on_boundaries = []
        all_metal_3D_tags = []
        if full_3D_params and full_3D_params['metal_thickness'] > 0:
            metal_thickness = full_3D_params['metal_thickness'] * 1e3

            for metal_name in metal_names:
                dict_all_entities[metal_name] = gmsh.model.occ.extrude(dict_all_entities[metal_name], 0, 0, metal_thickness)
                self._gmsh_sync()
        else:
            dict_all_entities['metal_extrusions'] = []
            dict_all_entities['metals_on_boundaries'] = []
        self._gmsh_sync()

        #Process trenching into dielectric
        if full_3D_params is not None and full_3D_params['substrate_trenching'] > 0:
            trench_depth = full_3D_params['substrate_trenching'] * 1e3

            trenches = gmsh.model.occ.extrude(dict_all_entities['dielectric_gaps'], 0, 0, -trench_depth)
            self._gmsh_sync()
            dict_all_entities['temp_trenches'] = [x for x in trenches if x[0] == 3]

            self._gmsh_cut_overlapping_parts(dict_all_entities, ['chip'], ['temp_trenches'], remove_tool=True)
        self._gmsh_sync()

        lePortPolys = {}

        #Process any wave-ports last
        uc = self._unit_conv
        waveport_surface_tags = []
        wave_port_names = []
        for cur_wave_port in wave_ports:
            if cur_wave_port['plane'] == 'x_neg':
                centre = self._extents[0][0], cur_wave_port['y']/uc, cur_wave_port['z']/uc
                v1,v2 = [0,cur_wave_port['size_y']/2/uc,0], [0,0,cur_wave_port['size_z']/2/uc]
            elif cur_wave_port['plane'] == 'x_pos':
                centre = self._extents[0][1], cur_wave_port['y']/uc, cur_wave_port['z']/uc
                v1,v2 = [0,cur_wave_port['size_y']/2/uc,0], [0,0,cur_wave_port['size_z']/2/uc]
            elif cur_wave_port['plane'] == 'y_neg':
                centre = cur_wave_port['x']/uc, self._extents[1][0], cur_wave_port['z']/uc
                v1,v2 = [cur_wave_port['size_x']/2/uc,0,0], [0,0,cur_wave_port['size_z']/2/uc]
            elif cur_wave_port['plane'] == 'y_pos':
                centre = cur_wave_port['x']/uc, self._extents[1][1], cur_wave_port['z']/uc
                v1,v2 = [cur_wave_port['size_x']/2/uc,0,0], [0,0,cur_wave_port['size_z']/2/uc]
            elif cur_wave_port['plane'] == 'z_neg':
                centre = cur_wave_port['x']/uc, cur_wave_port['y']/uc, self._extents[2][0]
                v1,v2 = [cur_wave_port['size_x']/2/uc,0,0], [0,cur_wave_port['size_y']/2/uc,0]
            elif cur_wave_port['plane'] == 'z_pos':
                centre = cur_wave_port['x']/uc, cur_wave_port['y']/uc, self._extents[2][1]
                v1,v2 = [cur_wave_port['size_x']/2/uc,0,0], [0,cur_wave_port['size_y']/2/uc,0]

            if cur_wave_port['shape'] == 'rectangle':
                new_wvprt = self._create_gmsh_plane(centre, v1, v2, message="Rectangular waveport")
            #
            wave_port_name = cur_wave_port['port_name']
            wave_port_names.append(wave_port_name)
            dict_all_entities[wave_port_name] = [(2,new_wvprt)]
            self._gmsh_sync()
            
            #TODO: Should this also include lumped ports - but they shouldn't overlap; so perhaps have a check that they don't intersect?
            # self._gmsh_fragment_make_conformal(dict_all_entities, metal_names+['dielectric_gaps', 'metal_extrusions', 'chip', 'airbox'], [wave_port_name])
            waveport_surface_tags += [x[1] for x in dict_all_entities[wave_port_name] if x[0] == 2]

        #Create Farfield Planes
        centre = self._extents[0][0], self._ff_centre[1], self._ff_centre[2]
        v1,v2 = [0,self._ff_sizes[1]/2,0], [0,0,self._ff_sizes[2]/2]
        dict_all_entities['ffxNeg'] = [(2,self._create_gmsh_plane(centre, v1, v2, message="Farfield"))]
        #
        centre = self._extents[0][1], self._ff_centre[1], self._ff_centre[2]
        v1,v2 = [0,self._ff_sizes[1]/2,0], [0,0,self._ff_sizes[2]/2]
        dict_all_entities['ffxPos'] = [(2,self._create_gmsh_plane(centre, v1, v2, message="Farfield"))]
        #
        centre = self._ff_centre[0], self._extents[1][0], self._ff_centre[2]
        v1,v2 = [self._ff_sizes[0]/2,0,0], [0,0,self._ff_sizes[2]/2]
        dict_all_entities['ffyNeg'] = [(2,self._create_gmsh_plane(centre, v1, v2, message="Farfield"))]
        #
        centre = self._ff_centre[0], self._extents[1][1], self._ff_centre[2]
        v1,v2 = [self._ff_sizes[0]/2,0,0], [0,0,self._ff_sizes[2]/2]
        dict_all_entities['ffyPos'] = [(2,self._create_gmsh_plane(centre, v1, v2, message="Farfield"))]
        #
        centre = self._ff_centre[0], self._ff_centre[1], self._extents[2][0]
        v1,v2 = [self._ff_sizes[0]/2,0,0], [0,self._ff_sizes[1]/2,0]
        dict_all_entities['ffzNeg'] = [(2,self._create_gmsh_plane(centre, v1, v2, message="Farfield"))]
        #
        centre = self._ff_centre[0], self._ff_centre[1], self._extents[2][1]
        v1,v2 = [self._ff_sizes[0]/2,0,0], [0,self._ff_sizes[1]/2,0]
        dict_all_entities['ffzPos'] = [(2,self._create_gmsh_plane(centre, v1, v2, message="Farfield"))]
        #
        ff_names = ['ffxNeg', 'ffxPos', 'ffyNeg', 'ffyPos', 'ffzNeg', 'ffzPos']

        #Consolidate all overlapping surfaces
        #NOTE: Strange embedding issue with fix here: https://gitlab.onelab.info/gmsh/gmsh/-/issues/2669
        # - Basically consider a port that's off the chip and only in the air-box
        # - The extrusions force remapping of the air-box and its adjacent surfaces - fine...
        # - But it fails to return the remapping of the port unless it's added first - i.e. 1D, then 2D, then 3D?!
        self._gmsh_fragment_make_conformal(dict_all_entities, ['chip'], ['airbox'], filter_tags=True)
        self._gmsh_fragment_make_conformal(dict_all_entities, ['airbox','chip'],metal_names)
        if full_3D_params is not None and full_3D_params['substrate_trenching'] > 0:
            # self._gmsh_fragment_make_conformal(dict_all_entities, ['airbox','chip']+metal_names,ff_names+wave_port_names,[])          
            # dict_all_entities['chip'] = [x for x in dict_all_entities['chip'] if not (x[0] == 2 and x[1] in all_dielectric_gaps)]
            all_surfaces = gmsh.model.occ.getEntities(2)
            all_dielectric_gaps = []
            dict_all_entities['dielectric_gaps'] = []
            for dim, tag in all_surfaces:
                com = gmsh.model.occ.getCenterOfMass(dim, tag)
                if com[2] < -0.25*trench_depth and com[2] > -1.25*trench_depth and np.min(np.abs(self._extents.T - np.array(com))) > 1e-6:
                    dict_all_entities['dielectric_gaps'].append((dim, tag))
                    all_dielectric_gaps.append(tag)
        self._gmsh_fragment_make_conformal(dict_all_entities, ['chip', 'airbox']+metal_names,ff_names+['dielectric_gaps']+sim_construct_names+wave_port_names,[])

        metal_tags_on_boundaries = []
        for metal_name in metal_names:
            all_metal_3D_tags = [x[1] for x in dict_all_entities[metal_name] if x[0] == 3]
            if len(all_metal_3D_tags) > 0:
                all_metal_2D_tags = []
                for cur_tag_3D in all_metal_3D_tags:
                    all_metal_2D_tags += gmsh.model.getAdjacencies(3, cur_tag_3D)[1].tolist()
                all_metal_2D_tags = [(2,x) for x in all_metal_2D_tags]
            else:
                all_metal_2D_tags = dict_all_entities[metal_name]
            #Remove planes touching the walls of the simulation
            allowed_tags = []
            for x in all_metal_2D_tags:
                com = gmsh.model.occ.getCenterOfMass(x[0], x[1])
                if np.min(np.abs(self._extents.T - np.array(com))) > 1e-6:
                    allowed_tags.append(x)
                else:
                    metal_tags_on_boundaries.append(x[1])
            dict_all_entities[metal_name] = allowed_tags

        all_metal_tags = []
        for metal_name in metal_names:
            all_metal_tags += [x[1] for x in dict_all_entities[metal_name]]

        # metal_tags_on_boundaries = [x[1] for x in dict_all_entities['metals_on_boundaries']]
        if len(wave_ports) > 0:
            #Anything on the metal overlapping with the waveport belongs to the metal
            for cur_wave_port in wave_port_names:
                dict_all_entities[cur_wave_port] = [x for x in dict_all_entities[cur_wave_port] if not (x[0]==2 and x[1] in metal_tags_on_boundaries)]
        # print(metal_tags_on_boundaries)
        waveport_surface_tags = []
        for wave_port_name in wave_port_names:
            waveport_surface_tags += [x[1] for x in dict_all_entities[wave_port_name] if x[0] == 2]

        #Create a physical group for each port/junction simulation construct...
        for cur_sim_constr in sim_construct_names:
            lePortPolys[cur_sim_constr] = gmsh.model.addPhysicalGroup(2, [x[1] for x in dict_all_entities[cur_sim_constr]], name = cur_sim_constr)
        #Create physical groups for metals
        lemetals_physgrps = []
        for metal_name in metal_names:
            metal_physical_group = gmsh.model.addPhysicalGroup(2, [x[1] for x in dict_all_entities[metal_name]], name = metal_name)
            lemetals_physgrps.append(metal_physical_group)
        leAirBox_entities = []
        chip_entities = [x[1] for x in dict_all_entities['chip']]
        for cur_entity in dict_all_entities['airbox']:
            if cur_entity[1] not in chip_entities:
                leAirBox_entities.append(cur_entity[1])
        leAirBox = gmsh.model.addPhysicalGroup(3, leAirBox_entities, name = 'air_box')
        leDielectric = gmsh.model.addPhysicalGroup(3, [x[1] for x in dict_all_entities['chip']], name = 'dielectric_substrate')       
        leDielectricGaps = gmsh.model.addPhysicalGroup(2, [x[1] for x in dict_all_entities['dielectric_gaps']], name = 'dielectric_gaps')
        for wave_port_name in wave_port_names:
            lePortPolys[wave_port_name] = gmsh.model.addPhysicalGroup(2, [x[1] for x in dict_all_entities[wave_port_name]], name = wave_port_name)

        ff_physgrps = {}
        for ff_name in ff_names:
            dict_all_entities[ff_name] = [x for x in dict_all_entities[ff_name] if not (x[0]==2 and (x[1] in metal_tags_on_boundaries or x[1] in waveport_surface_tags))]
            ff_physgrps[ff_name] = gmsh.model.addPhysicalGroup(2, [x[1] for x in dict_all_entities[ff_name]], name = ff_name)

        gmsh.model.occ.removeAllDuplicates()

        ret_dict = {
            'air_box'   : self._conv_to_lists(leAirBox),
            'metals'    : lemetals_physgrps,
            'metalsShapely'    : metals,
            'far_field' : ff_physgrps,
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

    def _gmsh_sync(self):
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()

    def _gmsh_fragment_make_conformal(self, all_entities:dict, main_objects:list, tool_objects:list, filter_tags=False):
        '''
        This function runs gmsh.model.occ.fragment in a nicer wrapper. The tool_objects are carved out of the main_objects so that the
        final geometry remains conformal.

        Inputs:
            all_entities - dictionary of all geometric objects labelled by their keys. The values are the Gmsh geometry (dim,tag) pairs
            main_objects - list of strings representing the keys of the objects from which to carve out
            tool_objects - list of strings representing the keys of the objects that are to be carved out
            filter_tags  - If True, it will go through the tool_objects to remove all dim,tag pairs found in common within main_objects.
                           Only really relevant if the dimensions are the same (e.g. carving a 3D chunk out of a 3D object etc.)
        
        Outputs:
            None
        
        The dictionary all_entities is modified with the new geometry tags
        '''
        main_dim_tags = []
        for x in main_objects:
            main_dim_tags += sorted(all_entities[x], key=lambda tup: tup[0])
        tool_dim_tags = []
        for x in tool_objects:
            tool_dim_tags += all_entities[x]

        #See here: https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/tutorials/python/t16.py for details on the mapping...
        main_objs, leMap = gmsh.model.occ.fragment(main_dim_tags, tool_dim_tags, removeObject=True, removeTool=True)

        #Using leMap (which gives the list of dim,tag pairs that the original dim,tag pairs in target_dim_tags map onto), redo the main_objects and
        #tool_objects while making sure to keep the appropriate names:
        cur_index = 0
        for cur_target in (main_objects + tool_objects):
            cur_target_dimTags = []
            for cur_dimTag in all_entities[cur_target]:
                cur_target_dimTags += leMap[cur_index]
                cur_index += 1
            all_entities[cur_target] = sorted(cur_target_dimTags, key=lambda tup: tup[0])

        if filter_tags:
            main_dim_tags = []
            for x in main_objects:
                main_dim_tags += all_entities[x]
            for cur_tool in tool_objects:
                filtered_dim_tags = []
                for cur_dim_tag in all_entities[cur_tool]:
                    if not cur_dim_tag in main_dim_tags:
                        filtered_dim_tags.append(cur_dim_tag)
                all_entities[cur_tool] = filtered_dim_tags

        self._gmsh_sync()

    def _gmsh_cut_overlapping_parts(self, all_entities:dict, target_objects:list, tool_objects:list, remove_tool=False):
        '''
        This function runs gmsh.model.occ.cut in a nicer wrapper. The tool_objects are carved out of the target_objects while the
        tool_objects remain.

        Inputs:
            all_entities   - dictionary of all geometric objects labelled by their keys. The values are the Gmsh geometry (dim,tag) pairs
            target_objects - list of strings representing the keys of the objects from which to cut out
            tool_objects   - list of strings representing the keys of the objects that are to cut out from the target_objects
        
        Outputs:
            None
        
        The dictionary all_entities is modified with the new geometry tags
        '''
        target_dim_tags = []
        for x in target_objects:
            target_dim_tags += all_entities[x]
        tool_dim_tags = []
        for x in tool_objects:
            tool_dim_tags += all_entities[x]

        dimTags,leMap = gmsh.model.occ.cut(target_dim_tags, tool_dim_tags, removeObject=True, removeTool=remove_tool)

        #Using leMap (which gives the list of dim,tag pairs that the original dim,tag pairs in target_dim_tags map onto), redo the target target_objects
        #making sure to keep the appropriate names
        cur_index = 0
        for cur_target in target_objects:
            cur_target_dimTags = []
            for cur_dimTag in all_entities[cur_target]:
                cur_target_dimTags += leMap[cur_index]
                cur_index += 1
            all_entities[cur_target] = cur_target_dimTags

        self._gmsh_sync()


    def _create_gmsh_plane(self, centre_vec, v1, v2, message=""):
        centre_vec = np.array(centre_vec)
        v1 = np.array(v1)
        v2 = np.array(v2)
        coords = [centre_vec-v1-v2, centre_vec+v1-v2, centre_vec+v1+v2, centre_vec-v1+v2]
        coords = np.array(coords)
        assert np.min(coords[:,0]) >= self._extents[0][0] and np.max(coords[:,0]) <= self._extents[0][1], f"Element {message} leaks outside x-axis boundary."
        assert np.min(coords[:,1]) >= self._extents[1][0] and np.max(coords[:,1]) <= self._extents[1][1], f"Element {message} leaks outside y-axis boundary."
        assert np.min(coords[:,2]) >= self._extents[2][0] and np.max(coords[:,2]) <= self._extents[2][1], f"Element {message} leaks outside z-axis boundary."
        return self._draw_polygon_in_GMSH_from_coords(coords)

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
        coords = np.array(coords)

        #define list to store lines
        lines = []
    
        #create 2D points in gmsh
        if coords.shape[1] == 2:
            points = [gmsh.model.occ.addPoint(coord[0], coord[1], self._geom_processor.chip_centre[2]) for coord in coords]
        else:
            points = [gmsh.model.occ.addPoint(coord[0], coord[1], coord[2]) for coord in coords]

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

        return chip_base
    

    def _draw_air_box(self, boundary_distances):
        '''Creates the airbox which the chip will be placed in.

        Inputs:
            boundary_distances - dictionary containing distances of the boundaries with respect to the chip. For example, the
                                 keys x_prop/y_prop set the boundaries as a proportion of the x/y chip sizes, whereas the keys
                                 x_neg/x_pos give the absolute distances along the x axes on the negative/positve sides (it's
                                 the same for the y and z axes). One exception is for proportions along the z-axis, where it's
                                 z_prop_top and z_prop_bottom (i.e. flexibility in having the chip grounded on its bottom etc.)

        Returns:
            Tuple of integers with the first value specifying the dimension and the second value the ID of the object in GMSH. 
        '''

        if 'x_prop' in boundary_distances:
            xMin = self.center_x - self.size_x/2 - boundary_distances['x_prop']*self.size_x
            xMax = self.center_x + self.size_x/2 + boundary_distances['x_prop']*self.size_x
        else:
            xMin = self.center_x - self.size_x/2 - boundary_distances['x_neg']/self._unit_conv
            xMax = self.center_x + self.size_x/2 + boundary_distances['x_pos']/self._unit_conv

        if 'y_prop' in boundary_distances:
            yMin = self.center_y - self.size_y/2 - boundary_distances['y_prop']*self.size_y
            yMax = self.center_y + self.size_y/2 + boundary_distances['y_prop']*self.size_y
        else:
            yMin = self.center_y - self.size_y/2 - boundary_distances['y_neg']/self._unit_conv
            yMax = self.center_y + self.size_y/2 + boundary_distances['y_pos']/self._unit_conv

        if 'z_prop_top' in boundary_distances:
            zMin = self.center_z - self.size_z - boundary_distances['z_prop_bottom']*self.size_z
            zMax = self.center_z + boundary_distances['z_prop_top']*self.size_z
        else:
            zMin = self.center_z - self.size_z - boundary_distances['z_neg']/self._unit_conv
            zMax = self.center_z + boundary_distances['z_pos']/self._unit_conv

        air_box = gmsh.model.occ.addBox(xMin, yMin, zMin, xMax-xMin, yMax-yMin, zMax-zMin)
        self._extents = np.array([(xMin,xMax), (yMin,yMax), (zMin,zMax)])
        self._ff_centre = np.array([(xMin+xMax)/2, (yMin+yMax)/2, (zMin+zMax)/2])
        self._ff_sizes = np.array([(xMax-xMin), (yMax-yMin), (zMax-zMin)])

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

