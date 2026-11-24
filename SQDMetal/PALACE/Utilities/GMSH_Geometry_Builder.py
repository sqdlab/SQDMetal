# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0

from SQDMetal.Utilities.ShapelyEx import ShapelyEx
from SQDMetal.Utilities.QUtilities import QUtilities
import gmsh
import numpy as np
import shapely

from SQDMetal.Utilities.GeometryProcessors.GeomBase import GeomBase
from SQDMetal.Utilities.GeometryProcessors.GeomQiskitMetal import GeomQiskitMetal
from SQDMetal.PALACE.Utilities.GMSH_WireBond_Builder import GMSH_WireBond_Builder
# import matplotlib.pyplot as plt
# import geopandas as gpd
class GMSH_Geometry_Builder:

    def __init__(self, geom_processor:GeomBase, fillet_resolution, gmsh_verbosity=1):
        # gmsh_verbosity: 0: silent except for fatal errors, 1: +errors, 2: +warnings, 3: +direct, 4: +information, 5: +status, 99: +debug

        self._geom_processor = geom_processor
        self.fillet_resolution = fillet_resolution

        # Initialize the GMSH API (from a clean state) and name the model. Only tear down a running GMSH session: calling
        # model.remove()/finalize() on a finalised session raises 'Could not get last error' (e.g. on a second
        # simulation in the same Jupyter kernel).
        if gmsh.isInitialized():
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

        # Handy tid-bits to use in debugging if using assert False and Notebooks...
        # import gmsh
        # gmsh.fltk.run()

        self._unit_conv = 1e-3

        # Do pre-processing in shapely to get metallic elements and dielectric cutouts ready to build in GMSH.
        # Note: metals list contains ground plane. Dielectric gaps are the difference between the dielectric cutout
        # and the metals
        metals, dielectric_gaps = self._geom_processor.process_layers(metallic_layers, ground_plane, fuse_threshold=fuse_threshold, fillet_resolution=self.fillet_resolution, unit_conv=self._unit_conv, **kwargs)    #It's in mm...
        metals_div = [{f'metal_{m}':x} for m,x in enumerate(metals)]
        gaps_div = [{f'gap_{m}':x} for m,x in enumerate(dielectric_gaps)]

        #If integration areas are specified, the metal plus dielectric-gaps are subdivided. Note that these subdivisions
        #should keep Polygons as Polygons - i.e. can't create holes from just partitioned cuts...
        if 'integration_areas' in kwargs:
            int_areas = []
            for m,cur_area in enumerate(kwargs['integration_areas']):
                cur_area_scaled = shapely.affinity.scale( cur_area, xfact=1/self._unit_conv, yfact=1/self._unit_conv, origin=(0, 0) )
                metals_div, gaps_div, int_areas = self._divide_surface_with_area(cur_area_scaled, metals_div, gaps_div, int_areas)

        # Get dimensions of chip base and convert to design units in 'mm'
        chip_centre = self._geom_processor.chip_centre
        self.center_x = chip_centre[0] / self._unit_conv
        self.center_y = chip_centre[1] / self._unit_conv
        self.center_z = chip_centre[2] / self._unit_conv
        self.size_x = self._geom_processor.chip_size_x / self._unit_conv
        self.size_y = self._geom_processor.chip_size_y / self._unit_conv
        self.size_z = self._geom_processor.chip_size_z / self._unit_conv

        # # Plot the shapely metals for user to see device
        # geoms = metals# + ground_plane
        # metal_names = [str(i) for i,_ in enumerate(geoms)]
        # gdf = gpd.GeoDataFrame({'names':metal_names}, geometry=geoms)
        # fig, ax = plt.subplots()
        # gdf.plot(ax = ax, column='names', cmap='tab10', alpha=0.75, categorical=True, legend=True)
        # plt.show()

        dict_all_entities = {}

        # Create substrate and air box
        dict_all_entities['chip'] = [(3,self._create_chip_base())]
        dict_all_entities['airbox'] = [(3,self._draw_air_box(kwargs.get('boundary_distances', {})))]
        # Draw shapely metal and dielectric gap polygons into GMSH. If the polygon has an interior, the interior sections are subtracted from
        # the exterior boundary
        metal_names = []
        for cur_metal in metals_div:
            for metal_name in cur_metal:
                cur_gmsh_poly = self._create_gmsh_geometry_from_shapely_polygons([cur_metal[metal_name]])
                metal_names.append(metal_name)
                dict_all_entities[metal_name] = cur_gmsh_poly
        gap_names = []
        for cur_gap in gaps_div:
            for gap_name in cur_gap:
                cur_gmsh_poly = self._create_gmsh_geometry_from_shapely_polygons([cur_gap[gap_name]])
                gap_names.append(gap_name)
                dict_all_entities[gap_name] = cur_gmsh_poly
        # Process simulation constructs - e.g. ports
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

        # Now cut out overlapping sections between simulation constructs (e.g. ports) and metals
        if len(sim_construct_names) > 0:
            cleanup_sim_construct_metal_overlap =  kwargs.get('cleanup_sim_construct_metal_overlap', 'preserve_sim_constructs')
            if cleanup_sim_construct_metal_overlap == 'preserve_sim_constructs':
                self._gmsh_cut_overlapping_parts(dict_all_entities, metal_names, sim_construct_names)
            elif cleanup_sim_construct_metal_overlap == 'preserve_metals':
                self._gmsh_cut_overlapping_parts(dict_all_entities, sim_construct_names, metal_names)
            else:
                assert cleanup_sim_construct_metal_overlap == '', "The argument cleanup_sim_construct_metal_overlap must be either '', 'preserve_sim_constructs' or 'preserve_metals'"
            self._gmsh_cut_overlapping_parts(dict_all_entities, gap_names, sim_construct_names)

        fine_mesh_elems = []
        for cur_fine_mesh in fine_meshes:
            if cur_fine_mesh['type'] == 'box':
                cur_mesh_attrb = {}
                x1, x2 = cur_fine_mesh['x_bnds']
                y1, y2 = cur_fine_mesh['y_bnds']
                # Get it in mm...
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
                # Get it in mm...
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

        # Process the metals into groups (extrude if required)
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

        # Process trenching into dielectric
        if full_3D_params is not None and full_3D_params['substrate_trenching'] > 0:
            trench_depth = full_3D_params['substrate_trenching'] * 1e3

            gap_entities = []
            for cur_gap in gap_names:
                gap_entities += dict_all_entities[cur_gap]
            # preserve ports over trenched regions
            port_entities = [x for cur_name in sim_construct_names for x in dict_all_entities[cur_name] if x[0] == 2]
            if len(port_entities) > 0:
                gap_entities = gap_entities + gmsh.model.occ.copy(port_entities)
            trenches = gmsh.model.occ.extrude(gap_entities, 0, 0, -trench_depth)
            self._gmsh_sync()
            dict_all_entities['temp_trenches'] = [x for x in trenches if x[0] == 3]

            self._gmsh_cut_overlapping_parts(dict_all_entities, ['chip'], ['temp_trenches'], remove_tool=True)
        self._gmsh_sync()

        # Build 3D wire-bonds (they sit in the air-box and land on the metals/substrate)
        wb_names, wb_meta = self._build_wirebonds(kwargs.get('wirebonds', []), metals, ground_plane, full_3D_params, dict_all_entities)

        lePortPolys = {}

        # Process any wave-ports last
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

            # TODO: Should this also include lumped ports - but they shouldn't overlap; so perhaps have a check that they don't intersect?
            # self._gmsh_fragment_make_conformal(dict_all_entities, metal_names+['dielectric_gaps', 'metal_extrusions', 'chip', 'airbox'], [wave_port_name])
            waveport_surface_tags += [x[1] for x in dict_all_entities[wave_port_name] if x[0] == 2]

        # Create Farfield Planes
        centre = self._extents[0][0], self._ff_centre[1], self._ff_centre[2]
        v1,v2 = [0,self._ff_sizes[1]/2,0], [0,0,self._ff_sizes[2]/2]
        dict_all_entities['x_neg'] = [(2,self._create_gmsh_plane(centre, v1, v2, message="Farfield"))]
        #
        centre = self._extents[0][1], self._ff_centre[1], self._ff_centre[2]
        v1,v2 = [0,self._ff_sizes[1]/2,0], [0,0,self._ff_sizes[2]/2]
        dict_all_entities['x_pos'] = [(2,self._create_gmsh_plane(centre, v1, v2, message="Farfield"))]
        #
        centre = self._ff_centre[0], self._extents[1][0], self._ff_centre[2]
        v1,v2 = [self._ff_sizes[0]/2,0,0], [0,0,self._ff_sizes[2]/2]
        dict_all_entities['y_neg'] = [(2,self._create_gmsh_plane(centre, v1, v2, message="Farfield"))]
        #
        centre = self._ff_centre[0], self._extents[1][1], self._ff_centre[2]
        v1,v2 = [self._ff_sizes[0]/2,0,0], [0,0,self._ff_sizes[2]/2]
        dict_all_entities['y_pos'] = [(2,self._create_gmsh_plane(centre, v1, v2, message="Farfield"))]
        #
        centre = self._ff_centre[0], self._ff_centre[1], self._extents[2][0]
        v1,v2 = [self._ff_sizes[0]/2,0,0], [0,self._ff_sizes[1]/2,0]
        dict_all_entities['z_neg'] = [(2,self._create_gmsh_plane(centre, v1, v2, message="Farfield"))]
        #
        centre = self._ff_centre[0], self._ff_centre[1], self._extents[2][1]
        v1,v2 = [self._ff_sizes[0]/2,0,0], [0,self._ff_sizes[1]/2,0]
        dict_all_entities['z_pos'] = [(2,self._create_gmsh_plane(centre, v1, v2, message="Farfield"))]
        #
        ff_names = ['x_neg', 'x_pos', 'y_neg', 'y_pos', 'z_neg', 'z_pos']

        # Consolidate all overlapping surfaces
        # NOTE: Strange embedding issue with fix here: https://gitlab.onelab.info/gmsh/gmsh/-/issues/2669
        # - Basically consider a port that's off the chip and only in the air-box
        # - The extrusions force remapping of the air-box and its adjacent surfaces - fine...
        # - But it fails to return the remapping of the port unless it's added first - i.e. 1D, then 2D, then 3D?!
        self._gmsh_fragment_make_conformal(dict_all_entities, ['chip'], ['airbox'], filter_tags=True)
        # Only the wire-bond solids are fragmented here; their 2D parts (sheets/pads) join in the final fragmentation
        wb_names_2D = self._split_wirebond_2D_parts(dict_all_entities, wb_names)
        self._gmsh_fragment_make_conformal(dict_all_entities, ['airbox','chip'],metal_names+wb_names)
        self._merge_wirebond_2D_parts(dict_all_entities, wb_names)
        self._filter_wirebond_volumes_from_airbox(dict_all_entities, wb_names)
        if full_3D_params is not None and full_3D_params['substrate_trenching'] > 0:
            # self._gmsh_fragment_make_conformal(dict_all_entities, ['airbox','chip']+metal_names,ff_names+wave_port_names,[])
            assert len(kwargs.get('integration_areas',[])) == 0, "Integration areas is not supported when using trenching..."
            all_surfaces = gmsh.model.occ.getEntities(2)
            dict_all_entities['dielectric_gaps'] = []
            for dim, tag in all_surfaces:
                com = gmsh.model.occ.getCenterOfMass(dim, tag)
                if com[2] < -0.25*trench_depth and com[2] > -1.25*trench_depth and np.min(np.abs(self._extents.T - np.array(com))) > 1e-6:
                    dict_all_entities['dielectric_gaps'].append((dim, tag))
            gap_names = ['dielectric_gaps']
        # 2D wire-bond parts (sheets, off-chip pads) go on the tool side (like ports); OCC otherwise drops their mapping
        wb_names_2D = self._split_wirebond_2D_parts(dict_all_entities, wb_names)
        self._gmsh_fragment_make_conformal(dict_all_entities, ['chip', 'airbox']+metal_names+wb_names,ff_names+gap_names+sim_construct_names+wave_port_names+wb_names_2D,[])
        self._merge_wirebond_2D_parts(dict_all_entities, wb_names)
        self._filter_wirebond_volumes_from_airbox(dict_all_entities, wb_names)

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
            # Remove planes touching the walls of the simulation
            allowed_tags = []
            for x in all_metal_2D_tags:
                com = gmsh.model.occ.getCenterOfMass(x[0], x[1])
                if np.min(np.abs(self._extents.T - np.array(com))) > 1e-6:
                    allowed_tags.append(x)
                else:
                    metal_tags_on_boundaries.append(x[1])
            dict_all_entities[metal_name] = allowed_tags

        # Convert wire-bond solids into their bounding PEC/conductor surfaces (the interior is not meshed)
        self._finalise_wirebond_surfaces(dict_all_entities, wb_names, metal_names, gap_names, sim_construct_names, wb_meta)

        # metal_tags_on_boundaries = [x[1] for x in dict_all_entities['metals_on_boundaries']]
        if len(wave_ports) > 0:
            # Anything on the metal overlapping with the waveport belongs to the metal
            for cur_wave_port in wave_port_names:
                dict_all_entities[cur_wave_port] = [x for x in dict_all_entities[cur_wave_port] if not (x[0]==2 and x[1] in metal_tags_on_boundaries)]
        # print(metal_tags_on_boundaries)
        waveport_surface_tags = []
        for wave_port_name in wave_port_names:
            waveport_surface_tags += [x[1] for x in dict_all_entities[wave_port_name] if x[0] == 2]

        # Create a physical group for each port/junction simulation construct...
        for cur_sim_constr in sim_construct_names:
            lePortPolys[cur_sim_constr] = gmsh.model.addPhysicalGroup(2, [x[1] for x in dict_all_entities[cur_sim_constr]], name = cur_sim_constr)
        # Create physical groups for metals
        phys_grp_ints = {}
        lemetals_physgrps = []
        contiguous_metal_mapping = []
        for m,cur_metal in enumerate(metals_div):
            cur_metal_pieces = []
            for cur_piece in cur_metal:
                cur_entities = [x[1] for x in dict_all_entities[cur_piece]]
                metal_physical_group = gmsh.model.addPhysicalGroup(2, cur_entities, name = cur_piece)
                lemetals_physgrps.append(metal_physical_group)
                cur_metal_pieces.append(metal_physical_group)
                phys_grp_ints[cur_piece] = metal_physical_group
            contiguous_metal_mapping.append(cur_metal_pieces)
        wirebond_attrs = []
        for wb_name in wb_names:
            meta = wb_meta[wb_name]
            surfs = [x[1] for x in dict_all_entities[wb_name] if x[0] == 2]
            if len(surfs) == 0:
                print(f"WARNING: Wire-bond {meta['name']} has no surfaces left after meshing preparation; it is ignored.")
                continue
            wb_physgrp = gmsh.model.addPhysicalGroup(2, surfs, name=meta['name'])
            lemetals_physgrps.append(wb_physgrp)    #So that wire-bond-unaware simulation types still treat it as PEC
            meta['physical_group'] = wb_physgrp
            meta['surfaces'] = surfs
            wirebond_attrs.append(meta)
        # Conductors joined by wire-bonds become a single conductor (e.g. a single capacitance terminal)
        metals_conductors = metals
        conductor_groups = [[m] for m in range(len(metals))]
        if len(wirebond_attrs) > 0:
            contiguous_metal_mapping, metals_conductors, conductor_groups = self._merge_wirebonded_conductors(contiguous_metal_mapping, metals, wirebond_attrs)
        leAirBox_entities = []
        chip_entities = [x[1] for x in dict_all_entities['chip']]
        for cur_entity in dict_all_entities['airbox']:
            if cur_entity[1] not in chip_entities:
                leAirBox_entities.append(cur_entity[1])
        leAirBox = gmsh.model.addPhysicalGroup(3, leAirBox_entities, name = 'air_box')
        leDielectric = gmsh.model.addPhysicalGroup(3, [x[1] for x in dict_all_entities['chip']], name = 'dielectric_substrate')       
        for wave_port_name in wave_port_names:
            lePortPolys[wave_port_name] = gmsh.model.addPhysicalGroup(2, [x[1] for x in dict_all_entities[wave_port_name]], name = wave_port_name)

        ff_physgrps = {}
        for ff_name in ff_names:
            dict_all_entities[ff_name] = [x for x in dict_all_entities[ff_name] if not (x[0]==2 and (x[1] in metal_tags_on_boundaries or x[1] in waveport_surface_tags or x[1] in self._wb_wall_faces))]
            ff_physgrps[ff_name] = gmsh.model.addPhysicalGroup(2, [x[1] for x in dict_all_entities[ff_name]], name = ff_name)

        gmsh.model.occ.removeAllDuplicates()

        #Partition the dielectric gaps as per the integration areas:
        #TODO: support overlapping areas (will require pre-paritioning into mutually exclusive areas and then combing them in post-processing).
        #TODO: support integration over metals (for superconductors, it's fine as it'll reject B-field anyway)
        #Process Dielectric Gaps
        leDielectricGaps = []
        if full_3D_params is not None and full_3D_params['substrate_trenching'] > 0:
            # With trenching, the z=0 gap surfaces were consumed when cutting the trenches; the substrate-air interface is
            # now the trench floors and walls (collected in 'dielectric_gaps'). Using the gap pieces here would give empty
            # physical groups (which Palace reports as 'Unknown boundary postprocessing attributes').
            leDielectricGaps.append(gmsh.model.addPhysicalGroup(2, [x[1] for x in dict_all_entities['dielectric_gaps'] if x[0] == 2], name='dielectric_gaps'))
        else:
            for cur_gap in gaps_div:
                for cur_piece in cur_gap:
                    cur_phys_grp = gmsh.model.addPhysicalGroup(2, [x[1] for x in dict_all_entities[cur_piece]], name=cur_piece)
                    phys_grp_ints[cur_piece] = cur_phys_grp
                    leDielectricGaps.append(cur_phys_grp)
        
        final_int_areas = []
        if 'integration_areas' in kwargs:
            for m,cur_area in enumerate(int_areas):
                cur_entities = []
                for cur_piece in cur_area:
                    cur_entities.append(phys_grp_ints[cur_piece])
                final_int_areas.append(cur_entities)

        for meta in wirebond_attrs:
            fm = meta['fine_mesh']
            if fm is not None:
                fine_mesh_elems.append({'region': [(2,x) for x in meta['surfaces']], **fm})

        ret_dict = {
            'wirebonds' : wirebond_attrs,
            'air_box'   : self._conv_to_lists(leAirBox),
            'metals'    : lemetals_physgrps,
            'contiguous_metal_mapping': contiguous_metal_mapping,
            'metalsShapely'    : metals_conductors,
            'conductor_groups': conductor_groups,
            'far_field' : ff_physgrps,
            'ports'     : lePortPolys,
            'dielectric': self._conv_to_lists(leDielectric),
            'dielectric_gaps': leDielectricGaps,
            'fine_mesh_elems': fine_mesh_elems,
            'integration_areas':final_int_areas
        }
        return ret_dict

    def _build_wirebonds(self, wirebonds, metals, ground_plane, full_3D_params, dict_all_entities):
        '''
        Classifies and builds the wire-bonds in GMSH. The wire-bond specifications are given in metres (as stored by
        the PALACE_Model classes) and converted into the GMSH units here.
        '''
        if len(wirebonds) == 0:
            return [], {}
        uc = self._unit_conv
        t_metal = full_3D_params['metal_thickness']/uc if (full_3D_params and full_3D_params['metal_thickness'] > 0) else 0.0
        chip_bnds = (self.center_x-self.size_x/2, self.center_y-self.size_y/2, self.center_x+self.size_x/2, self.center_y+self.size_y/2)
        wbb = GMSH_WireBond_Builder(metals, chip_bnds, self.center_z, t_metal, self._extents, not ground_plane.get('omit', True))
        len_keys = ['width', 'thickness', 'height', 'foot_length', 'foot_width', 'foot_thickness', 'offchip_z']
        mesh_keys = ['mesh_min', 'mesh_max', 'taper_dist_min', 'taper_dist_max']
        wb_names, wb_meta = [], {}
        for cur_wb in wirebonds:
            wb = dict(cur_wb)
            wb['pos_start'] = np.array(wb['pos_start'][:2], dtype=float)/uc
            wb['pos_end'] = np.array(wb['pos_end'][:2], dtype=float)/uc
            for k in len_keys:
                if wb.get(k, None) is not None:
                    wb[k] = wb[k]/uc
            wb['offchip_pad_size'] = [x/uc for x in wb['offchip_pad_size']]
            wb['profile_params'] = dict(wb['profile_params'])
            if wb['profile_params'].get('corner_radius', None) is not None:
                wb['profile_params']['corner_radius'] = wb['profile_params']['corner_radius']/uc
            if callable(wb['profile']):
                f = cur_wb['profile']
                wb['profile'] = lambda t, f=f: f(t)/uc
            elif isinstance(wb['profile'], (list, tuple, np.ndarray)):
                arr = np.array(wb['profile'], dtype=float)
                wb['profile'] = np.c_[arr[:,0], arr[:,1]/uc]
            #
            cls = wbb.classify(wb)
            candidates = wbb.check_footprint_overlap(wb, cls)
            dim_tags = wbb.build_airbridge(wb, cls) if wb.get('kind') == 'airbridge' else wbb.build(wb, cls)
            self._gmsh_sync()
            name = 'wirebond::' + wb['name']   #Namespaced to avoid clashing with the other entity keys (metal_0, x_neg, ports...)
            assert name not in dict_all_entities, f"Wire-bond name '{wb['name']}' is used more than once."
            for other in candidates:
                GMSH_WireBond_Builder.check_3D_intersection(wb['name'], dim_tags, other, dict_all_entities['wirebond::' + other])
            wb_names.append(name)
            dict_all_entities[name] = dim_tags
            #
            fine_mesh = None
            if wb['fine_mesh']:
                fine_mesh = {k: wb[k]/uc for k in mesh_keys}
            wb_meta[name] = {
                'name': wb['name'],
                'kind': wb.get('kind', 'wirebond'),
                'classification': cls,
                'category': cls['category'],
                'bond_type': cls['bond_type'],
                'boundary': wb['boundary'],
                'material': wb['material_conductor'],
                'epr_metal_air': wb['epr_metal_air'],
                'fine_mesh': fine_mesh,
                'is_3D': any(x[0] == 3 for x in dim_tags),
                'pos_start': wb['pos_start'], 'pos_end': wb['pos_end'], 'width': wb['width']
            }
        return wb_names, wb_meta

    def _merge_wirebonded_conductors(self, contiguous_metal_mapping, metals, wirebond_attrs):
        '''
        Union-find over the metallic conductors joined by wire-bonds. Returns the merged contiguous_metal_mapping (lists of
        physical groups per conductor, including the wire-bond groups), the merged shapely conductors (in the same order)
        and the list of original metal indices forming each conductor. Each wire-bond metadata dictionary gets the
        'conductor_index' key.
        '''
        num_metals = len(metals)
        parent = list(range(num_metals))
        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i
        for meta in wirebond_attrs:
            idx = [f['metal_index'] for f in meta['classification']['feet'] if f['metal_index'] is not None]
            for i in idx[1:]:
                ra, rb = find(idx[0]), find(i)
                if ra != rb:
                    parent[max(ra, rb)] = min(ra, rb)
        roots = []
        for m in range(num_metals):
            if find(m) not in roots:
                roots.append(find(m))   #Ordered by the lowest original metal index in each conductor
        conductor_groups = [[m for m in range(num_metals) if find(m) == r] for r in roots]
        new_mapping = [sum([contiguous_metal_mapping[m] for m in grp], []) for grp in conductor_groups]
        new_polys = [metals[grp[0]] if len(grp) == 1 else shapely.unary_union([metals[m] for m in grp]) for grp in conductor_groups]
        for meta in wirebond_attrs:
            idx = [f['metal_index'] for f in meta['classification']['feet'] if f['metal_index'] is not None]
            if len(idx) > 0:
                c = roots.index(find(idx[0]))
            else:
                #A bond touching no conductor (dielectric foot to off-chip) is its own floating conductor
                new_mapping.append([])
                new_polys.append(shapely.LineString([meta['pos_start'], meta['pos_end']]).buffer(meta['width']/2))
                conductor_groups.append([])
                c = len(new_mapping) - 1
            new_mapping[c].append(meta['physical_group'])
            meta['conductor_index'] = c
        return new_mapping, new_polys, conductor_groups

    def _split_wirebond_2D_parts(self, dict_all_entities, wb_names):
        names_2D = []
        for x in wb_names:
            dict_all_entities[x + '__2D'] = [y for y in dict_all_entities[x] if y[0] == 2]
            dict_all_entities[x] = [y for y in dict_all_entities[x] if y[0] == 3]
            names_2D.append(x + '__2D')
        return names_2D

    def _merge_wirebond_2D_parts(self, dict_all_entities, wb_names):
        for x in wb_names:
            dict_all_entities[x] = sorted(dict_all_entities[x] + dict_all_entities.pop(x + '__2D'), key=lambda tup: tup[0])

    def _filter_wirebond_volumes_from_airbox(self, dict_all_entities, wb_names):
        if len(wb_names) == 0:
            return
        wb_vols = set()
        for x in wb_names:
            wb_vols.update([y for y in dict_all_entities[x] if y[0] == 3])
        dict_all_entities['airbox'] = [x for x in dict_all_entities['airbox'] if x not in wb_vols]

    def _finalise_wirebond_surfaces(self, dict_all_entities, wb_names, metal_names, gap_names, sim_construct_names, wb_meta):
        '''
        Replaces each wire-bond volume by its boundary surfaces and deletes the volume (leaving a hole in the air-box).
        Surfaces shared with the planar metals (i.e. the landing faces) stay with the metals, while surfaces on the
        air-box walls (e.g. bonds terminating on the enclosure floor) stay with the far-field boundaries.
        '''
        self._wb_wall_faces = set()
        if len(wb_names) == 0:
            return
        metal_surfs = set()
        for x in metal_names:
            metal_surfs.update([y[1] for y in dict_all_entities[x] if y[0] == 2])
        port_surfs = set()
        for x in sim_construct_names:
            port_surfs.update([y[1] for y in dict_all_entities[x] if y[0] == 2])
        all_wb_surfs = set()
        for wb_name in wb_names:
            vols = [x[1] for x in dict_all_entities[wb_name] if x[0] == 3]
            surfs = [x[1] for x in dict_all_entities[wb_name] if x[0] == 2]
            for v in vols:
                surfs += gmsh.model.getAdjacencies(3, v)[1].tolist()
            if len(vols) > 0:
                gmsh.model.occ.remove([(3,v) for v in vols], recursive=False)
            final = []
            on_port = False
            wb_meta[wb_name]['grounds_to_enclosure'] = False
            for s in dict.fromkeys(surfs):
                if s in metal_surfs:
                    continue
                if s in port_surfs:
                    on_port = True
                    continue
                com = gmsh.model.occ.getCenterOfMass(2, s)
                if np.min(np.abs(self._extents.T - np.array(com))) <= 1e-6:
                    self._wb_wall_faces.add(s)  #Lies on the air-box wall inside the hollowed bond - i.e. not part of the domain
                    wb_meta[wb_name]['grounds_to_enclosure'] = True
                    continue
                final.append((2,s))
            if on_port:
                print(f"WARNING: Wire-bond {wb_meta[wb_name]['name']} lands on a lumped port/junction surface; that surface is left to the port.")
            dict_all_entities[wb_name] = final
            all_wb_surfs.update([x[1] for x in final])
        self._gmsh_sync()
        #A foot overhanging a dielectric gap leaves its bottom face coincident with the gap surface - it belongs to the bond
        for gap_name in gap_names:
            dict_all_entities[gap_name] = [x for x in dict_all_entities[gap_name] if not (x[0] == 2 and x[1] in all_wb_surfs)]

    def _conv_to_lists(self, data):
        if isinstance(data, list):
            return data
        else:
            return [data]

    def _divide_surface_with_area(self, new_area, metals:list, dielectric_gaps:list, areas: list):
        #Try to keep the contiguous metals contiguous - hence metals is a list of dictionaries; e.g. [{'metal0':poly1, 'metal1':poly2}, {'metal3':poly3}] means poly1/poly2 comprise the first metal while
        #poly3 is another metallic piece. Same with the gaps. The list `areas` is a list of lists; e.g. [['metal0i', 'metal1ii'], [metal2i]] implies that the first area comprises of the pieces given
        #by those IDs...
        final_new_area = []
        final_areas = areas
        final_metals = []
        for cur_metal in metals:
            new_metal = {}
            for cur_piece in cur_metal:
                intersecting_piece = shapely.intersection(cur_metal[cur_piece], new_area)
                left_over_piece = shapely.difference(cur_metal[cur_piece], new_area)
                left_over = False
                if left_over_piece.area > 1e-15:
                    #An edge-case where the intersection splits the polygon into multiple separated pieces...
                    left_over_piece = ShapelyEx.shapely_to_list(left_over_piece)
                    if len(left_over_piece) == 1:
                        new_metal[cur_piece] = left_over_piece[0]
                    else:
                        for m in range(len(left_over_piece)):
                            new_metal[cur_piece+f'_{m}'] = left_over_piece[m]
                    left_over = True
                if intersecting_piece.area > 1e-15:
                    new_name = cur_piece+'i'
                    new_metal[new_name] = intersecting_piece
                    final_new_area.append(new_name)
                    for m in range(len(areas)):
                        if cur_piece in areas[m]:
                            if not left_over:
                                areas[m].pop(cur_piece)
                            areas[m].append(new_name)
            final_metals.append(new_metal)
        #
        final_gaps = []
        for cur_gap in dielectric_gaps:
            new_gap = {}
            for cur_piece in cur_gap:
                intersecting_piece = shapely.intersection(cur_gap[cur_piece], new_area)
                left_over_piece = shapely.difference(cur_gap[cur_piece], new_area)
                left_over = False
                if left_over_piece.area > 1e-15:
                    #An edge-case where the intersection splits the polygon into multiple separated pieces...
                    left_over_piece = ShapelyEx.shapely_to_list(left_over_piece)
                    if len(left_over_piece) == 1:
                        new_gap[cur_piece] = left_over_piece[0]
                    else:
                        for m in range(len(left_over_piece)):
                            new_gap[cur_piece+f'_{m}'] = left_over_piece[m]
                    left_over = True
                if intersecting_piece.area > 1e-15:
                    new_name = cur_piece+'i'
                    new_gap[new_name] = intersecting_piece
                    final_new_area.append(new_name)
                    for m in range(len(areas)):
                        if cur_piece in areas[m]:
                            if not left_over:
                                areas[m].pop(cur_piece)
                            areas[m].append(new_name)
            final_gaps.append(new_gap)
        #
        areas.append(final_new_area)
        return final_metals, final_gaps, areas


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

        # See here: https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/tutorials/python/t16.py for details on the mapping...
        main_objs, leMap = gmsh.model.occ.fragment(main_dim_tags, tool_dim_tags, removeObject=True, removeTool=True)

        # Using leMap (which gives the list of dim,tag pairs that the original dim,tag pairs in main_objects map onto), redo the main_objects and
        # tool_objects while making sure to keep the appropriate names:
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

        # Using leMap (which gives the list of dim,tag pairs that the original dim,tag pairs in target_dim_tags map onto), redo the target target_objects
        # making sure to keep the appropriate names
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
        assert (
            np.min(coords[:, 0]) >= self._extents[0][0]-1e-9
            and np.max(coords[:, 0]) <= self._extents[0][1]+1e-9
        ), f"Element {message} leaks outside x-axis boundary: ({np.min(coords[:,0])} >= {self._extents[0][0]})."
        assert (
            np.min(coords[:, 1]) >= self._extents[1][0]-1e-9
            and np.max(coords[:, 1]) <= self._extents[1][1]+1e-9
        ), f"Element {message} leaks outside y-axis boundary: ({np.min(coords[:,1])} >= {self._extents[1][0]})."
        assert (
            np.min(coords[:, 2]) >= self._extents[2][0]-1e-9
            and np.max(coords[:, 2]) <= self._extents[2][1]+1e-9
        ), f"Element {message} leaks outside z-axis boundary: ({np.min(coords[:,2])} >= {self._extents[2][0]})."

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

        # define list to store lines
        lines = []

        # create 2D points in gmsh
        if coords.shape[1] == 2:
            points = [gmsh.model.occ.addPoint(coord[0], coord[1], self._geom_processor.chip_centre[2]) for coord in coords]
        else:
            points = [gmsh.model.occ.addPoint(coord[0], coord[1], coord[2]) for coord in coords]

        # update model with the created geometry items
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()

        # draw lines from points
        for j,value in enumerate(points):
            if(j<len(points)-1):
                # skip lines on problematic points
                try:
                    line = gmsh.model.occ.add_line(points[j], points[j+1])
                    lines.append(line)
                except:
                    print(f"GMSH: skipped line ({coords[j]}, {coords[j+1]})")

                # # original
                # line = gmsh.model.occ.add_line(points[j], points[j+1])
                # lines.append(line)

        line = gmsh.model.occ.add_line(points[len(points)-1], points[0])
        lines.append(line)      

        # create curved loop
        curve_loop = gmsh.model.occ.add_curve_loop(lines)

        # create_surface
        surface = gmsh.model.occ.add_plane_surface([curve_loop])

        # update model with the created geometry items
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

        # define region of device for simulation

        # surface of chip base rests at z = 0
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
            zMin = self.center_z - np.abs(self.size_z) - boundary_distances['z_prop_bottom']*np.abs(self.size_z)
            zMax = self.center_z + boundary_distances['z_prop_top']*np.abs(self.size_z)
        else:
            zMin = self.center_z - np.abs(self.size_z) - boundary_distances['z_neg']/self._unit_conv
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
