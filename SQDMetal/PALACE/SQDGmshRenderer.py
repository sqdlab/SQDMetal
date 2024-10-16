from qiskit_metal.qgeometries.qgeometries_handler import QGeometryTables
from qiskit_metal.renderers.renderer_mpl.mpl_renderer import QMplRenderer
import gmsh
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import shapely
import qiskit_metal 
from qiskit_metal import designs, draw
from qiskit_metal import MetalGUI, Dict, open_docs
from qiskit_metal import qgeometries
from qiskit_metal.toolbox_metal import math_and_overrides
from qiskit_metal.qlibrary.core import QComponent
from collections import OrderedDict
from SQDMetal.Utilities.QUtilities import QUtilities
from SQDMetal.Utilities.QiskitShapelyRenderer import QiskitShapelyRenderer
from SQDMetal.Utilities.ShapelyEx import ShapelyEx

class Palace_Gmsh_Renderer:
    #meshing parameter
    lc = None

    ###GeoPandas GeoSeries###
    gpd_polys = []

    ###Lists for Gmsh Geometries###
    gmsh_metal_surface = None
    gmsh_dielectric_gaps = []
    gmsh_gap_list = []

    #constructor takes qiskit-metal design
    def __init__(self, design):
        self.design = design

        #reset gmsh model
        # gmsh.initialize()
        # gmsh.model.remove()
        # gmsh.finalize()

        #TODO use SQDmetal utilities parse value
        #get dimensions of chip base and convert to design units in 'mm'
        self.center_x = self.design.parse_value(self.design.chips['main'].size.center_x)
        self.center_y = self.design.parse_value(self.design.chips['main'].size.center_y)
        self.center_z = self.design.parse_value(self.design.chips['main'].size.center_z)
        self.size_x = self.design.parse_value(self.design.chips['main'].size.size_x)
        self.size_y = self.design.parse_value(self.design.chips['main'].size.size_y)
        self.size_z = self.design.parse_value(self.design.chips['main'].size.size_z)


    def _prepare_design(self, metallic_layers, ground_plane, sim_constructs, fillet_resolution, simulation, **kwargs):
        '''Takes qiskit metal design geometries from shapely and turns all geometries into polygons. Some elements such as
        launch pads are already polygons however other elements such as meander resonators are linestrings and need to be
        converted to polygons in shapely. 
        
        Qiskit metal usually has two polygons associated with an element if it is not a linestring geometry. These are
        (1) metal and (2) a dielectric gap surrounding element. This method also separates the metal and dielectric gap components
        and stores them in their respective lists. '''

        fuse_threshold = kwargs.get('fuse_threshold', 1e-12)

        #Remove any previous model and add new gmsh model
        gmsh.model.remove()
        gmsh.finalize()
        gmsh.initialize()
        gmsh.model.add('qiskit_to_gmsh')

        unit_conv = 1   #TODO: Because it is stuck in mm?!

        #Handle the ground plane...
        qmpl = QiskitShapelyRenderer(None, self.design, None)
        gsdf = qmpl.get_net_coordinates(fillet_resolution)
        filt = gsdf.loc[gsdf['subtract'] == True]
        space_polys = ShapelyEx.fuse_polygons_threshold(filt['geometry'].buffer(0), fuse_threshold)
        if isinstance(space_polys, shapely.geometry.multipolygon.MultiPolygon):
            space_polys = [x for x in space_polys.geoms]
        else:
            space_polys = [space_polys] #i.e. it's just a lonely Polygon object...
        space_polys = [shapely.affinity.scale(x, xfact=unit_conv, yfact=unit_conv, origin=(0,0)) for x in space_polys]
        #cutout dielectric gaps from metal ground plain <-- haha ground "plain"...
        #create chip surface where structures will sit
        metal_surface = shapely.geometry.box(self.center_x - 0.5*self.size_x, self.center_y - 0.5*self.size_y,
                             self.center_x + 0.5*self.size_x, self.center_y + 0.5*self.size_y)
        ground_plain = shapely.difference(metal_surface, shapely.geometry.multipolygon.MultiPolygon(space_polys))
        #
        #Gather all polygons into contiguous groups
        metal_polys = {}
        num_polys = 0
        if not ground_plain.is_empty:
            metal_polys[0] = ground_plain
            num_polys += 1
        for cur_layer in metallic_layers:
            if cur_layer['type'] == 'design_layer':
                cur_layer['unit_conv'] = unit_conv
                cur_layer['resolution'] = fillet_resolution
                metal_polys_all, metal_sel_ids = QUtilities.get_metals_in_layer(self.design, **cur_layer)
                unique_ids = np.unique(metal_sel_ids)
                for m in range(unique_ids.size):
                    metal_polys[m+num_polys] = []
                for m in range(len(metal_polys_all)):
                    metal_polys[metal_sel_ids[m]+num_polys] += [metal_polys_all[m]]
                for m in range(unique_ids.size):
                    metal_polys[m+num_polys] = shapely.geometry.multipolygon.MultiPolygon(metal_polys[m+num_polys])
                num_polys += unique_ids.size
            elif cur_layer['type'] == 'Uclip':
                if cur_layer['clip_type'] == 'inplane':
                    Uclip = QUtilities.get_RFport_CPW_groundU_Launcher_inplane(self.design, cur_layer['qObjName'], cur_layer['thickness_side'], cur_layer['thickness_back'], cur_layer['separation_gap'], cur_layer['unit_conv_extra']*1e3) #TODO: Stuck in mm?!
                metal_polys[num_polys] = shapely.Polygon(Uclip)
                num_polys += 1
        #Fuse all contiguous polygons into single groups
        def get_poly_cardinality(poly):
            if isinstance(poly, shapely.geometry.multipolygon.MultiPolygon):
                return len(poly.geoms)
            else:
                return 1
        new_polys = {}
        new_polys[0] = metal_polys[0]
        for m in range(1, num_polys):
            found_match = False
            for n in range(len(new_polys)):
                cur_check = ShapelyEx.fuse_polygons_threshold([metal_polys[m], new_polys[n]], fuse_threshold)
                if get_poly_cardinality(cur_check) < get_poly_cardinality(metal_polys[m]) + get_poly_cardinality(new_polys[n]):
                    new_polys[n] = cur_check
                    found_match = True
                    break
            if not found_match:
                new_polys[len(new_polys)] = metal_polys[m]

        return self._sim_prep(simulation.replace('_',' '), new_polys, sim_constructs)

    def _sim_prep(self, type_name, metal_polys, sim_constructs):
        '''Prepare simulation by fusing dielectric gaps and then cutting these out from the chip's metal surface.
            All metal components are fused including the ground plane. At this stage all geometries are still in
            shapely.'''

        metal_list = []
        for m in metal_polys:
            cur_poly = metal_polys[m]
            if isinstance(cur_poly, shapely.geometry.multipolygon.MultiPolygon):
                lePolys = list(cur_poly.geoms)
            else:
                lePolys = [cur_poly]
            for cur_poly in lePolys:

                p = gpd.GeoSeries(cur_poly)
                p.plot()
                metal_simplified = cur_poly.simplify(1e-6)
                gmsh_exterior = self._draw_polygon_from_coords(metal_simplified.exterior.coords[:])
                if len(metal_simplified.interiors) >= 1:
                    interiors = []
                    for _,interior in enumerate(metal_simplified.interiors):
                        gmsh_interior = self._draw_polygon_from_coords(interior.coords[:])
                        interiors.append((2,gmsh_interior))
                    gmsh_surface, gmsh_surface_map = gmsh.model.occ.cut([(2,gmsh_exterior)], 
                                                                        interiors,
                                                                        removeObject=True, removeTool=True)
                    metal_list.append((m, 2,gmsh_surface[0][1]))
                else:
                    metal_list.append((m, 2,gmsh_exterior))

        #update geometries
        gmsh.model.geo.synchronize()
        gmsh.model.occ.synchronize()

        #create volume of dielectric base
        dielectric_box = self._draw_chip_base()

        #update geometries
        gmsh.model.geo.synchronize()
        gmsh.model.occ.synchronize()

        #add physical groups for metal list
        leMetals = []
        for m in metal_polys:
            cur_metals = [x[2] for x in metal_list if x[0] == m]
            leMetals.append( gmsh.model.addPhysicalGroup(2, cur_metals, name = 'metal'+str(m)) )
        

        #create list to fragment with chip base
        fragment_list = [(x[1], x[2]) for x in metal_list]

        #Add in polygons for simulation constructs...
        lePortPolys = {}
        for cur_sim_poly in sim_constructs:
            sim_poly = self._draw_polygon_from_coords(np.array(cur_sim_poly[1]))
            lePortPoly = gmsh.model.addPhysicalGroup(2, [sim_poly], name = cur_sim_poly[0])
            fragment_list.append((2, sim_poly))
            lePortPolys[cur_sim_poly[0]] = lePortPoly

        gmsh.model.geo.synchronize()
        gmsh.model.occ.synchronize()

        #fragment the newly created elements with the dielectric volume
        chip, chip_map = gmsh.model.occ.fragment([dielectric_box[1]], fragment_list, removeObject=True, removeTool=True)
        
        #update geometries
        gmsh.model.geo.synchronize()
        gmsh.model.occ.synchronize()

        #add dielectric volume as dielectric base
        leDielectric = gmsh.model.addPhysicalGroup(3, [chip[0][1]], name = 'dielectric_base')
        Palace_Gmsh_Renderer.frag_dielectric_vol = chip[0][1]

        # #add dielectric gaps
        # start_index = len(Palace_Gmsh_Renderer.gmsh_metals) + 1 + 4 
        # end_index = len(chip) - 1
        # gaps = chip[start_index:end_index]

        # Palace_Gmsh_Renderer.gmsh_gap_list = []
        # for _,gap in enumerate(gaps):
        #     Palace_Gmsh_Renderer.gmsh_gap_list.append(gap[1])
        
        # gmsh.model.addPhysicalGroup(2, Palace_Gmsh_Renderer.gmsh_gap_list, name = 'dielectric_gaps')

        #draw airbox surrounding chip
        leFarField, leAirBox = self._draw_air_box()

        #get physical group indentifiers to use in config file
        ret_dict = {
            'air_box'   : leAirBox,
            'metals'    : leMetals,
            'far_field' : leFarField,
            'ports'     : lePortPolys,
            'dielectric': leDielectric
        }

        # #plot fused metal polygons and add fused metal polygons to list
        # fig2, ax2 = plt.subplots()  # a figure with a single Axes
        # for i, geom in enumerate(fused_metals.geoms):
        #     Palace_Gmsh_Renderer.fused_metals.append(geom)
        #     poly = gpd.GeoSeries([geom])
        #     Palace_Gmsh_Renderer.gpd_polys.append(poly)
        #     poly.plot(ax=ax2)
        
        # #add title and axes for plot
        # plt.title(f'Geometry for {type_name} Simulation')
        # plt.xlabel('x (mm)')
        # plt.ylabel('y (mm)')
        # plt.show()

        return ret_dict
         

    def _draw_components(self, component_list):
        '''Takes a list of Comoponents and draws them based on their type: path, polygon or junction'''
        
        for component_name in component_list:
            if(not self.design.qgeometry.get_component(component_name)['path'].empty):
                path_to_draw = self.design.qgeometry.get_component(component_name)['path']
                for index in path_to_draw.index:
                    self.draw_path(path_to_draw.loc[index], component_name, index, flag = 'path')
            if(not self.design.qgeometry.get_component(component_name)['poly'].empty):
                poly_to_draw = self.design.qgeometry.get_component(component_name)['poly']
                for index in poly_to_draw.index:
                    self._draw_polygon(poly_to_draw.loc[index], component_name, index, flag = 'poly')
            if(not self.design.qgeometry.get_component(component_name)['junction'].empty):
                junc_to_draw = self.design.qgeometry.get_component(component_name)['junction']
                for index in junc_to_draw.index:
                    self.draw_path(junc_to_draw.loc[index], component_name, index, flag = 'junction')


    def _draw_chip_base(self):
        '''This method draws the chip base into gmsh given the dimensions defined by the user in qiskit metal'''

        #half values of the sizes
        half_size_x = self.size_x/2
        half_size_y = self.size_y/2

        #store coordinates for the top surface of the chip
        surface_1 = {'point1': [self.center_x + half_size_x, self.center_y + half_size_y, self.center_z],
                     'point2': [self.center_x + half_size_x, self.center_y - half_size_y, self.center_z],
                     'point3': [self.center_x - half_size_x, self.center_y - half_size_y, self.center_z],
                     'point4': [self.center_x - half_size_x, self.center_y + half_size_y, self.center_z]}
        
        #define lists to store points, lines and surfaces
        points = []
        lines = []
        surfaces = []

        #add points for chip base
        for i,value in enumerate(surface_1):
            point = gmsh.model.occ.addPoint(surface_1[value][0], surface_1[value][1], surface_1[value][2])
            points.append(point)

        #draw lines for chip base
        for j,value in enumerate(points):
            if(j<len(points)-1):
                line = gmsh.model.occ.add_line(points[j], points[j+1])
                lines.append(line)        
        line = gmsh.model.occ.add_line(points[len(points)-1], points[0])
        lines.append(line)        

        #create curved loop
        curve_loop = gmsh.model.occ.add_curve_loop(lines)
        
        #create_surface
        base_surface = gmsh.model.occ.add_plane_surface([curve_loop])

        #create volume from top surface using extrude function in gmsh
        dielectric_box = gmsh.model.occ.extrude([(2, base_surface)],0,0,self.size_z)

        #update model with the created geometry items
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()

        return dielectric_box

        
    def _draw_polygon(self, polygon):
        '''takes a shapely polygon object or pandas series 
            in as an argument and then draws it in Gmsh and returns the surface ID''' 

        #simplify polygon, remove points that are too close together
        new_polygon = polygon.simplify(1e-6)

        #gets polygon coordinates depending on the type of argument fed to the function
        x_coords = new_polygon.exterior.coords.xy[0]
        y_coords = new_polygon.exterior.coords.xy[1]
       
        #define lists to store points and lines
        points = []
        lines = []
    
        #create 2D points in gmsh
        for i,coord in enumerate(x_coords):
            point = gmsh.model.occ.addPoint(x_coords[i], y_coords[i], 
                                            self.design.parse_value(self.design.chips['main'].size.center_z))
            points.append(point)

        points.pop(-1) #remove last point as it is repeated

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


    def _draw_polygon_from_coords(self, coords):
        '''draws polygon in GMSH from x and y coordinates and returns the surface ID''' 

        #gets polygon coordinates depending on the type of argument fed to the function
        x_coords = []
        y_coords = []

        for _,coordintes in enumerate(coords):
            x_coords.append(coordintes[0])
            y_coords.append(coordintes[1])

        #define lists to store points and lines
        points = []
        lines = []
    
        #create 2D points in gmsh
        for i,coord in enumerate(x_coords):
            point = gmsh.model.occ.addPoint(x_coords[i], y_coords[i], 
                                            self.design.parse_value(self.design.chips['main'].size.center_z))
            points.append(point)

        points.pop(-1) #remove last point as it is repeated

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


    def _draw_path(self, path: pd.Series):
        '''Takes a pandas series in as an argument and then draws the path as a polygon in shapely. Returns the polygon generated 
        from the path'''

        #get width and buffer amount for path
        width = path.width
        buffer_amt = width/2

        #fillet path using QMplRenderer
        qmpl = QMplRenderer(None, self.design, None)
        path_filleted = qmpl.fillet_path(path)

        #buffer the path by the user defined width
        poly_path = shapely.buffer(path_filleted, distance = buffer_amt, cap_style = 'flat')

        return poly_path


    def _draw_air_box(self):
        '''draw airbox into gmsh'''

        #for air box choose increase in dimensions in 'mm'
        air_box_delta_x = (1/4) * self.size_x
        air_box_delta_y = (1/4) * self.size_y
        air_box_delta_z = 2 * self.size_z

        #half values of the sizes plus add increase for air box
        half_size_x = self.size_x/2 + air_box_delta_x/2
        half_size_y = self.size_y/2 + air_box_delta_y/2

        #store coordinates for surfaces of chip
        air_box_surface = {'point1': [self.center_x + half_size_x, self.center_y + half_size_y, (self.center_z - self.size_z)],
                            'point2': [self.center_x + half_size_x, self.center_y - half_size_y, (self.center_z - self.size_z)],
                            'point3': [self.center_x - half_size_x, self.center_y - half_size_y, (self.center_z - self.size_z)],
                            'point4': [self.center_x - half_size_x, self.center_y + half_size_y, (self.center_z - self.size_z)]}
        
        #define lists to store points, lines and surfaces
        points = []
        lines = []
        surfaces = []

        #add points for air_box
        for i,value in enumerate(air_box_surface):
            point = gmsh.model.occ.addPoint(air_box_surface[value][0], air_box_surface[value][1], 
                                    air_box_surface[value][2])
            points.append(point)

        #draw lines for airbox
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

        #create volume from top surface of airbox using extrude function in gmsh
        air_box = gmsh.model.occ.extrude([(2, surface)],0,0,air_box_delta_z)

        #cut out chip base from airbox
        air_box_cutout, air_box_cutout_map = gmsh.model.occ.fragment([(3, Palace_Gmsh_Renderer.frag_dielectric_vol)], [air_box[1]],  
                                                 removeObject=True, removeTool=True)

        #update model with the created geometry items
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()

        #add physical group for dielectric chip base
        airBox = gmsh.model.addPhysicalGroup(3, [air_box_cutout[1][1]], name = 'air_box')

        #add physical group for the surfaces of the air box which will represent the far field bondary conditions
        far_field = []
        for i, value in enumerate(gmsh.model.getEntities(2)[-7:]): #last 7 surface entities represnt the far field
            far_field.append(value[1])
        return gmsh.model.addPhysicalGroup(2, far_field, name = 'far_field'), airBox


    def view_design_components(self):
        '''view the wireframe with the physical groups of the design displayed'''

        gmsh.option.setNumber('Geometry.CurveWidth',0.1)
        gmsh.option.setNumber('Geometry.LabelType', 4)
        gmsh.option.setNumber('Geometry.Points', 0)
        gmsh.option.setNumber('Geometry.SurfaceLabels', 1)
        gmsh.option.setNumber('General.GraphicsFontSize', 14)
        gmsh.option.setColor('Geometry.Color.Surfaces',178,37,23, 190)
        gmsh.fltk.run()
    

    def print_physical_groups(self):
        '''print out the physical groups which are the boundary conditions for the simulation'''

        physical_groups = gmsh.model.get_physical_groups()
        print('\nPhysical Groups')
        for i, value in enumerate(physical_groups):
            phys_group_name = gmsh.model.get_physical_name(value[0], value[1])
            print('name:', phys_group_name + '\t\t', 'identifier:', value[1])
        

    def _box_mesh_geometry(self):

        lc = 100e-3

        #get the coordinates for each object in the design
        mesh_fields = []       
        for i,polygon in enumerate(Palace_Gmsh_Renderer.fused_metals):
            mesh_fields.append((i+1))
            gmsh.model.mesh.field.add("Box", i+1)
            gmsh.model.mesh.field.setNumber(i+1, "VIn", lc / 10)
            gmsh.model.mesh.field.setNumber(i+1, "VOut", lc)
            #bounds returns ((minx, miny, maxx, maxy))
            gmsh.model.mesh.field.setNumber(i+1, "XMin", polygon.bounds[0])
            gmsh.model.mesh.field.setNumber(i+1, "XMax", polygon.bounds[2])
            gmsh.model.mesh.field.setNumber(i+1, "YMin", polygon.bounds[1])
            gmsh.model.mesh.field.setNumber(i+1, "YMax", polygon.bounds[3])
            gmsh.model.mesh.field.setNumber(i+1, "Thickness", lc / 10)
        
        #create min field and set min field as background mesh
        field_no = len(mesh_fields) + 1
        gmsh.model.mesh.field.add("Min", field_no)
        gmsh.model.mesh.field.setNumbers(field_no, "FieldsList", mesh_fields)
        gmsh.model.mesh.field.setAsBackgroundMesh(field_no)

        #create mesh
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
        gmsh.model.mesh.generate(dim=3)
        gmsh.fltk.run()


    def _create_metals_for_fine_mesh(self):
        '''draw metal surface into gmsh in order to mesh finer in these regions. Note that these surfaces are just used to tell gmsh
        where to mesh finer, they are not included in the output of the design.'''

        #get metal polygons
        metal_fuse_list = []
        for i, metal_poly in enumerate(Palace_Gmsh_Renderer.metal_polygons):
            #mesh all metals finely except launchpads
            if (type(self.design.components.__getitem__(metal_poly[0])) != qiskit_metal.qlibrary.terminations.launchpad_wb.LaunchpadWirebond):
                metal_fuse_list.append(metal_poly[1])
        
        #Fuse the metal elements using Qutilities function (unary_union)
        fused_metals = ShapelyEx.fuse_polygons_threshold(metal_fuse_list)

        #list to store geopandas polygons
        mesh_gpd_polys = []

        # #plot fused metal polygons and add fused metal polygons to list
        # fig2, ax2 = plt.subplots()  # a figure with a single Axes
        # for _, geom in enumerate(fused_metals.geoms):
        #     poly = gpd.GeoSeries([geom])
        #     mesh_gpd_polys.append(poly)
        #     poly.plot(ax=ax2)
        
        #add title and axes for plot
        plt.title('Metal Surfaces for Fine Meshing')
        plt.xlabel('x (mm)')
        plt.ylabel('y (mm)')
        
        #create metal surface in Gmsh
        gmsh_metals_for_meshing = []
        for i, metal in enumerate(mesh_gpd_polys):
            metal_simplified = metal.simplify(1e-6)
            gmsh_poly = self._draw_polygon_from_coords(metal_simplified.exterior[0].coords[:])
            gmsh_metals_for_meshing.append(gmsh_poly)
        
        return gmsh_metals_for_meshing

    def _set_field_params(self, field_name, index, **kwargs):
        gmsh.model.mesh.field.add(field_name, index)
        for cur_key in kwargs:
            cur_val = kwargs.get(cur_key)
            if isinstance(cur_val, (list, tuple)):
                gmsh.model.mesh.field.setNumbers(index, cur_key, cur_val)
            else:
                gmsh.model.mesh.field.setNumber(index, cur_key, cur_val)
    
    def _set_field_fine_mesh_path(self, index, path, mesh_sampling, min_size, max_size):
        #Path given as (x,y) coordinates
        lePoints = [gmsh.model.geo.addPoint(x[0],x[1], 0) for x in path]
        leLines = [gmsh.model.geo.addLine(lePoints[m-1], lePoints[m]) for m in range(1, len(lePoints))]
        gmsh.model.geo.synchronize()
        self._set_field_params("Distance", index, Sampling=mesh_sampling, CurvesList=leLines)
        self._set_field_params("Threshold", index+1, InField=index, SizeMin=min_size, SizeMax=max_size, DistMin=max_size/4, DistMax=3*max_size)
        return index+1
    
    def _set_field_fine_mesh_box(self, index, x_bnds, y_bnds, mesh_sampling, min_size, max_size):
        self._set_field_params("Box", index, XMin=min(x_bnds), XMax=max(x_bnds), YMin=min(y_bnds), YMax=max(y_bnds), #ZMin=-0.01, ZMax=0.01,
                                             Thickness = 3*max_size, VIn=min_size, VOut=max_size)
        return index

    def _intelligent_mesh(self, simulation_type, min_size, max_size, mesh_sampling):
        '''mesh chip geometry intelligently by meshing finer in regions where the geometry is smaller.'''

        #turn off meshing parameters that are not required
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

        #set mesh algorithm Mesh.Algorithm3D to HXT (option - 10) or Delaunay 3D (option - 1)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1) 

        #define distance field
        gmsh.model.mesh.field.add("Distance", 1)
        gmsh.model.mesh.field.setNumber(1, "Sampling", mesh_sampling)

        if simulation_type == 'capacitance_simulation':
            gmsh.model.mesh.field.setNumbers(1, "SurfacesList", self._create_metals_for_fine_mesh()) #Palace_Gmsh_Renderer.gmsh_metals
        if simulation_type == 'eigenmode_simulation' or 'driven_simulation':
            gmsh.model.mesh.field.setNumbers(1, "SurfacesList", [2])#self._create_metals_for_fine_mesh()) #Palace_Gmsh_Renderer.gmsh_gap_list

        #define mesh field
        self._set_field_params("Threshold", 2, InField=1, SizeMin=min_size, SizeMax=max_size, DistMin=max_size/4, DistMax=3*max_size)

        self._set_field_params("Box", 3, XMin=-0.3, XMax=0.3, YMin=-1.5, YMax=0.0, ZMin=-0.01, ZMax=0.01, VIn=0.001, VOut=0.1)
        self._set_field_params("Threshold", 4, InField=3, SizeMin=min_size, SizeMax=max_size, DistMin=max_size/4, DistMax=3*max_size)

        self._set_field_params("Min", 5, FieldsList=[2, 4])

        self._set_field_fine_mesh_path(6, [(-0.3,-0.5), (0.3,-1.0), (-0.3,-1.5)], mesh_sampling, min_size, max_size)

        #set background field and render mesh
        gmsh.model.mesh.field.setAsBackgroundMesh(7)
        gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
        gmsh.model.mesh.generate(dim = 3)
    
    def fine_mesh(self, fine_mesh_params):
        '''mesh chip geometry intelligently by meshing finer in regions where the geometry is smaller.'''

        #turn off meshing parameters that are not required
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

        #set mesh algorithm Mesh.Algorithm3D to HXT (option - 10) or Delaunay 3D (option - 1)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1) 

        cur_field_index = 1
        global_field_inds = []
        for cur_fine in fine_mesh_params:
            if cur_fine['type'] == 'path':
                new_ind = self._set_field_fine_mesh_path(cur_field_index, cur_fine['path'], cur_fine['mesh_sampling'], cur_fine['min_size'], cur_fine['max_size'])
            if cur_fine['type'] == 'box':
                new_ind = self._set_field_fine_mesh_box(cur_field_index, cur_fine['x_bnds'], cur_fine['y_bnds'], cur_fine['mesh_sampling'], cur_fine['min_size'], cur_fine['max_size'])
            global_field_inds.append(new_ind)
            cur_field_index = new_ind + 1

        #Combine all field constraints
        self._set_field_params("Min", cur_field_index, FieldsList=global_field_inds)
        #set background field and render mesh
        gmsh.model.mesh.field.setAsBackgroundMesh(cur_field_index)
        gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
        gmsh.model.mesh.generate(dim = 3)





        