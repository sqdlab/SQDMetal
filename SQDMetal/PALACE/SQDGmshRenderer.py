from qiskit_metal.qgeometries.qgeometries_handler import QGeometryTables
from qiskit_metal.renderers.renderer_mpl.mpl_renderer import QMplRenderer
import gmsh
import pandas as pd
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import shapely
import qiskit_metal as metal
from qiskit_metal import designs, draw
from qiskit_metal import MetalGUI, Dict, open_docs
from qiskit_metal import qgeometries
from qiskit_metal.toolbox_metal import math_and_overrides
from qiskit_metal.qlibrary.core import QComponent
from collections import OrderedDict
from SQDMetal.Utilities.QUtilities import QUtilities


class Palace_Gmsh_Renderer:

    #lists to store entites which comprise the dielectric gaps and metals in the design
    dielectric_gap_surface = []
    metal_surfaces = []

    #class variable to store ID (int) representing the chip base
    dielectric_box = None
    
    #fragmented dielectric vol
    frag_dielectric_vol = None

    #list to store launch pads to create ports on
    launch_pads_list = []

    #lumped element ports
    ports_list = []

    #meshing parameter
    lc = None

    ###Lists for Shapely Polygons###
    metals = []
    metal_polygons = []
    fused_metals = []
    dielectric_gaps = []
    dielectric_gap_polygons = []
    fused_dielectric_gaps = []
    metal_surface = None
    ground_plane_cut = None
    ground_plane_fused_metal = None
    fused_cutouts = None

    ###GeoPandas GeoSeries###
    gpd_polys = []

    ###Lists for Gmsh Geometries###
    gmsh_metals = []
    gmsh_metal_surface = None
    gmsh_dielectric_gaps = []
    gmsh_gap_list = []

    ###Physical Groups/Boundary Conditions###
    config_metals_rf = []
    config_metals_cap = []
    config_ground_plane = None
    config_far_field = None
    config_dielectric_base = None
    config_dielectric_gaps = []
    config_air_box = None
    config_ports = []

    #constructor takes qiskit-metal design
    def __init__(self, design):
        self.design = design

        #reset gmsh model
        gmsh.model.remove()
        gmsh.finalize()

        #TODO use SQDmetal utilities parse value
        #get dimensions of chip base and convert to design units in 'mm'
        self.center_x = self.design.parse_value(self.design.chips['main'].size.center_x)
        self.center_y = self.design.parse_value(self.design.chips['main'].size.center_y)
        self.center_z = self.design.parse_value(self.design.chips['main'].size.center_z)
        self.size_x = self.design.parse_value(self.design.chips['main'].size.size_x)
        self.size_y = self.design.parse_value(self.design.chips['main'].size.size_y)
        self.size_z = self.design.parse_value(self.design.chips['main'].size.size_z)



    def _prepare_design(self, simulation):
        '''convert qiskit metal design to geometry in Gmsh'''

        #Flag to define the simulation type
        simulation_type = None

        if simulation == 'capacitance_simulation':
            simulation_type = 'cap_sim'
        elif simulation == 'eigenmode_simulation':
            simulation_type = 'eigen_sim'
        elif simulation == 'frequency_driven_simulation':
            simulation_type = 'freq_sim'
        else:
            raise Exception("Invalid simulation type entered.")

        #Start gmsh and add model
        gmsh.initialize()
        gmsh.model.add('qiskit_to_gmsh')

        #create list with component names
        component_list = self.design.all_component_names_id()
        component_names = []
        for i in range(len(component_list)):
            component_names.append(component_list[i][0])

        #qiskit metal element shapely types
        element_types = ['path', 'poly', 'junction']

        #Determine whether component is a metal or dielectric gap and add it to its corresponding list
        for i, component_name in enumerate(component_names):
            for j, element_type in enumerate(element_types):
                if(not self.design.qgeometry.get_component(component_name)[element_type].empty):
                    component_parts = self.design.qgeometry.get_component(component_name)[element_type]
                    for index in component_parts.index:
                        geom_pd_series = self.design.qgeometry.get_component(component_name)[element_type].loc[index]
                        if self.design.qgeometry.get_component(component_name)[element_type].loc[index].at['subtract'] == True:
                            Palace_Gmsh_Renderer.dielectric_gaps.append((component_name, element_type, geom_pd_series))
                        else:
                            Palace_Gmsh_Renderer.metals.append((component_name, element_type, geom_pd_series))

        #For metal geometries convert linestrings to polygons, all objects need to be polygons for plotting in gmsh
        for i, entry in enumerate(Palace_Gmsh_Renderer.metals):
            if isinstance(entry[2].geometry, shapely.geometry.linestring.LineString):
                linestring_poly =  self._draw_path(entry[2])
                Palace_Gmsh_Renderer.metal_polygons.append((entry[0], linestring_poly))
            else:
                Palace_Gmsh_Renderer.metal_polygons.append((entry[0], entry[2].geometry))

        #For dielectric gap geometries convert linestrings to polygons, all objects need to be polygons for plotting in gmsh
        for i, entry in enumerate(Palace_Gmsh_Renderer.dielectric_gaps):
            if isinstance(entry[2].geometry, shapely.geometry.linestring.LineString):
                linestring_poly =  self._draw_path(entry[2])
                Palace_Gmsh_Renderer.dielectric_gap_polygons.append((entry[0], linestring_poly))
            else:
                Palace_Gmsh_Renderer.dielectric_gap_polygons.append((entry[0], entry[2].geometry))

        #create chip surface
        Palace_Gmsh_Renderer.metal_surface = shapely.geometry.box(self.center_x - 0.5*self.size_x, self.center_y - 0.5*self.size_y,
                             self.center_x + 0.5*self.size_x, self.center_y + 0.5*self.size_y)

        #Prep the design based on simulation type
        if simulation_type == 'cap_sim':
            self._cap_sim_prep()
            self._convert_design_to_gmsh_cap()
        if simulation_type == 'eigen_sim' or 'freq_sim':
            self._rf_sim_prep()
            self._convert_design_to_gmsh_rf()


    def _cap_sim_prep(self):   
        
        #get metal polygons
        metal_fuse_list = []
        for i, value in enumerate(Palace_Gmsh_Renderer.metal_polygons):
            metal_fuse_list.append(value[1])
        
        #Fuse the metal elements using Qutilities function (unary_union)
        fused_metals = QUtilities.fuse_polygons_threshold(metal_fuse_list)

        #plot fused metal polygons and add fused metal polygons to list
        fig, ax = plt.subplots()  # a figure with a single Axes
        
        #uncomment plot for debugging purposes
        for i, geom in enumerate(fused_metals.geoms):
            Palace_Gmsh_Renderer.fused_metals.append(geom)
            #poly = gpd.GeoSeries([geom])
            #poly.plot(ax=ax)

        #get dielectric gap polygons
        dielectric_fuse_list = []
        for i, value in enumerate(Palace_Gmsh_Renderer.dielectric_gap_polygons):
            dielectric_fuse_list.append(value[1])

        #Fuse the dielectric gap elements using Qutilities function (unary_union)
        fused_dielectric_gaps = QUtilities.fuse_polygons_threshold(dielectric_fuse_list)
       
        #plot fused dielectric polygons
        for i, geom in enumerate(fused_dielectric_gaps.geoms):
            Palace_Gmsh_Renderer.fused_dielectric_gaps.append(geom)
            #poly = gpd.GeoSeries([geom])
            #poly.boundary.plot(ax=ax, color='red', linewidth=0.3)
        
        poly_chip = gpd.GeoSeries([Palace_Gmsh_Renderer.metal_surface])
        poly_chip.boundary.plot(ax=ax, color='blue', linewidth=0.3)

        plt.show()

    
    def _rf_sim_prep(self):
        
        #get dielectric gap polygons
        dielectric_fuse_list = []
        for i, value in enumerate(Palace_Gmsh_Renderer.dielectric_gap_polygons):
            dielectric_fuse_list.append(value[1])

        #Fuse the dielectric gap elements using Qutilities function (unary_union)
        fused_dielectric_gaps = QUtilities.fuse_polygons_threshold(dielectric_fuse_list)

        #cutout dielectric gaps from metal ground plain
        ground_plain = shapely.difference(Palace_Gmsh_Renderer.metal_surface, fused_dielectric_gaps)

        #get metal polygons
        metal_fuse_list = []
        for i, value in enumerate(Palace_Gmsh_Renderer.metal_polygons):
            metal_fuse_list.append(value[1])
        
        #add ground plane to the metal fuse list
        metal_fuse_list.append(ground_plain)

        #Fuse the metal elements using Qutilities function (unary_union)
        fused_metals = QUtilities.fuse_polygons_threshold(metal_fuse_list)

        #plot fused metal polygons and add fused metal polygons to list
        
        fig2, ax2 = plt.subplots()  # a figure with a single Axes
        for i, geom in enumerate(fused_metals.geoms):
            Palace_Gmsh_Renderer.fused_metals.append(geom)
            poly = gpd.GeoSeries([geom])
            Palace_Gmsh_Renderer.gpd_polys.append(poly)
            poly.plot(ax=ax2)

        plt.show()
       

    def _convert_design_to_gmsh_rf(self):
        
        metal_list = []
        for i, metal in enumerate(Palace_Gmsh_Renderer.gpd_polys):
            metal_simplified = metal.simplify(1e-9)
            gmsh_exterior = self.draw_polygon_from_coords(metal_simplified.exterior[0].coords[:])
            if len(metal_simplified.interiors[0]) >= 1:
                interiors = []
                for _,interior in enumerate(metal_simplified.interiors[0]):
                    gmsh_interior = self.draw_polygon_from_coords(interior.coords[:])
                    interiors.append((2,gmsh_interior))
                gmsh_surface, gmsh_surface_map = gmsh.model.occ.cut([(2,gmsh_exterior)], 
                                                                    interiors,
                                                                    removeObject=True, removeTool=True)   
                Palace_Gmsh_Renderer.gmsh_metals.append(gmsh_surface[0][1])
                metal_list.append((2,gmsh_surface[0][1]))
                interiors.clear()
            else:
                Palace_Gmsh_Renderer.gmsh_metals.append(gmsh_exterior)
                metal_list.append((2,gmsh_exterior))

        #update geometries
        gmsh.model.geo.synchronize()
        gmsh.model.occ.synchronize()
        
        #create volume of dielectric base
        self.draw_chip_base()

        #update geometries
        gmsh.model.geo.synchronize()
        gmsh.model.occ.synchronize()
        
        #add launchpads
        for i,launch_pad in enumerate(Palace_Gmsh_Renderer.launch_pads_list):
            self.create_ports_on_launchpad(launch_pad)
   
        #lumped_ports to fragment with dielectric volume
        lumped_ports = []
        for i,value in enumerate(Palace_Gmsh_Renderer.ports_list):
            lumped_ports.append((2,value))

        elements_to_fragment = metal_list
        elements_to_fragment.extend(lumped_ports)

        #update geometries
        gmsh.model.geo.synchronize()
        gmsh.model.occ.synchronize()

        #fragment the newly created elements with the dielectric volume
        chip, chip_map = gmsh.model.occ.fragment([Palace_Gmsh_Renderer.dielectric_box[1]], 
                                                                 elements_to_fragment,
                                                                removeObject=True, removeTool=True)
        
        #update geometries
        gmsh.model.geo.synchronize()
        gmsh.model.occ.synchronize()

        #add physical group for ground plane
        gmsh.model.addPhysicalGroup(2, Palace_Gmsh_Renderer.gmsh_metals, name = 'metals')
        
        #add dielectric volume as dielectric base
        gmsh.model.addPhysicalGroup(3, [chip[0][1]], name = 'dielectric_base')
        Palace_Gmsh_Renderer.frag_dielectric_vol = chip[0][1]

        #add dielectric gaps
        start_index = len(Palace_Gmsh_Renderer.gmsh_metals) + 1 + 4 + 2*len(Palace_Gmsh_Renderer.launch_pads_list)
        end_index = len(chip) - 1
        gaps = chip[start_index:end_index]

        Palace_Gmsh_Renderer.gmsh_gap_list = []
        for _,gap in enumerate(gaps):
            Palace_Gmsh_Renderer.gmsh_gap_list.append(gap[1])
        
        gmsh.model.addPhysicalGroup(2, Palace_Gmsh_Renderer.gmsh_gap_list, name = 'dielectric_gaps')

        #draw airbox surrounfding chip
        self.draw_air_box()

        #get physical group indentifiers to use in config file
        self._identify_physical_groups()


    def _convert_design_to_gmsh_cap(self):
        
        #draw metal polygons into Gmsh
        metal_list = []
        for i, metal in enumerate(Palace_Gmsh_Renderer.fused_metals):
            gmsh_surface = self.draw_polygon(metal)
            Palace_Gmsh_Renderer.gmsh_metals.append(gmsh_surface)
            metal_list.append((2,gmsh_surface))

        #draw metal surface into Gmsh
        Palace_Gmsh_Renderer.gmsh_metal_surface = self.draw_polygon(Palace_Gmsh_Renderer.metal_surface)

        #draw dielectric gap polygons into Gmsh
        dielectric_gap_list = []
        for i, dielectric_gap in enumerate(Palace_Gmsh_Renderer.fused_dielectric_gaps):
            gmsh_surface = self.draw_polygon(dielectric_gap)
            Palace_Gmsh_Renderer.gmsh_dielectric_gaps.append(gmsh_surface)
            dielectric_gap_list.append((2,gmsh_surface))
        
        #update geometries
        gmsh.model.geo.synchronize()
        gmsh.model.occ.synchronize()

        #cut out dielectric gaps from metal surface
        ground_plane, ground_plane_map = gmsh.model.occ.cut([(2,Palace_Gmsh_Renderer.gmsh_metal_surface)], 
                                                                  dielectric_gap_list,
                                                                 removeObject=True, removeTool=True)

        #update geometries
        gmsh.model.geo.synchronize()
        gmsh.model.occ.synchronize()

        #add physical group for ground plane
        ground_plane_pieces = []
        ground_plane_pieces_to_fragment = []
        for i, value in enumerate(ground_plane):
            ground_plane_pieces.append(value[1])
            ground_plane_pieces_to_fragment.append((2,value[1]))
        
        #add metal elements (eg: resonator, feed line, etc) to ground plane
        metal, metal_map = gmsh.model.occ.fragment(ground_plane, 
                                                                 metal_list,
                                                                removeObject=True, removeTool=True)

        #create volume of dielectric base
        chip_base = self.draw_chip_base()
        
        #update geometries
        gmsh.model.geo.synchronize()
        gmsh.model.occ.synchronize()

        for i,value in enumerate(Palace_Gmsh_Renderer.gmsh_metals):
            gmsh.model.addPhysicalGroup(2, [value], name = 'metal_'+str(i))

        #metal elements to fragment with dielectric volume
        lumped_ports = []
        for i,value in enumerate(Palace_Gmsh_Renderer.ports_list):
            lumped_ports.append((2,value))

        elements_to_fragment = metal_list
        elements_to_fragment.extend(ground_plane)
        elements_to_fragment.extend(lumped_ports)

        #update geometries
        gmsh.model.geo.synchronize()
        gmsh.model.occ.synchronize()

        #fragment the newly created elements with the dielectric volume
        chip, chip_map = gmsh.model.occ.fragment([Palace_Gmsh_Renderer.dielectric_box[1]], 
                                                                 elements_to_fragment,
                                                                removeObject=True, removeTool=True)
        
        #update geometries
        gmsh.model.geo.synchronize()
        gmsh.model.occ.synchronize()

        #add physical group for ground plane
        gmsh.model.addPhysicalGroup(2, ground_plane_pieces, name = 'ground_plane')
        
        #add dielectric volume as dielectric base
        gmsh.model.addPhysicalGroup(3, [chip[0][1]], name = 'dielectric_base')
        Palace_Gmsh_Renderer.frag_dielectric_vol = chip[0][1]

        #draw airbox surrounfding chip
        self.draw_air_box()

        #get physical group indentifiers to use in config file
        self._identify_physical_groups()
         

    def _draw_components(self, component_list):
        '''Takes a list of Comoponents and '''
        
        for component_name in component_list:
            if(not self.design.qgeometry.get_component(component_name)['path'].empty):
                path_to_draw = self.design.qgeometry.get_component(component_name)['path']
                for index in path_to_draw.index:
                    self.draw_path(path_to_draw.loc[index], component_name, index, flag = 'path')
            if(not self.design.qgeometry.get_component(component_name)['poly'].empty):
                poly_to_draw = self.design.qgeometry.get_component(component_name)['poly']
                for index in poly_to_draw.index:
                    self.draw_polygon(poly_to_draw.loc[index], component_name, index, flag = 'poly')
            if(not self.design.qgeometry.get_component(component_name)['junction'].empty):
                junc_to_draw = self.design.qgeometry.get_component(component_name)['junction']
                for index in junc_to_draw.index:
                    self.draw_path(junc_to_draw.loc[index], component_name, index, flag = 'junction')



    def draw_chip_base(self):
        '''This method draws the chip base given the dimensions defined by the user'''

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
        Palace_Gmsh_Renderer.dielectric_box = gmsh.model.occ.extrude([(2, base_surface)],0,0,self.size_z)

        #update model with the created geometry items
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()

        return Palace_Gmsh_Renderer.dielectric_box

        

    def draw_polygon(self, polygon):
        '''takes a shapely polygon object or pandas series 
            in as an argument and then draws it in Gmsh and returns the surface ID''' 

        #simplify polygon, remove points that are too close together
        new_polygon = polygon.simplify(1e-9)

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


    def draw_polygon_from_coords(self, coords):
        '''draw polygon from x and y coordinates''' 

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
        '''takes a pandas series in as an argument and then draws it in Gmsh'''

        #get width and buffer amount for path
        width = path.width
        buffer_amt = width/2

        #fillet path using QMplRenderer
        qmpl = QMplRenderer(None, self.design, None)
        path_filleted = qmpl.fillet_path(path)

        #buffer the path by the user defined width
        poly_path = shapely.buffer(path_filleted, distance = buffer_amt, cap_style = 'flat')

        return poly_path


    def draw_air_box(self):
        
        #for air box choose increase in dimensions in 'mm'
        air_box_delta_x = (1/3) * self.size_x
        air_box_delta_y = (1/3) * self.size_y
        air_box_delta_z = 3 * self.size_z

        #half values of the sizes plus add increase for air box
        half_size_x = self.size_x/2 + air_box_delta_x/2
        half_size_y = self.size_y/2 + air_box_delta_y/2

        #store coordinates for surfaces of chip
        air_box_surface = {'point1': [self.center_x + half_size_x, self.center_y + half_size_y, (self.center_z + self.size_z) - air_box_delta_z],
                            'point2': [self.center_x + half_size_x, self.center_y - half_size_y, (self.center_z + self.size_z) - air_box_delta_z],
                            'point3': [self.center_x - half_size_x, self.center_y - half_size_y, (self.center_z + self.size_z) - air_box_delta_z],
                            'point4': [self.center_x - half_size_x, self.center_y + half_size_y, (self.center_z + self.size_z) - air_box_delta_z]}
        
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
        air_box = gmsh.model.occ.extrude([(2, surface)],0,0,(air_box_delta_z + 3*self.center_z))

        #cut out chip base from airbox
        air_box_cutout, air_box_cutout_map = gmsh.model.occ.fragment([(3, Palace_Gmsh_Renderer.frag_dielectric_vol)], [air_box[1]],  
                                                 removeObject=True, removeTool=True)

        #update model with the created geometry items
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()

        #add physical group for dielectric chip base
        gmsh.model.addPhysicalGroup(3, [air_box_cutout[1][1]], name = 'air_box')

        #add physical group for the surfaces of the air box which will represent the far field bondary conditions
        far_field = []
        for i, value in enumerate(gmsh.model.getEntities(2)[-7:]): #last 7 surface entities represnt the far field
            far_field.append(value[1])
        gmsh.model.addPhysicalGroup(2, far_field, name = 'far_field')



    def add_ports_on_launchpad(self, launch_pad):
        '''store launch pad objects to have ports created on'''
        Palace_Gmsh_Renderer.launch_pads_list.append(launch_pad)



    def create_ports_on_launchpad(self, launch_pad):
        '''create lumped port on launch pad for RF simulation'''

        #get indices for lp and pocket
        lp_index = self.design.qgeometry.get_component(launch_pad.name)['poly'].index[0]
        pocket_index = self.design.qgeometry.get_component(launch_pad.name)['poly'].index[1]

        #get pandas geodataframe for the launch pad
        lp_coords = self.design.qgeometry.get_component(launch_pad.name)['poly'].loc[lp_index].geometry.exterior.coords.xy
        pocket_coords = self.design.qgeometry.get_component(launch_pad.name)['poly'].loc[pocket_index].geometry.exterior.coords.xy

        #coordinates for first port
        start_point_lp_x_1 = lp_coords[0][2]
        start_point_lp_y_1 = lp_coords[1][2]
        start_point_poc_x = pocket_coords[0][2]
        start_point_poc_y = pocket_coords[1][2]
        x_diff = start_point_poc_x - start_point_lp_x_1
        y_diff = start_point_poc_y - start_point_lp_y_1

        #coordinates for second port
        start_point_lp_x_2 = lp_coords[0][3]
        start_point_lp_y_2 = lp_coords[1][3]

        #coordinates for second port
        scale = 0.25

        #Check orientation of the launchpad to define direction of lumped element ports
        if(str(launch_pad.options['orientation']) == '0' or str(launch_pad.options['orientation']) == '360'):
            dx1 = scale*-x_diff
            dy1 = y_diff
            dx2 = scale*-x_diff
            dy2 = -y_diff
        elif (str(launch_pad.options['orientation']) == '180' or str(launch_pad.options['orientation']) == '-180'):
            dx1 = scale*-x_diff
            dy1 = y_diff
            dx2 = scale*-x_diff
            dy2 = -y_diff
        elif (str(launch_pad.options['orientation'] == '270') or str(launch_pad.options['orientation']) == '-90'):
            dx1 = x_diff
            dy1 = scale*-y_diff
            dx2 = -x_diff
            dy2 = scale*-y_diff
        elif (str(launch_pad.options['orientation']) == '90' or str(launch_pad.options['orientation']) == '-270'):
            dx1 = x_diff
            dy1 = -scale*y_diff
            dx2 = -x_diff
            dy2 = -scale*-y_diff
        
        #draw the ports
        port1 = gmsh.model.occ.add_rectangle(start_point_lp_x_1, start_point_lp_y_1, self.center_z, dx1, dy1)
        port2 = gmsh.model.occ.add_rectangle(start_point_lp_x_2, start_point_lp_y_2, self.center_z, dx2, dy2)

        Palace_Gmsh_Renderer.ports_list.append(port1)
        Palace_Gmsh_Renderer.ports_list.append(port2)

        gmsh.model.addPhysicalGroup(2, [port1], name = 'port_'+launch_pad.name+'a')
        gmsh.model.addPhysicalGroup(2, [port2], name = 'port_'+launch_pad.name+'b')

        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()



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



    def _identify_physical_groups(self):

        physical_groups = gmsh.model.get_physical_groups()

        for i, value in enumerate(physical_groups):
            phys_group_name = gmsh.model.get_physical_name(value[0], value[1])
            if phys_group_name == 'metals':
                    Palace_Gmsh_Renderer.config_metals_rf.append(value[1])
            elif phys_group_name == 'ground_plane':
                    Palace_Gmsh_Renderer.config_ground_plane = value[1]
                    Palace_Gmsh_Renderer.config_metals_rf.append(value[1])
                    Palace_Gmsh_Renderer.config_metals_cap.append(value[1])
            elif phys_group_name == 'far_field':
                    Palace_Gmsh_Renderer.config_far_field = value[1]
            elif phys_group_name == 'dielectric_base':
                    Palace_Gmsh_Renderer.config_dielectric_base = value[1]
            elif phys_group_name == 'dielectric_gaps':
                    Palace_Gmsh_Renderer.config_dielectric_gaps = value[1]
            elif phys_group_name == 'air_box':
                    Palace_Gmsh_Renderer.config_air_box = value[1]
            elif phys_group_name[0:4] == 'port':
                    Palace_Gmsh_Renderer.config_ports.append(value[1])
            elif phys_group_name[0:6] == 'metal_':
                    Palace_Gmsh_Renderer.config_metals_cap.append(value[1])
        



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

    
    def _intelligent_mesh(self, simulation_type, min_size, max_size, mesh_sampling):

        #turn off meshing parameters that are not required
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 5e-3)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 5e-3)

        #set mesh algorithm Mesh.Algorithm3D to HXT
        gmsh.option.setNumber("Mesh.Algorithm3D", 10) 

        #define distance field
        gmsh.model.mesh.field.add("Distance", 1)
        gmsh.model.mesh.field.setNumber(1, "Sampling", mesh_sampling)

        if simulation_type == 'capacitance_simulation':
            gmsh.model.mesh.field.setNumbers(1, "SurfacesList", Palace_Gmsh_Renderer.gmsh_metals)
        if simulation_type == 'eigenmode_simulation' or 'frequency_driven_simulation':
            gmsh.model.mesh.field.setNumbers(1, "SurfacesList", Palace_Gmsh_Renderer.gmsh_gap_list)

        #define mesh field
        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", min_size)
        gmsh.model.mesh.field.setNumber(2, "SizeMax", max_size)
        gmsh.model.mesh.field.setNumber(2, "DistMin", 5*min_size)
        gmsh.model.mesh.field.setNumber(2, "DistMax", 7*min_size)

        #set background field and render mesh
        gmsh.model.mesh.field.setAsBackgroundMesh(2)
        gmsh.option.setNumber('Mesh.MshFileVersion', 2.2)
        gmsh.model.mesh.generate(dim = 3)
        




        