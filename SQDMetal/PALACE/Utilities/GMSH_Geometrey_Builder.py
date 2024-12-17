from qiskit_metal.qgeometries.qgeometries_handler import QGeometryTables
from qiskit_metal.renderers.renderer_mpl.mpl_renderer import QMplRenderer
from Utilities.QiskitShapelyRenderer import QiskitShapelyRenderer
from Utilities.ShapelyEx import ShapelyEx
from qiskit_metal.qlibrary.terminations.launchpad_wb import LaunchpadWirebond
from Utilities.QUtilities import QUtilities
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

    def __init__(self, design, simulation_type, ports):
        
        self.design = design 
        self.simulation_type = simulation_type
        self.ports = ports

        #Get dimensions of chip base and convert to design units in 'mm'
        self.center_x = self.design.parse_value(self.design.chips['main'].size.center_x)
        self.center_y = self.design.parse_value(self.design.chips['main'].size.center_y)
        self.center_z = self.design.parse_value(self.design.chips['main'].size.center_z)
        self.length_x = self.design.parse_value(self.design.chips['main'].size.size_x)
        self.length_y = self.design.parse_value(self.design.chips['main'].size.size_y)
        self.length_z = self.design.parse_value(self.design.chips['main'].size.size_z)

        #Initialize the GMSH API and name the model
        gmsh.initialize()

      
    def construct_geometry_in_GMSH(self):
        '''This function takes the existing geometry in the Qiskit Metal design and constructs them in GMSH.
        
        Args:
            None.

        Returns:
            
        '''

        #Do pre-processing in shapely to get metallic elements and dielectric cutouts ready to build in GMSH.
        #Note: metals list contains ground plane. Dielectric gaps are the difference between the dielectric cutout
        #and the metals
        metals, ground_plane, dielectric_gaps, dielectric_cutouts = self._process_qiskit_geometries_in_shapely()

        #Draw shapely metal and dielectric gap polygons into GMSH. If the polygon has an interior, the interior sections are subtracted from
        #the exterior boundary
        metal_list = self._create_gmsh_geometry_from_shapely_polygons(metals)                           #create all metal traces
        metal_meshing_list = self._create_gmsh_geometry_from_shapely_polygons(metals)                   #this list is used for meshing
        dielectric_cutout_list = self._create_gmsh_geometry_from_shapely_polygons(dielectric_cutouts)   #dielectric cutouts to be used for meshing 
        dielectric_gap_list = self._create_gmsh_geometry_from_shapely_polygons(dielectric_gaps)         #create all the gaps between the metal traces and the ground plane 
        ground_plane_list = self._create_gmsh_geometry_from_shapely_polygons(ground_plane)              #creates all the pieces of the ground plane

        #Plot the shapely metals for user to see device
        geoms = metals + ground_plane
        metal_names = [str(i) for i,_ in enumerate(geoms)]
        gdf = gpd.GeoDataFrame({'names':metal_names}, geometry=geoms)
        fig, ax = plt.subplots()
        gdf.plot(ax = ax, column='names', cmap='tab10', alpha=0.75, categorical=True, legend=True)
        plt.show()

        #Create the chip base (dielectric substrate) in Gmsh
        chip_base = self._create_chip_base()
        
        #Create the airbox which houses the chip base
        bottom_grounded = True
        air_box = self._draw_air_box(bottom_grounded)

        print('Creating geometry for', self.simulation_type, 'Simulation.')

        #process ports
        ports_list, ports_dict = self._process_ports()

        #list for metals for capacitance simulation
        metal_cap_physical_group = []
        metal_cap_names = []

        #The geometry of the design needs to be altered depending on what time of simulation is being run
        if self.simulation_type == 'Eigenmode' or self.simulation_type == 'Driven':
            
            #Fuse the ground plane and any metallic element which are shorted to ground plane i.e. lambda/4 resonators
            metal_pieces = [(x[1],x[2]) for x in metal_list]
            ground_plane_pieces = [(x[1],x[2]) for x in ground_plane_list]
            metals_and_gp, metals_and_gp_map = gmsh.model.occ.fuse(ground_plane_pieces, metal_pieces, removeObject=True, removeTool=True)
            
            #create JJ's as lumped ports, if they exist
            junctions = []
            if not self.design.qgeometry.tables['junction'].empty:
                junctions, jj_inductance = self._create_JJ_lumped_ports()
            
            #create junction list for boolean operations in gmsh
            junction_list = [(2,x) for x in junctions]

            #if there are junctions in the design, create the physical group
            jj_dict = {}
            if junctions:
                for i,junction in enumerate(junctions):
                    jj_name = 'jj_' + str(i)
                    jj_physical_group = gmsh.model.addPhysicalGroup(2, [junction], name = jj_name)
                    jj_dict[jj_name] = (jj_physical_group, jj_inductance)

            #dielectric gap list needs to be updated to account for the presence of ports and junctions
            cut_list = ports_list + junction_list
            dielectric_gap_list_for_cut = [(x[1],x[2]) for x in dielectric_gap_list]
            dielectric_gap_list,map = gmsh.model.occ.cut(dielectric_gap_list_for_cut, cut_list, removeObject=True, removeTool=False)
            
            gmsh.model.occ.synchronize()
            gmsh.model.geo.synchronize()

            #Fragment the newly created fused elements with the dielectric volume
            fragment_list = metals_and_gp + dielectric_gap_list + ports_list + junction_list
            chip, chip_map = gmsh.model.occ.fragment([(3,chip_base)], fragment_list, removeObject=True, removeTool=True)

            #Create a single physical group for all the metals in the design
            metals_group = [x[1] for x in metals_and_gp]
            gmsh.model.addPhysicalGroup(2, metals_group, name = 'metals')

            #create physical group for dielectric gaps
            gmsh.model.addPhysicalGroup(2, [x[1] for x in dielectric_gap_list], name = 'dielectric_gaps')

           
        elif self.simulation_type == 'Capacitance':

            #Fragment the newly created elements with the dielectric volume
            fragment_list = metal_list + dielectric_gap_list + ground_plane_list
            fragment_list = [(x[1],x[2]) for x in fragment_list]
            chip, chip_map = gmsh.model.occ.fragment([(3,chip_base)], fragment_list, removeObject=True, removeTool=True)

            #Create a physical group for each metal in the design. This is important for capacitance simulations because we need to uniquely
            #design the
            metal_list = metal_list + ground_plane_list 
            
            for i,metal in enumerate(metal_list):
                metal_name = 'metal_' + str(i)
                metal_physical_group = gmsh.model.addPhysicalGroup(2, [metal[2]], name = metal_name)
                metal_cap_physical_group.append(metal_physical_group)
                metal_cap_names.append(metal_name)

            #add physicla group for dielectric gap list
            gmsh.model.addPhysicalGroup(2, [x[2] for x in dielectric_gap_list], name = 'dielectric_gaps')

            #no JJs used in capacitance simulation
            jj_dict = {}
            
        gmsh.model.geo.synchronize()
        gmsh.model.occ.synchronize()

        #Add in physical groups for each part of the device
        gmsh.model.addPhysicalGroup(3, [chip_base], name = 'dielectric_substrate')
        #gmsh.model.addPhysicalGroup(2, [x[2] for x in ground_plane_list], name = 'ground_plane')

        #Fragment the dielectric volume with the airbox
        chip_and_air_box, chip_and_air_box_map = gmsh.model.occ.fragment([(3,chip_base)], [(3, air_box)], removeObject=True, removeTool=True)
        
        gmsh.model.occ.synchronize()
        gmsh.model.geo.synchronize()
        
        gmsh.model.addPhysicalGroup(3, [air_box], name = 'air_box')
        
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

        gmsh.model.addPhysicalGroup(2, far_field_surfaces, name = 'far_field')
        
        #synchronise
        gmsh.model.geo.synchronize()
        gmsh.model.occ.synchronize()

        print('Geometry successfully built in Gmsh.')

        return metal_list, metal_meshing_list, dielectric_gap_list, dielectric_cutout_list, ports_dict, metal_cap_physical_group, metal_cap_names, jj_dict


    def _process_qiskit_geometries_in_shapely(self):
        '''This function takes the existing geometry in the Qiskit Metal design and processes them using shapely
            to get them ready to build in GMSH. Processing includes fusing metallic elements together such as fusing the launch
            pads and the CPW transmission line together.
        
        Args:
            None.

        Returns:
            List of polygons for the metallic elements and dielectric cutouts in the design.
        '''

        unit_conv = 1 #convert all units into metres
        fillet_resolution = 15 #improve resolution around curves i.e. bends in CPW resonators

        #Use renderer to get all shapely objects in Qiskit Metal design
        QSR = QiskitShapelyRenderer(None, self.design, None)

        #Get the coordinates of all the objects in the Qiskit Metal design
        design_objects = QSR.get_net_coordinates(fillet_resolution)
        
        #Filter design objects to get the shapely objects which correspond to cutouts in the ground plane and metals
        filtered_cutouts = design_objects.loc[design_objects['subtract'] == True]
        filtered_metals = design_objects.loc[design_objects['subtract'] == False]

        #Get dielectric cutouts and metals - buffer 0 trick used
        dielectric_cutouts = ShapelyEx.fuse_polygons_threshold(filtered_cutouts['geometry'].buffer(0), 1e-12)
        metals = ShapelyEx.fuse_polygons_threshold(filtered_metals['geometry'].buffer(0), 1e-12)

        #Make sure all dielectric cutouts are polygons
        if isinstance(dielectric_cutouts, shapely.geometry.multipolygon.MultiPolygon):
            dielectric_cutouts = [x for x in dielectric_cutouts.geoms]
        else:
            dielectric_cutouts = [dielectric_cutouts] #i.e. it's just a lonely Polygon object...
        dielectric_cutouts = [shapely.affinity.scale(x, xfact=unit_conv, yfact=unit_conv, origin=(0,0)) for x in dielectric_cutouts]

        #Make sure all metals are polygons
        if isinstance(metals, shapely.geometry.multipolygon.MultiPolygon):
            metals = [x for x in metals.geoms]
        else:
            metals = [metals] #i.e. it's just a lonely Polygon object...
        metals = [shapely.affinity.scale(x, xfact=unit_conv, yfact=unit_conv, origin=(0,0)) for x in metals]

        #Create metal surface as basis for ground plane
        metal_surface = shapely.geometry.box(self.center_x - 0.5*self.length_x, self.center_y - 0.5*self.length_y,
                             self.center_x + 0.5*self.length_x, self.center_y + 0.5*self.length_y)
        
        #Cut dielectric cutouts from metal surface to create the ground plane for the chip and then add to metals
        ground_plane = shapely.difference(metal_surface, shapely.geometry.multipolygon.MultiPolygon(dielectric_cutouts))
        
        #Make ground plane elements are polygons
        if isinstance(ground_plane, shapely.geometry.multipolygon.MultiPolygon):
            ground_plane = [x for x in ground_plane.geoms]
        else:
            ground_plane = [ground_plane] #i.e. it's just a lonely Polygon object...
        ground_plane = [shapely.affinity.scale(x, xfact=unit_conv, yfact=unit_conv, origin=(0,0)) for x in ground_plane]

        #Cut the metals from the dielectric cutouts to produce dielectric gaps
        dielectric_gaps = shapely.difference(dielectric_cutouts, shapely.geometry.multipolygon.MultiPolygon(metals))
        
        return metals, ground_plane, dielectric_gaps, dielectric_cutouts


    def _create_gmsh_geometry_from_shapely_polygons(self, polygons):

        polygons_list = []
        for m,poly in enumerate(polygons):
            p = gpd.GeoSeries(poly) #convert to geoseries
            poly_simplified = poly.simplify(1e-6) #this removes points that are spaced too closely together
            gmsh_exterior = self._draw_polygon_in_GMSH_from_coords(poly_simplified.exterior.coords[:-1]) #remove last coord from list because it's repeated
            if len(poly_simplified.interiors) >= 1:
                interiors = []
                for _,interior in enumerate(poly_simplified.interiors):
                    gmsh_interior = self._draw_polygon_in_GMSH_from_coords(interior.coords[:-1]) #remove last coord from list because it's repeated
                    interiors.append((2,gmsh_interior))
                gmsh_surface, gmsh_surface_map = gmsh.model.occ.cut([(2,gmsh_exterior)], 
                                                                    interiors,
                                                                    removeObject=True, removeTool=True)
                polygons_list.append((m, 2,gmsh_surface[0][1]))
            else:
                polygons_list.append((m, 2,gmsh_exterior))
        
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
        chip_base = gmsh.model.occ.addBox(self.center_x-self.length_x/2, self.center_y-self.length_y/2, self.center_z, self.length_x, self.length_y, -np.abs(self.length_z))
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
        air_box_delta_x = increase * self.length_x
        air_box_delta_y = increase * self.length_y
        
        #bottom left corner of air box
        x_point = self.center_x-self.length_x/2-air_box_delta_x/2 
        y_point = self.center_y-self.length_y/2-air_box_delta_y/2

        #Check to place ground plane on back side by choosing statring z point
        if bottom_grounded == True: 
            z_point = self.center_z - np.abs(self.length_z)
            air_box_delta_z = 2 * np.abs(self.length_z)
        else:
            z_point = self.center_z - 2 * np.abs(self.length_z)
            air_box_delta_z = 3 * np.abs(self.length_z)

        air_box = gmsh.model.occ.addBox(x_point, y_point, z_point, self.length_x + air_box_delta_x, self.length_y + air_box_delta_y, air_box_delta_z)
        gmsh.model.occ.synchronize()

        return air_box
    
    def _process_ports(self):
        
        ports_list = [] #list to store Gmsh identifier after port is drawn into Gmsh
        ports_dict ={} #dictionary to store port subcomopentes, including name, physical group (boundary condition) and orientation

        #If there are ports start processing
        if len(self.ports) != 0:
            
            #get component object from qiskit metal
            qObj = self.design.components[self.ports[0]]

            if isinstance(qObj, LaunchpadWirebond):
                for i,port in enumerate(self.ports):
                    
                    #for each launch pad draw in the ports list there are two terminations to ground
                    launchesA, launchesB, vec_perp = QUtilities.get_RFport_CPW_coords_Launcher(self.design, port, 20e-3, 1e3)
                    
                    #check port orientation
                    port_orientation = self._check_port_orientation(vec_perp)

                    #draw in first component
                    portA = self._draw_polygon_in_GMSH_from_coords(launchesA)
                    port_name_a = port + 'a'
                    port_phys_group_a = gmsh.model.addPhysicalGroup(2, [portA], name = port_name_a)
                    ports_list.append((2,portA))

                    #draw in second component
                    portB = self._draw_polygon_in_GMSH_from_coords(launchesB)
                    port_name_b = port + 'b'
                    port_phys_group_b = gmsh.model.addPhysicalGroup(2, [portB], name = port_name_b)
                    ports_list.append((2,portB))

                    ports_dict['port_' + str(i+1)] = {port_name_a: (port_phys_group_a, port_orientation[0]),
                                                port_name_b: (port_phys_group_b, port_orientation[1])}
                    
        else: 
            print('No ports for processing.')

        return ports_list, ports_dict
    

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
    
