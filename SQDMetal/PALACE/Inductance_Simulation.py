# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0

from SQDMetal.PALACE.Model import PALACE_Model_Base
from SQDMetal.COMSOL.Model import COMSOL_Model
from SQDMetal.COMSOL.SimCapacitance import COMSOL_Simulation_CapMats
from SQDMetal.Utilities.Materials import Material
from SQDMetal.Utilities.QubitDesigner import FloatingTransmonDesigner
from SQDMetal.Utilities.QUtilities import QUtilities
from SQDMetal.Utilities.ShapelyEx import ShapelyEx
import matplotlib.pyplot as plt
import json
import os
import gmsh
import shapely
import numpy as np
import pandas as pd
import geopandas as gpd
from SQDMetal.PALACE.Utilities.GMSH_Geometry_Builder import GMSH_Geometry_Builder
from SQDMetal.PALACE.Utilities.GMSH_Mesh_Builder import GMSH_Mesh_Builder
from SQDMetal.PALACE.PVDVTU_Viewer import PVDVTU_Viewer

class PALACE_Inductance_Simulation(PALACE_Model_Base):

    #Class Variables
    default_user_options = {
                 "solns_to_save": -1,
                 "solver_order": 2,
                 "solver_tol": 1.0e-8,
                 "solver_maxits": 100,
                 "HPC_Parameters_JSON": ""
    }

    # Parent Directory path
    HPC_parent_simulation_dir = "C:/PALACE_Simulations/"
    simPC_parent_simulation_dir = "/home/experiment/PALACE/Simulations/input"

    #constructor
    def __init__(self, name, sim_parent_directory, mode, meshing, user_options = {}, 
                                    view_design_gmsh_gui = False, create_files = False, **kwargs):
        self.name = name
        self.sim_parent_directory = sim_parent_directory
        self.mode = mode
        self.meshing = meshing
        self.user_options = {}
        for key in PALACE_Inductance_Simulation.default_user_options:
            self.user_options[key] = user_options.get(key, PALACE_Inductance_Simulation.default_user_options[key])
        self.view_design_gmsh_gui = view_design_gmsh_gui
        self.create_files = create_files
        self._cur_cap_terminals = []
        self.cap_matrix = None
        self._ports = []
        self._int_areas = []
        super().__init__(meshing, mode, user_options, **kwargs)

    def _create_directory(self, directory_name):
        '''create a directory to hold the simulation files'''

        if self.create_files:
            parent_simulation_dir = self._check_simulation_mode()

            # Directory
            directory = directory_name
    
            # Path
            path = os.path.join(parent_simulation_dir, directory)
    
            # Create the directory
            if not os.path.exists(path):
                os.mkdir(path)
                # print("Directory '% s' created" % directory)

    def _save_mesh_gmsh(self):
        parent_simulation_dir = self._check_simulation_mode()

        # file_name
        file_name = self.name + "/" + self.name + ".msh"
  
        # Path
        path = os.path.join(parent_simulation_dir, file_name)
        gmsh.write(path)
        
    def _save_mesh_comsol(self, comsol_obj):
        '''function used to save the comsol mesh file'''
        if self.create_files:
            parent_simulation_dir = self._check_simulation_mode()

            # file_name
            file_name = self.name + "/" + self.name + ".mphbin"
    
            # Path
            path = os.path.join(parent_simulation_dir, file_name)

            #COMSOL export commands
            comsol_obj._model.java.component("comp1").mesh("mesh1").export().set("filename", path)
            comsol_obj._model.java.component("comp1").mesh("mesh1").export(path)

    def _create_config_file(self, **kwargs):
        '''create the configuration file which specifies the simulation type and the parameters'''    

        if self.meshing == 'GMSH':

            #GMSH renderer object needed to get boundary conditions for config file
            gmsh_render_attrs = kwargs['gmsh_render_attrs']

            #GMSH config file variables
            material_air = gmsh_render_attrs['air_box']
            material_dielectric = gmsh_render_attrs['dielectric']
            far_field = [gmsh_render_attrs['far_field'][x] for x in gmsh_render_attrs['far_field']]
            ports = gmsh_render_attrs['ports']
            PEC_metals = gmsh_render_attrs['metals']

            int_areas = gmsh_render_attrs['integration_areas']

            #define length scale
            l0 = 1e-3

            #file extension
            file_ext = '.msh'
        
        if self.meshing == 'COMSOL':
            assert False, "Not supported."

        #get material parameters
        dielectric = Material(self.user_options["dielectric_material"])

        #Process current sources
        current_sources = []
        for m, cur_port in enumerate(self._ports):
            current_sources.append({
                'Index': m+1,
                'Attributes': [ports[cur_port['port_name']]],
                'Direction': cur_port['vec_field'] + [0]
            })
        
        #Process integration areas
        post_procs = []
        for m,cur_area in enumerate(int_areas):
            post_procs.append({
                'Index': m+1,
                'Attributes': cur_area,
                'Type': "Magnetic",
	            'Center': [0, 0, -1]    #Positive in +Z Not relevant to calculation, but just to make sure normal is pointing up (i.e. it can also be down for a 2D surface)
            })

        #Define python dictionary to convert to json file
        if self._output_subdir == "":
            self.set_local_output_subdir("")
        filePrefix = self.hpc_options["input_dir"]  + self.name + "/" if self.hpc_options["input_dir"] != "" else ""
        self._mesh_name = filePrefix + self.name + file_ext
        config = {
                    "Problem":
                    {
                        "Type": "Magnetostatic",
                        "Verbose": 2,
                        "Output": filePrefix  + "outputFiles"
                    },
                    "Model":
                    {
                        "Mesh":  self._mesh_name,
                        "L0": l0,  # mm
                        "Refinement": self._mesh_refinement
                    },
                    "Domains":
                    {
                        "Materials":
                        [
                            {
                                "Attributes": material_air,  # Air
                                "Permeability": 1.0,
                                "Permittivity": 1.0,
                                "LossTan": 0.0
                            },
                            {
                                "Attributes": material_dielectric,  # Dielectric
                                "Permeability": dielectric.permeability,
                                "Permittivity": dielectric.permittivity,
                                "LossTan": dielectric.loss_tangent
                            }
                        ]
                    },
                    "Boundaries":
                    {
                        "PEC": {
                            'Attributes': far_field + PEC_metals
                        },
                        "SurfaceCurrent": current_sources,
                        "Postprocessing":
                        {
                        "SurfaceFlux": post_procs,
                        }
                    },
                    "Solver":
                    {
                        "Order": self.user_options["solver_order"],
                        "Magnetostatic":
                        {
                        "Save": len(self._ports) if self.user_options["solns_to_save"] == -1 else self.user_options["solns_to_save"]
                        },
                        "Linear":
                        {
                        "Type": "AMS",
                        "KSPType": "CG",
                        "Tol": self.user_options["solver_tol"],
                        "MaxIts": self.user_options["solver_maxits"]
                        }
                    }
                }

        #check simulation mode and return appropriate parent directory 
        parent_simulation_dir = self._check_simulation_mode()

        #destination for config file
        simulation_dir = parent_simulation_dir + str(self.name)

        #simulation file name
        sim_file_name = self.name+'.json'

        #save to created directory
        file = os.path.join(simulation_dir, sim_file_name)

        #write to file
        with open(file, "w+") as f:
            json.dump(config, f, indent=2)
        self.path_mesh = os.path.join(simulation_dir, self._mesh_name)
        self._sim_config = file
        self.set_local_output_subdir(self._output_subdir)

    def prepare_simulation(self):
        '''Creates and saves the GMSH (.msh) and Palace configuration file (.json). 
        '''

        if self.meshing == 'GMSH':
            lePorts = []
            for cur_port in self._ports:
                if cur_port['elem_type'] == 'single':
                    lePorts.append({'type':'lumped', 'name':cur_port['port_name'], 'coords':cur_port['portCoords']})

            ggb = GMSH_Geometry_Builder(self._geom_processor, self.user_options['fillet_resolution'], self.user_options['gmsh_verbosity'])
            gmsh_render_attrs = ggb.construct_geometry_in_GMSH(self._metallic_layers, self._ground_plane, lePorts,
                                                               self._fine_meshes, self.user_options["fuse_threshold"],
                                                               threshold=self.user_options["threshold"],
                                                               simplify_edge_min_angle_deg=self.user_options["simplify_edge_min_angle_deg"],
                                                               full_3D_params = self._full_3D_params,
                                                               boundary_distances = self._boundary_distances,
                                                               integration_areas = self._int_areas)
            #
            gmb = GMSH_Mesh_Builder(gmsh_render_attrs['fine_mesh_elems'], self.user_options)
            gmb.build_mesh()

            if self.create_files:
                #create directory to store simulation files
                self._create_directory(self.name)

                #create config file
                self._create_config_file(gmsh_render_attrs = gmsh_render_attrs)

                #create batch file
                # if self.mode == 'HPC':
                #     self.create_batch_file()

                self._save_mesh_gmsh()
            
        if self.meshing == 'COMSOL':
            assert False, "COMSOL Meshing not supported here."

    def create_current_source_with_Uclip_on_Launcher(self, qObjName:str, thickness_side:float=20e-6, thickness_back:float=20e-6, separation_gap:float=20e-6):
        """
        *Applies Qiskit-Metal designs only*

        Create a U-Clip that adjoins the two ground planes on a CPW launcher feeding from the edge of the chip. That
        is, a metallic piece (that is coplanar with the CPW) goes off the chip and comes around to join the two ground
        planes. It is a simple concave octogan with 90° corners. It can also have an optional RF port.

        This function must be run before calling :func:`~SQDMetal.PALACE.Model.PALACE_Model_Base.prepare_simulation`.

        .. figure:: /_static/palace_sim_Uclip_dims.drawio.svg
            :alt: Parameters used in the U-clip
            :align: center
            :scale: 100%

            Parameters used in the U-clip. The orange portion is the optional single-element RF port (arrow indicating excitation direction).

        Parameters
        ----------
        qObjName : float
            Name of the LauncherWB Q-Component in the Qiskit-Metal design.
        thickness_side : float
            Distance the U-clip, in metres, goes into the ground plane from the edge. Defaults to 20e-6.
        thickness_back : float
            Thickness of the furthest section from the chip (but also parallel with the chip) in metres. Defaults to 20e-6
        separation_gap : float
            Distance of the section parallel with the chip from the edge of the chip in metres. The default value is set
            to zero whereupon, it will use the gap distance of the CPW.
        """
        self._metallic_layers += [{
            'type': 'Uclip',
            'clip_type':'inplaneLauncher',
            'qObjName':qObjName,
            'thickness_side':thickness_side,
            'thickness_back':thickness_back,
            'separation_gap':separation_gap
        }]
        #
        vec_ori, vec_launch, cpw_wid, cpw_gap = QUtilities._get_LauncherWB_params(self._geom_processor.design, qObjName)
        unit_conv = QUtilities.get_units(self._geom_processor.design)
        pos1 = vec_ori - vec_launch*separation_gap
        pos2 = vec_ori
        rect_width = cpw_wid

        v_parl = pos2-pos1
        v_parl /= np.linalg.norm(v_parl)

        portCoords = ShapelyEx.rectangle_from_line(pos1, pos2, rect_width, False)
        portCoords = [x for x in portCoords]
        portCoords = portCoords + [portCoords[0]]   #Close loop...

        port_name = "I_port_" + str(len(self._ports))
        self._ports += [{'port_name':port_name, 'type':'single_rect', 'elem_type':'single', 'pos1':pos1, 'pin2':pos2, 'rect_width': rect_width,
                'vec_field': v_parl.tolist(),
                'portCoords': portCoords}]

    def create_current_source_with_Uclip_on_Route(self, route_name:str, pin_name:str, thickness_side:float=20e-6, thickness_back:float=20e-6, separation_gap:float=20e-6):
        """
        *Applies Qiskit-Metal designs only*

        Create a U-Clip that adjoins the two ground planes on a CPW route feeding from the edge of the chip. That is,
        a metallic piece (that is coplanar with the CPW) goes off the chip and comes around to join the two ground
        planes. It is a simple concave octogan with 90° corners. It can also have an optional RF port.

        This function must be run before calling :func:`~SQDMetal.PALACE.Model.PALACE_Model_Base.prepare_simulation`.

        .. figure:: /_static/palace_sim_Uclip_dims.drawio.svg
            :alt: Parameters used in the U-clip
            :align: center
            :scale: 100%

            Parameters used in the U-clip. The orange portion is the optional single-element RF port (arrow indicating excitation direction).

        Parameters
        ----------
        route_name : str
            Name of the routing component (e.g. a CPW) in the Qiskit-Metal design
        pin_name : str
            Name of the pin within the routing component in the Qiskit-Metal design. Typically a routing component has,
            at least, two pins named something like 'start' and 'end'.
        thickness_side : float
            Distance the U-clip, in metres, goes into the ground plane from the edge. Defaults to 20e-6.
        thickness_back : float
            Thickness of the furthest section from the chip (but also parallel with the chip) in metres. Defaults to 20e-6
        separation_gap : float
            Distance of the section parallel with the chip from the edge of the chip in metres. The default value is set
            to zero whereupon, it will use the gap distance of the CPW.
        """
        self._metallic_layers += [{
            'type': 'Uclip',
            'clip_type':'inplaneRoute',
            'route_name':route_name,
            'pin_name':pin_name,
            'thickness_side':thickness_side,
            'thickness_back':thickness_back,
            'separation_gap':separation_gap
        }]
        #
        vec_ori, vec_launch, cpw_wid, cpw_gap = QUtilities._get_Route_params(self._geom_processor.design, route_name, pin_name)
        unit_conv = QUtilities.get_units(self._geom_processor.design)
        pos1 = vec_ori - vec_launch*separation_gap
        pos2 = vec_ori
        rect_width = cpw_wid

        v_parl = pos2-pos1
        v_parl /= np.linalg.norm(v_parl)

        portCoords = ShapelyEx.rectangle_from_line(pos1, pos2, rect_width, False)
        portCoords = [x for x in portCoords]
        portCoords = portCoords + [portCoords[0]]   #Close loop...

        port_name = "I_port_" + str(len(self._ports))
        self._ports += [{'port_name':port_name, 'type':'single_rect', 'elem_type':'single', 'pos1':pos1, 'pin2':pos2, 'rect_width': rect_width,
                'vec_field': v_parl.tolist(),
                'portCoords': portCoords}]
        #TODO: Do this in general for GDS as well?


    def add_integration_area(self, integration_area:shapely.Polygon|shapely.geometry.multipolygon.MultiPolygon):
        self._int_areas.append(integration_area)

    def retrieve_data(self):
        '''
        Retrieves flux per 
        
        Creates and saves B-field plots due to the different terminal injections (input1_normB.png and input1_normBz.png),
        and mesh (mesh.png)

        This function must be run after calling :func:`~SQDMetal.PALACE.Model.PALACE_Model_Base.run`.

        Returns:
            The flux per unit current (that is, Wb/A) through the different integration surfaces due to the different
            terminal current injections. If no surfaces are specified, then it just returns the dictionary specified
            in the function :func:`~SQDMetal.PALACE.PALACE_Inductance_Simulation.retrieve_magnetostatic_data_from_file`.
        '''
        raw_data = self.retrieve_magnetostatic_data()
        if raw_data['surface_Flux'].size > 0:
            raw_data = raw_data['surface_Flux'] / raw_data['terminal_I']

        lePlots = self._output_data_dir + '/paraview/magnetostatic/magnetostatic.pvd'
        if os.path.exists(lePlots):
            leView = PVDVTU_Viewer(lePlots)
            for m in range(leView.num_datasets):
                try:
                    leSlice = leView.get_data_slice(m, slice_plane_origin=np.array([0,0,1e-6]))
                    if 'B' in leSlice.get_params():
                        fig = leSlice.plot(np.linalg.norm(leSlice.get_data('B'), axis=1), 'coolwarm', True)
                        fig.savefig(self._output_data_dir + f'/input{m}_normB.png')
                        plt.close(fig)
                        fig = leSlice.plot(leSlice.get_data('B')[:,2], 'coolwarm', True)
                        fig.savefig(self._output_data_dir + f'/input{m}_normBz.png')
                        plt.close(fig)
                except Exception as e:
                    print(f"Error in plotting: {e}")
            try:
                fig = leSlice.plot_mesh()
                fig.savefig(self._output_data_dir + '/mesh.png')
                plt.close(fig)
            except Exception as e:
                print(f"Error in plotting: {e}")

        return raw_data

    @staticmethod
    def retrieve_magnetostatic_data_from_file(output_directory):
        raw_data = pd.read_csv(output_directory + '/surface-F.csv')
        surf_Flux = raw_data.to_numpy()[:,1:]    #First column is just the indices...
        raw_data = pd.read_csv(output_directory + '/terminal-I.csv')
        term_I = raw_data.to_numpy()[:,1:]    #First column is just the indices...
        raw_data = pd.read_csv(output_directory + '/terminal-M.csv')
        term_M = raw_data.to_numpy()[:,1:]    #First column is just the indices...

        return {
            'surface_Flux': surf_Flux,
            'terminal_I': term_I,
            'inductance_matrix': term_M
        }


    def retrieve_magnetostatic_data(self):
        return PALACE_Inductance_Simulation.retrieve_magnetostatic_data_from_file(self._output_data_dir)
