from SQDMetal.PALACE.Model import PALACE_Model
from SQDMetal.PALACE.SQDGmshRenderer import Palace_Gmsh_Renderer
from SQDMetal.COMSOL.Model import COMSOL_Model
from SQDMetal.COMSOL.SimCapacitance import COMSOL_Simulation_CapMats
from SQDMetal.Utilities.Materials import Material
from SQDMetal.Utilities.QUtilities import QUtilities
import matplotlib.pyplot as plt
import numpy as np
import json
import os
import io
import gmsh

class PALACE_Capacitance_Simulation(PALACE_Model):

    #Class Variables
    default_user_options = {
                 "fillet_resolution": 4,
                 "mesh_refinement":  0,
                 "dielectric_material": "silicon",
                 "solns_to_save": 3,
                 "solver_order": 2,
                 "solver_tol": 1.0e-8,
                 "solver_maxits": 100,
                 "mesh_max": 120e-3,                                
                 "mesh_min": 10e-3,                                 
                 "mesh_sampling": 130,                              
                 "fillet_resolution":12,
                 "HPC_Parameters_JSON": ""
    }

    # Parent Directory path
    HPC_parent_simulation_dir = "C:/PALACE_Simulations/"
    simPC_parent_simulation_dir = "/home/experiment/PALACE/Simulations/input"

    #constructor
    def __init__(self, name, metal_design, sim_parent_directory, mode, meshing, user_options = {}, 
                                    view_design_gmsh_gui = False, create_files = False):
        self.name = name
        self.metal_design = metal_design
        self.sim_parent_directory = sim_parent_directory
        self.mode = mode
        self.meshing = meshing
        self.user_options = {}
        for key in PALACE_Capacitance_Simulation.default_user_options:
            self.user_options[key] = user_options.get(key, PALACE_Capacitance_Simulation.default_user_options[key])
        self.view_design_gmsh_gui = view_design_gmsh_gui
        self.create_files = create_files
        super().__init__(meshing, mode, user_options)


    def prepare_simulation(self):
        '''set-up the simulation'''

        if self.meshing == 'GMSH':

            #Create the gmsh renderer to convert qiskit metal geomentry to gmsh geometry
            pgr = Palace_Gmsh_Renderer(self.metal_design)

            #prepare design by converting shapely geometries to Gmsh geometries
            gmsh_render_attrs = pgr._prepare_design(self._metallic_layers, self._ground_plane, [], self.user_options['fillet_resolution'], 'capacitance_simulation')

            if self.create_files == True:
                #create directory to store simulation files
                self._create_directory(self.name)

                #create config file
                self.create_config_file(gmsh_render_attrs = gmsh_render_attrs)

                #create batch file
                if self.mode == 'HPC':
                    self.create_batch_file()
                
                #create mesh
                # pgr._intelligent_mesh('capacitance_simulation', 
                #                 min_size = self.user_options['mesh_min'], 
                #                 max_size = self.user_options['mesh_max'], 
                #                 mesh_sampling = self.user_options['mesh_sampling'])
                pgr.fine_mesh(self._fine_meshes)

                self._save_mesh_gmsh()

            if self.view_design_gmsh_gui == True:
                #plot design in gmsh gui
                pgr.view_design_components()
            
                   
        if self.meshing == 'COMSOL':

            #initialise the COMSOL engine
            COMSOL_Model.init_engine()
            cmsl = COMSOL_Model('res')

            #Create COMSOL capacitance sim object
            simCapMats = COMSOL_Simulation_CapMats(cmsl)
            cmsl.initialize_model(self.metal_design, [simCapMats], bottom_grounded = True, resolution = 10)

            #Add metallic layers
            cmsl.add_metallic(1, threshold=1e-10, fuse_threshold=1e-10)
            cmsl.add_ground_plane(threshold=1e-10)
            #cmsl.fuse_all_metals()

            #build model
            cmsl.build_geom_mater_elec_mesh(mesh_structure = self.user_options["comsol_meshing"])

            #plot model
            cmsl.plot()

            #save comsol file
            #cmsl.save(self.name)

            if self.create_files == True:
                #create directory to store simulation files
                self._create_directory(self.name)

                #create config file
                self.create_config_file(comsol_obj = cmsl, simCap_object = simCapMats)

                #create batch file
                if self.mode == 'HPC':
                    self.create_batch_file()

            #save mesh
            self._save_mesh_comsol(comsol_obj = cmsl)

    def _create_directory(self, directory_name):
        '''create a directory to hold the simulation files'''

        if self.create_files == True:
            parent_simulation_dir = self._check_simulation_mode()

            # Directory
            directory = directory_name
    
            # Path
            path = os.path.join(parent_simulation_dir, directory)
    
            # Create the directory
            if not os.path.exists(path):
                os.mkdir(path)
                print("Directory '% s' created" % directory)



    def _save_mesh_gmsh(self):

        parent_simulation_dir = self._check_simulation_mode()

        # file_name
        file_name = self.name + "/" + self.name + ".msh"
  
        # Path
        path = os.path.join(parent_simulation_dir, file_name)
        gmsh.write(path)
        


    def _save_mesh_comsol(self, comsol_obj):
        '''function used to save the comsol mesh file'''
        if self.create_files == True:
            parent_simulation_dir = self._check_simulation_mode()

            # file_name
            file_name = self.name + "/" + self.name + ".mphbin"
    
            # Path
            path = os.path.join(parent_simulation_dir, file_name)

            #COMSOL export commands
            comsol_obj._model.java.component("comp1").mesh("mesh1").export().set("filename", path)
            comsol_obj._model.java.component("comp1").mesh("mesh1").export(path)

    def fine_mesh_in_rectangle(self, x1, y1, x2, y2, **kwargs):
        self._fine_meshes.append({
            'type': 'box',
            'x_bnds': (x1 * 1e3, x2 * 1e3), #Get it into mm...
            'y_bnds': (y1 * 1e3, y2 * 1e3), #Get it into mm...
            'mesh_sampling': kwargs.get('mesh_sampling', self.user_options['mesh_sampling']),
            'min_size': kwargs.get('min_size', self.user_options['mesh_min']),
            'max_size': kwargs.get('max_size', self.user_options['mesh_max'])
        })

    def fine_mesh_along_path(self, dist_resolution, qObjName, trace_name='', **kwargs):
        leUnits = QUtilities.get_units(self.metal_design)
        lePath = QUtilities.calc_points_on_path(dist_resolution/leUnits, self.metal_design, qObjName, trace_name)[0] * leUnits
        self._fine_meshes.append({
            'type': 'path',
            'path': lePath * 1e3, #Get it into mm...
            'mesh_sampling': kwargs.get('mesh_sampling', self.user_options['mesh_sampling']),
            'min_size': kwargs.get('mesh_min', self.user_options['mesh_min']),
            'max_size': kwargs.get('mesh_max', self.user_options['mesh_max'])
        })

    def create_config_file(self, **kwargs):
        '''create the configuration file which specifies the simulation type and the parameters'''    

        if self.meshing == 'GMSH':

            #GMSH renderer object needed to get boundary conditions for config file
            gmsh_render_attrs = kwargs['gmsh_render_attrs']

            #GMSH config file variables
            material_air = [gmsh_render_attrs['air_box']]
            material_dielectric = [gmsh_render_attrs['dielectric']]
            far_field = [gmsh_render_attrs['far_field']]
            
            #metals to compute capacitances for
            Terminal = []
            for i,value in enumerate(gmsh_render_attrs['metals']):
                metal = {"Index": i+1,
                        "Attributes": [value]}
                Terminal.append(metal)

            #define length scale
            l0 = 1e-3

            #file extension
            file_ext = '.msh'
        
        if self.meshing == 'COMSOL':

            #COMSOL object and RF Sim object needed to get boundary conditions for config file
            comsol = kwargs['comsol_obj']  
            sCapMats = kwargs['simCap_object']

            #COMSOL config file variables
            material_air = list(comsol._model.java.component("comp1").material("Vacuum").selection().entities())
            material_dielectric = list(comsol._model.java.component("comp1").material("Substrate").selection().entities())
            far_field = list(comsol._model.java.component("comp1").physics(sCapMats.phys_es).feature("gnd1").selection().entities())

            #metals to compute capacitances for
            Terminal = []
            for i,term_name in enumerate(sCapMats.terminal_names):
                term_assigned_value = list(comsol._model.java.component("comp1").physics(sCapMats.phys_es).feature(term_name).selection().entities())
                metal = {"Index": i+1,
                        "Attributes": term_assigned_value}
                Terminal.append(metal)

            #define length scale
            l0 = 1

            #file extension
            file_ext = '.mphbin'

        #get material parameters
        dielectric = Material(self.user_options["dielectric_material"])

        
        #Define python dictionary to convert to json file
        filePrefix = self.hpc_options["input_dir"]  + self.name + "/" if self.hpc_options["input_dir"] != "" else ""
        config = {
                    "Problem":
                    {
                        "Type": "Electrostatic",
                        "Verbose": 2,
                        "Output": filePrefix  + "outputFiles"
                    },
                    "Model":
                    {
                        "Mesh":  filePrefix + self.name + file_ext,
                        "L0": l0,  # mm
                        "Refinement":
                        {
                            "UniformLevels": self.user_options["mesh_refinement"]
                        },
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
                                "LossTan": 1.2e-5
                            }
                        ]
                    },
                    "Boundaries":
                    {
                        "Ground":
                        {
                        "Attributes": far_field
                        },
                        "Terminal":
                            Terminal,
                        "Postprocessing":  # Capacitance from charge instead of energy
                        {
                        "Capacitance":
                            Terminal,
                        }
                    },
                    "Solver":
                    {
                        "Order": self.user_options["solver_order"],
                        "Electrostatic":
                        {
                        "Save": self.user_options["solns_to_save"]
                        },
                        "Linear":
                        {
                        "Type": "BoomerAMG",
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
        self._sim_config = file
