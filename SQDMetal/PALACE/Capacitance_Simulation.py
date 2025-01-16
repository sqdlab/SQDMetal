from SQDMetal.PALACE.Model import PALACE_Model
from SQDMetal.COMSOL.Model import COMSOL_Model
from SQDMetal.COMSOL.SimCapacitance import COMSOL_Simulation_CapMats
from SQDMetal.Utilities.Materials import Material
import matplotlib.pyplot as plt
import numpy as np
import json
import os
import io
import gmsh
from SQDMetal.PALACE.Utilities.GMSH_Geometry_Builder import GMSH_Geometry_Builder
from SQDMetal.PALACE.Utilities.GMSH_Mesh_Builder import GMSH_Mesh_Builder

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
                 "mesh_max": 100e-6,
                 "mesh_min": 10e-6,
                 "taper_dist_min": 30e-6,
                 "taper_dist_max": 200e-6,              
                 "gmsh_dist_func_discretisation": 130,
                 "HPC_Parameters_JSON": "",
                 "fuse_threshold": 1e-9
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
            ggb = GMSH_Geometry_Builder(self.metal_design, self.user_options['fillet_resolution'])
            gmsh_render_attrs = ggb.construct_geometry_in_GMSH(self._metallic_layers, self._ground_plane, [], self._fine_meshes, self.user_options["fuse_threshold"])
            #
            gmb = GMSH_Mesh_Builder(gmsh_render_attrs['fine_mesh_elems'], self.user_options)
            gmb.build_mesh()

            if self.create_files == True:
                #create directory to store simulation files
                self._create_directory(self.name)

                #create config file
                self.create_config_file(gmsh_render_attrs = gmsh_render_attrs)

                #create batch file
                # if self.mode == 'HPC':
                #     self.create_batch_file()

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



    def create_config_file(self, **kwargs):
        '''create the configuration file which specifies the simulation type and the parameters'''    

        if self.meshing == 'GMSH':

            #GMSH renderer object needed to get boundary conditions for config file
            gmsh_render_attrs = kwargs['gmsh_render_attrs']

            #GMSH config file variables
            material_air = gmsh_render_attrs['air_box']
            material_dielectric = gmsh_render_attrs['dielectric']
            far_field = gmsh_render_attrs['far_field']
            
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
        if self._output_subdir == "":
            self.set_local_output_subdir("", False)
        filePrefix = self.hpc_options["input_dir"]  + self.name + "/" if self.hpc_options["input_dir"] != "" else ""
        self._mesh_name = filePrefix + self.name + file_ext
        post_procs = []
        for x in Terminal:
            post_procs.append({"Index":x["Index"], "Attributes":x["Attributes"], "Type": "Electric"})
        config = {
                    "Problem":
                    {
                        "Type": "Electrostatic",
                        "Verbose": 2,
                        "Output": filePrefix  + "outputFiles"
                    },
                    "Model":
                    {
                        "Mesh":  self._mesh_name,
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
                                "LossTan": dielectric.loss_tangent
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
                        "SurfaceFlux": post_procs,
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
        self.path_mesh = os.path.join(simulation_dir, self._mesh_name)
        self._sim_config = file
        self.set_local_output_subdir(self._output_subdir)
