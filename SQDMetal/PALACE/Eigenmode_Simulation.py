from SQDMetal.PALACE.Model import PALACE_Simulation_Base
from SQDMetal.PALACE.SQDGmshRenderer import Palace_Gmsh_Renderer
from SQDMetal.COMSOL.Model import COMSOL_Model
from SQDMetal.COMSOL.SimRFsParameter import COMSOL_Simulation_RFsParameters
from SQDMetal.Utilities.Materials import Material
import matplotlib.pyplot as plt
import numpy as np
import json
import os
import gmsh

class PALACE_Eigenmode_Simulation(PALACE_Simulation_Base):

    #Class Variables
    default_user_options = {
                 "mesh_refinement":  0,
                 "dielectric_material": "silicon",
                 "starting_freq": 5.5,
                 "number_of_freqs": 5,
                 "solns_to_save": 5,
                 "solver_order": 2,
                 "solver_tol": 1.0e-8,
                 "solver_maxits": 100,
                 "mesh_max": 100e-3,
                 "mesh_min": 10e-3,
                 "mesh_sampling": 120,
                 "sim_memory": '300G',
                 "sim_time": '14:00:00',
                 "HPC_nodes": '4',
                }

    #Parent Directory path
    simPC_parent_simulation_dir = "/home/experiment/PALACE/Simulations/input"

    #constructor
    def __init__(self, name, metal_design, sim_parent_directory, mode, meshing, launchpad_list, user_options = {}, 
                 view_design_gmsh_gui = False, create_files = False):
        self.name = name
        self.metal_design = metal_design
        self.sim_parent_directory = sim_parent_directory
        self.mode = mode
        self.meshing = meshing
        self.launchpad_list = launchpad_list
        self.user_options = {}
        for key in default_user_options:
            self.user_options[key] = user_options.get(key, default_user_options[key])
        self.view_design_gmsh_gui = view_design_gmsh_gui
        self.create_files = create_files
        
        
    def run(self):
        pass  


    def prepare_simulation(self):
        '''set-up the simulation'''
        
        if self.meshing == 'GMSH':

            #Create the gmsh renderer to convert qiskit metal geomentry to gmsh geometry
            pgr = Palace_Gmsh_Renderer(self.metal_design)

            #add lumped ports on launch pads
            for _,launchpad in enumerate(self.launchpad_list):
                pgr.add_ports_on_launchpad(launchpad)

            #prepare design by converting shapely geometries to Gmsh geometries
            pgr._prepare_design('eigenmode_simulation')

            if self.create_files == True:
                #create directory to store simulation files
                self._create_directory(self.name)

                #create config file
                self.create_config_file(Palace_Gmsh_Renderer_object = pgr)

                #create batch file
                if self.mode == 'HPC':
                    self.create_batch_file()
                
                #create mesh
                pgr._intelligent_mesh('eigenmode_simulation', 
                                min_size = self.user_options['mesh_min'], 
                                max_size = self.user_options['mesh_max'], 
                                mesh_sampling = self.user_options['mesh_sampling'])

                self._save_mesh_gmsh()

            if self.view_design_gmsh_gui == True:
                #plot design in gmsh gui
                pgr.view_design_components()
                    
        if self.meshing == 'COMSOL':

            #initialise the COMSOL engine
            COMSOL_Model.init_engine()
            cmsl = COMSOL_Model('res')

            #Create COMSOL RF sim object
            sim_sParams = COMSOL_Simulation_RFsParameters(cmsl, adaptive = 'None')
            cmsl.initialize_model(self.metal_design, [sim_sParams], bottom_grounded = True, resolution = 10)

            #Add metallic layers
            cmsl.add_metallic(1, threshold=1e-10, fuse_threshold=1e-10)
            cmsl.add_ground_plane(threshold=1e-10)
            #cmsl.fuse_all_metals()

            #assign ports
            for _,launchpad in enumerate(self.launchpad_list):
                sim_sParams.create_port_CPW_on_Launcher(launchpad.name)
            
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
                self.create_config_file(comsol_obj = cmsl, simRF_object = sim_sParams)

                #create batch file
                if self.mode == 'HPC':
                    self.create_batch_file()

            #save mesh
            self._save_mesh_comsol(comsol_obj = cmsl)


    def _check_simulation_mode(self):
        '''method to check the type of simualtion being run and return
            the parent directory to store the simulation files'''

        parent_simulation_dir = None

        if self.mode == "HPC":
            parent_simulation_dir = self.sim_parent_directory
        elif self.mode == "simPC":
            parent_simulation_dir = self.simPC_parent_simulation_dir
        else:
            Exception('Invalid simulation mode entered.')
        
        return parent_simulation_dir



    def _create_directory(self, directory_name):
        '''create a directory to hold the simulation files'''

        parent_simulation_dir = self._check_simulation_mode()

        # Directory
        directory = directory_name
  
        # Path
        path = os.path.join(parent_simulation_dir, directory)
  
        # Create the directory
        os.mkdir(path)
        print("Directory '% s' created" % directory)


    def _save_mesh_gmsh(self):
        '''function used to save the gmsh mesh file'''
        if self.create_files == True:
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
            PGR = kwargs['Palace_Gmsh_Renderer_object']

            #GMSH config file variables
            material_air = [PGR.config_air_box]
            material_dielectric = [PGR.config_dielectric_base]
            PEC_metals = PGR.config_metals_rf
            far_field = [PGR.config_far_field]
            port1_pos_y = [PGR.config_ports[0]]
            port1_neg_y = [PGR.config_ports[1]]
            port2_pos_y = [PGR.config_ports[3]]
            port2_neg_y = [PGR.config_ports[2]]

            #define length scale
            l0 = 1e-3

            #file extension
            file_ext = '.msh'
        
        if self.meshing == 'COMSOL':

            #COMSOL object and Cap Sim object needed to get boundary conditions for config file
            comsol = kwargs['comsol_obj']  
            sParams_sim = kwargs['simRF_object']

            #COMSOL config file variables
            material_air = list(comsol._model.java.component("comp1").material("Vacuum").selection().entities())
            material_dielectric = list(comsol._model.java.component("comp1").material("Substrate").selection().entities())
            PEC_metals = list(comsol._model.java.component("comp1").material("Metal").selection().entities())
            far_field = list(comsol._model.java.component("comp1").physics(sParams_sim.phys_emw).feature("pec1").selection().entities())
            port1_pos_y = [list(comsol._model.java.component("comp1").physics("emw").feature("lport0").selection().entities())[1]]
            port1_neg_y = [list(comsol._model.java.component("comp1").physics("emw").feature("lport0").selection().entities())[0]]
            port2_pos_y = [list(comsol._model.java.component("comp1").physics("emw").feature("lport1").selection().entities())[1]]
            port2_neg_y = [list(comsol._model.java.component("comp1").physics("emw").feature("lport1").selection().entities())[0]]

            #define length scale
            l0 = 1

            #file extension
            file_ext = '.mphbin'

        #get material parameters
        dielectric = Material(self.user_options["dielectric_material"])

        #Define python dictionary to convert to json file
        config = {
            "Problem":
            {
                "Type": "Eigenmode",
                "Verbose": 2,
                "Output": "/scratch/project/palace-sqdlab/outputFiles/" + self.name
            },
            "Model":
            {
                "Mesh": "/scratch/project/palace-sqdlab/inputFiles/"  + self.name + "/" + self.name + file_ext,
                "L0": l0,  
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
                "PEC":
                {
                    "Attributes": PEC_metals,  # Metal trace
                },
                    "Absorbing":
                {
                    "Attributes": far_field,
                    "Order": 1
                },
                "LumpedPort":
                [
                {
                    "Index": 1,
                    "R": 50.00,  # Ω, 2-element uniform
                    "Excitation": True,
                    "Elements":
                    [
                        {
                        "Attributes": port1_pos_y,
                        "Direction": "+Y"
                        },
                        {
                        "Attributes": port1_neg_y,
                        "Direction": "-Y"
                        }
                    ]
                },
                {
                    "Index": 2,
                    "R": 50.00,
                    "Elements":
                    [
                        {
                        "Attributes": port2_pos_y,
                        "Direction": "+Y"
                        },
                        {
                        "Attributes": port2_neg_y,
                        "Direction": "-Y"
                        }
                    ]
                }
                ]
            },
            "Solver":
            {
                "Order": self.user_options["solver_order"],
                "Eigenmode":
                {
                    "N": self.user_options["number_of_freqs"],  # number of eigenfrequencies
                    "Tol": self.user_options["solver_tol"],  # solver tolerance
                    "Target": self.user_options["starting_freq"],  # GHz - starting point
                    "Save": self.user_options["solns_to_save"] # Number of computed field modes to save to disk for visualization with ParaView
                },
                "Linear":
                {
                    "Type": "SuperLU",
                    "KSPType": "FGMRES",
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
    


    def create_batch_file(self):
        
        #note: I have disabled naming the output file by setting '# SBATCH' instead of '#SBATCH' 
        #so I can get the slurm job number to use for testing

        sbatch = {
                "header": "#!/bin/bash --login",
                "job_name": "#SBATCH --job-name=" + self.name,
                "output_loc": "# SBATCH --output=" + self.name + ".out",
                "error_out": "#SBATCH --error=" + self.name + ".err",
                "partition": "#SBATCH --partition=general",
                "nodes": "#SBATCH --nodes=" + self.user_options["HPC_nodes"],
                "tasks": "#SBATCH --ntasks-per-node=20",
                "cpus": "#SBATCH --cpus-per-task=1",
                "memory": "#SBATCH --mem=" + self.user_options["sim_memory"],
                "time": "#SBATCH --time=" + self.user_options['sim_time'],
                "account": "#SBATCH --account=a_fedorov",
                "foss": "module load foss/2021a",
                "cmake": "module load cmake/3.20.1-gcccore-10.3.0",
                "pkgconfig": "module load pkgconfig/1.5.4-gcccore-10.3.0-python",
                "run_command": "srun /scratch/project/palace-sqdlab/Palace-Project/palace/build/bin/palace-x86_64.bin " +
                            "/scratch/project/palace-sqdlab/inputFiles/" + self.name + "/" + self.name + ".json"
        }
    
        #check simulation mode and return appropriate parent directory 
        parent_simulation_dir = self._check_simulation_mode()

        #create sbatch file name
        sim_file_name = self.name + '.sbatch'

        #destination for config file
        simulation_dir = parent_simulation_dir + str(self.name)

        #save to created directory
        file = os.path.join(simulation_dir, sim_file_name)

        #write sbatch dictionary to file
        with open(file, "w+", newline = '\n') as f:
            for value in sbatch.values():
                f.write('{}\n'.format(value))
