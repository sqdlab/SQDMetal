from SQDMetal.PALACE.Model import PALACE_Simulation_Base
from SQDMetal.PALACE.SQDGmshRenderer import Palace_Gmsh_Renderer
import matplotlib.pyplot as plt
import numpy as np
import json
import os
import gmsh

class PALACE_Frequency_Driven_Simulation(PALACE_Simulation_Base):

    #Class Variables
    default_user_options = {
                 "mesh_refinement":  0,
                 "dielectric_material": "silicon",
                 "min_freq": 6, #GHz
                 "max_freq": 9, #GHz
                 "freq_step": 0.1, #GHz
                 "save_step": 0,
                 "solver_order": 2,
                 "solver_maxits": 100,
                 "sim_memory": '300G',
                 "sim_time": '8:00:00',
    }

    # Parent Directory path
    HPC_parent_simulation_dir = "C:/PALACE_Simulations/"
    simPC_parent_simulation_dir = "/home/experiment/PALACE/Simulations/input"

    #constructor
    def __init__(self, name, Gmsh_Renderer, mode, launchpad_list, user_options = default_user_options, createFiles = False):
        self.name = name
        self.Gmsh_Renderer = Gmsh_Renderer
        self.mode = mode
        self.launchpad_list = launchpad_list
        self.user_options = user_options
        self.createFiles = createFiles
        
        
    def run(self):
        pass  


    def prepare_simulation(self):
        '''set-up the simulation'''

        #add lumped ports on launch pads
        for _,launchpad in enumerate(self.launchpad_list):
            self.Gmsh_Renderer.add_ports_on_launchpad(launchpad)

        #prepare design by converting shapely geometries to Gmsh geometries
        self.Gmsh_Renderer._prepare_design('frequency_driven_simulation')

        if self.createFiles == True:
            #create directory to store simulation files
            self._create_directory(self.name)

            if self.mode == 'HPC':
                #create batch file
                self.create_batch_file()

            #create config file
            self.create_config_file(self.Gmsh_Renderer)

        #plot design in gmsh gui
        self.Gmsh_Renderer.view_design_components()

        #create mesh
        self.Gmsh_Renderer._intelligent_mesh('frequency_driven_simulation')

        #save
        self._save_mesh()



    def _check_simulation_mode(self):
        '''method to check the type of simualtion being run and return
            the parent directory to store the simulation files'''

        parent_simulation_dir = None

        if self.mode == "HPC":
            parent_simulation_dir = self.HPC_parent_simulation_dir
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


    def _save_mesh(self):

        if self.createFiles == True:
            parent_simulation_dir = self._check_simulation_mode()

            # file_name
            file_name = self.name + "/" + self.name + ".msh"
    
            # Path
            path = os.path.join(parent_simulation_dir, file_name)
            gmsh.write(path)


    def create_config_file(self, Palace_Gmsh_Renderer_object):
        '''create the configuration file which specifies the simulation type and the parameters'''    


        PGR = Palace_Gmsh_Renderer_object

        #Define python dictionary to convert to json file
        config = {
            "Problem":
            {
                "Type": "Driven",
                "Verbose": 2,
                "Output": "/scratch/user/uqdsomm1/Palace-Project2/outputFiles/" + self.name
            },
            "Model":
            {
                "Mesh": "/scratch/user/uqdsomm1/Palace-Project2/inputFiles/"  + self.name + "/" + self.name + ".msh",
                "L0": 1.0e-3,  # mm
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
                        "Attributes": [PGR.config_air_box],  # Air
                        "Permeability": 1.0,
                        "Permittivity": 1.0,
                        "LossTan": 0.0
                    },
                    {
                        "Attributes": [PGR.config_dielectric_base],  # Dielectric
                        "Permeability": 1.0,
                        "Permittivity": 11.7,
                        "LossTan": 1.2e-5
                    }
                ]
            },
            "Boundaries":
            {
                "PEC":
                {
                    "Attributes": PGR.config_metals_rf,  # Metal trace
                },
                    "Absorbing":
                {
                    "Attributes": [PGR.config_far_field],
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
                        "Attributes": [PGR.config_ports[1]],
                        "Direction": "+Y"
                        },
                        {
                        "Attributes": [PGR.config_ports[0]],
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
                        "Attributes": [PGR.config_ports[2]],
                        "Direction": "+Y"
                        },
                        {
                        "Attributes": [PGR.config_ports[3]],
                        "Direction": "-Y"
                        }
                    ]
                }
                ]
            },
            "Solver":
            {
                "Order": self.user_options["solver_order"],
                "Driven":
                {
                    "MinFreq": self.user_options["min_freq"],  # GHz
                    "MaxFreq": self.user_options["max_freq"],  # GHz
                    "FreqStep": self.user_options["freq_step"],  # GHz
                    "SaveStep": self.user_options["save_step"] #Controls how often, in number of frequency steps, to save computed fields to disk for visualization with ParaView
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
        
        sbatch = {
                "header": "#!/bin/bash --login",
                "job_name": "#SBATCH --job-name=" + self.name,
                "output_loc": "#SBATCH --output=" + self.name + ".out",
                "error_out": "#SBATCH --error=" + self.name + ".err",
                "partition": "#SBATCH --partition=general",
                "nodes": "#SBATCH --nodes=1",
                "tasks": "#SBATCH --ntasks-per-node=10",
                "cpus": "#SBATCH --cpus-per-task=1",
                "memory": "#SBATCH --mem=" + self.user_options["sim_memory"],
                "time": "#SBATCH --time=" + self.user_options['sim_time'],
                "account": "#SBATCH --account=a_fedorov",
                "foss": "module load foss/2021a",
                "cmake": "module load cmake/3.20.1-gcccore-10.3.0",
                "pkgconfig": "module load pkgconfig/1.5.4-gcccore-10.3.0-python",
                "run_command": "srun /scratch/user/uqdsomm1/Palace-Project2/palace/build/bin/palace-x86_64.bin " +
                            "/scratch/user/uqdsomm1/Palace-Project2/inputFiles/" + self.name + "/" + self.name + ".json"
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
