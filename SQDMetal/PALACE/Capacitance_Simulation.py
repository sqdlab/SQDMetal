from SQDMetal.PALACE.Model import PALACE_Simulation_Base
from SQDMetal.PALACE.SQDGmshRenderer import Palace_Gmsh_Renderer
import matplotlib.pyplot as plt
import numpy as np
import json
import os
import io
import gmsh

class PALACE_Capacitance_Simulation(PALACE_Simulation_Base):

    #Class Variables
    default_user_options = {
                 "mesh_refinement":  0,
                 "dielectric_material": "silicon",
                 "solns_to_save": 3,
                 "solver_order": 2,
                 "solver_tol": 1.0e-8,
                 "solver_maxits": 100,
                 "sim_memory": '200G',
                 "sim_time": '1:00:00',
    }

    # Parent Directory path
    HPC_parent_simulation_dir = "C:/PALACE_Simulations/"
    simPC_parent_simulation_dir = "/home/experiment/PALACE/Simulations/input"

    #constructor
    def __init__(self, name, Gmsh_Renderer, mode, user_options = default_user_options, createFiles = False):
        self.name = name
        self.Gmsh_Renderer = Gmsh_Renderer
        self.mode = mode
        self.user_options = user_options
        self.createFiles = createFiles


    def run(self):
        pass  


    def prepare_simulation(self):
        '''set-up the simulation'''

        #prepare design in Gmsh
        self.Gmsh_Renderer._prepare_design('capacitance_simulation')
        
        #Create simulation files
        if self.createFiles == True:
            #create directory to store simulation files
            self._create_directory(self.name)

            #if using HPC create batch file
            if self.mode == 'HPC':
                self.create_batch_file()

            #create config file
            self.create_config_file(self.Gmsh_Renderer)

        #open gmsh gui to show design
        self.Gmsh_Renderer.view_design_components()

        #create mesh
        self.Gmsh_Renderer._intelligent_mesh('capacitance_simulation')

        #save mesh
        self._save_mesh()



    def _check_simulation_mode(self):
        
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

        if self.createFiles == True:
            parent_simulation_dir = self._check_simulation_mode()

            # Directory
            directory = directory_name
    
            # Path
            path = os.path.join(parent_simulation_dir, directory)
    
            # Create the directory
            os.mkdir(path)
            print("Directory '% s' created" % directory)


    def _save_mesh(self):

        parent_simulation_dir = self._check_simulation_mode()

        # file_name
        file_name = self.name + "/" + self.name + ".msh"
  
        # Path
        path = os.path.join(parent_simulation_dir, file_name)
        gmsh.write(path)
        


    def create_config_file(self, Palace_Gmsh_Renderer_object):
        '''create the configuration file which specifies the simulation type and the parameters'''    

        PGR = Palace_Gmsh_Renderer_object

        Terminal = []
        for i,value in enumerate(PGR.config_metals_cap):
            metal = {"Index": i+1,
                     "Attributes": [value]}
            Terminal.append(metal)

        #Define python dictionary to convert to json file
        config = {
                    "Problem":
                    {
                        "Type": "Electrostatic",
                        "Verbose": 2,
                        "Output": "/scratch/user/uqdsomm1/Palace-Project2/outputFiles/" + self.name
                    },
                    "Model":
                    {
                        "Mesh": "/scratch/user/uqdsomm1/Palace-Project2/inputFiles/" + self.name + "/" + self.name + ".msh",
                        "L0": 0.001,  # mm
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
                        "Ground":
                        {
                        "Attributes": [PGR.config_far_field]
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
                        "Order": 2,
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
        with open(file, "w+", newline='\n') as f:
            for value in sbatch.values():
                f.write('{}\n'.format(value))
