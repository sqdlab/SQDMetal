import os, subprocess
import json
import gmsh

class Simulation_Files_Builder:

    def __init__(self, name, user_options, sim_config, hpc_options):
        self.name = name
        self.user_options = user_options
        self.sim_config = sim_config
        self.hpc_options = hpc_options
    
    def create_simulation_files(self):
        
        #create directory to store simulation files
        self._create_directory()

        #save mesh to new directory
        self._save_mesh_gmsh()

        #write sim_config to json file and save to new directory
        self._save_config_file_as_json()

        if self.hpc_options:
            #create hpc batch file for simulations using the
            self._create_hpc_batch_file()

        print('Simulation files created.')


    def _save_config_file_as_json(self):
        
        #simulation file name
        json_file_name = self.name + "/" + self.name + ".json"
        
        #save to created directory
        file = os.path.join(self.user_options['sim_directory'], json_file_name)

        #write to file
        with open(file, "w+") as f:
            json.dump(self.sim_config, f, indent=2)


    def _create_directory(self):
        '''create a directory to hold the simulation files'''
  
        # Path
        path = os.path.join(self.user_options['sim_directory'], self.name)
  
        # Create the directory
        if not os.path.exists(path):
            os.makedirs(path)
            print("Directory '% s' created at '% s'." % (self.name, path))
        else:
            Exception('Path already exists. Create a new path.')


    def _save_mesh_gmsh(self):
        '''function used to save the gmsh mesh file'''

        mesh_file_name = self.name + "/" + self.name + ".msh"
        path = os.path.join(self.user_options['sim_directory'], mesh_file_name)
        gmsh.write(path)


    def _create_hpc_batch_file(self):
        
    
        #note: I have disabled naming the output file by setting '# SBATCH' instead of '#SBATCH' 
        #so I can get the slurm job number to use for testing
        sbatch = {
                "header": "#!/bin/bash --login",
                "job_name": "#SBATCH --job-name=" + self.name,
                "output_loc": "# SBATCH --output=" + self.name + ".out",
                "error_out": "#SBATCH --error=" + self.name + ".err",
                "partition": "#SBATCH --partition=general",
                "nodes": "#SBATCH --nodes=" + self.hpc_options["hpc_nodes"],
                "tasks": "#SBATCH --ntasks-per-node=" + self.hpc_options['cpus_per_node'],
                "cpus": "#SBATCH --ntasks-per-core=1",
                "memory": "#SBATCH --mem=" + self.hpc_options["sim_memory"],
                "time": "#SBATCH --time=" + self.hpc_options['sim_time'],
                "account": "#SBATCH --account=" + self.hpc_options['account_name'],
                "easy_build": "module use /scratch/project_mnt/palace-sqdlab/Palace-Project-NEW/EasyBuild/modules/all",
                "foss": "module load foss/2023a",
                "cmake": "module load cmake/3.26.3-gcccore-12.3.0",
                #"pkgconfig": "module load pkgconfig/1.5.5-gcccore-12.3.0-python",
                "run_command": f"srun {self.hpc_options['palace_location']} " + self.hpc_options['input_files_location'] + self.name + "/" + self.name + ".json"
        }

        #simulation file name
        file_name = self.name + "/" + self.name + ".sbatch"
        
        #save to created directory
        file = os.path.join(self.user_options['sim_directory'], file_name)

        #write sbatch dictionary to file
        with open(file, "w+", newline = '\n') as f:
            for value in sbatch.values():
                f.write('{}\n'.format(value))

