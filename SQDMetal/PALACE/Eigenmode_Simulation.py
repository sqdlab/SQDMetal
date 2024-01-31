from SQDMetal.PALACE.Model import PALACE_Model_RF_Base
from SQDMetal.PALACE.SQDGmshRenderer import Palace_Gmsh_Renderer
from SQDMetal.COMSOL.Model import COMSOL_Model
from SQDMetal.COMSOL.SimRFsParameter import COMSOL_Simulation_RFsParameters
from SQDMetal.Utilities.Materials import Material
from SQDMetal.Utilities.QUtilities import QUtilities
import matplotlib.pyplot as plt
import numpy as np
import json
import os
import gmsh

class PALACE_Eigenmode_Simulation(PALACE_Model_RF_Base):

    #Class Variables
    default_user_options = {
                 "fillet_resolution": 4,
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
    def __init__(self, name, metal_design, sim_parent_directory, mode, meshing, user_options = {}, 
                 view_design_gmsh_gui = False, create_files = False):
        self.name = name
        self.metal_design = metal_design
        self.sim_parent_directory = sim_parent_directory
        self.mode = mode
        self.user_options = {}
        for key in PALACE_Eigenmode_Simulation.default_user_options:
            self.user_options[key] = user_options.get(key, PALACE_Eigenmode_Simulation.default_user_options[key])
        self.view_design_gmsh_gui = view_design_gmsh_gui
        self.create_files = create_files
        self._ports = []
        super().__init__(meshing)
        
        
    def run(self):
        pass  


    def create_config_file(self, **kwargs):
        '''create the configuration file which specifies the simulation type and the parameters'''    

        if self.meshing == 'GMSH':

            #GMSH renderer object needed to get boundary conditions for config file
            gmsh_render_attrs = kwargs['gmsh_render_attrs']

            #GMSH config file variables
            material_air = [gmsh_render_attrs['air_box']]
            material_dielectric = [gmsh_render_attrs['dielectric']]
            PEC_metals = gmsh_render_attrs['metals']
            far_field = [gmsh_render_attrs['far_field']]
            ports = gmsh_render_attrs['ports']

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
            ports = {}
            for m, cur_port in enumerate(self._ports):
                port_name, qObjName, launchesA, launchesB, vec_perps = cur_port
                ports[port_name + 'a'] = list(comsol._model.java.component("comp1").physics("emw").feature(f"lport{m}").selection().entities())[1]
                ports[port_name + 'b'] = list(comsol._model.java.component("comp1").physics("emw").feature(f"lport{m}").selection().entities())[0]

            #define length scale
            l0 = 1

            #file extension
            file_ext = '.mphbin'

        #get material parameters
        dielectric = Material(self.user_options["dielectric_material"])

        #Process Ports
        config_ports = []
        for m, cur_port in enumerate(self._ports):
            port_name, qObjName, launchesA, launchesB, vec_perps = cur_port
            config_ports.append({
                    "Index": m+1,
                    "R": 50.00,  # Ω, 2-element uniform
                    "Elements":
                    [
                        {
                        "Attributes": [ports[port_name + 'a']],
                        "Direction": vec_perps[0]
                        },
                        {
                        "Attributes": [ports[port_name + 'b']],
                        "Direction": vec_perps[1]
                        }
                    ]
                })
        config_ports[0]["Excitation"] = True

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
                "LumpedPort": config_ports
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
