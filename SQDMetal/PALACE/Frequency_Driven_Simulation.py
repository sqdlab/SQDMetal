from SQDMetal.PALACE.Model import PALACE_Model_RF_Base
from SQDMetal.PALACE.SQDGmshRenderer import Palace_Gmsh_Renderer
from SQDMetal.COMSOL.Model import COMSOL_Model
from SQDMetal.COMSOL.SimRFsParameter import COMSOL_Simulation_RFsParameters
from SQDMetal.Utilities.Materials import Material
import matplotlib.pyplot as plt
import numpy as np
import json
import os
import gmsh
import pandas as pd

class PALACE_Driven_Simulation(PALACE_Model_RF_Base):

    #Class Variables
    default_user_options = {
                 "fillet_resolution": 4,
                 "mesh_refinement":  0,
                 "dielectric_material": "silicon",
                 "solns_to_save": 4,
                 "solver_order": 2,
                 "solver_tol": 1.0e-8,
                 "solver_maxits": 300,
                 "mesh_max": 100e-3,
                 "mesh_min": 10e-3,
                 "mesh_sampling": 120,
                 "comsol_meshing": "Extremely fine",
                 "HPC_Parameters_JSON": ""
                }

    #constructor
    def __init__(self, name, metal_design, sim_parent_directory, mode, meshing, user_options = {}, 
                 view_design_gmsh_gui = False, create_files = False):
        self.name = name
        self.metal_design = metal_design
        self.sim_parent_directory = sim_parent_directory
        self.mode = mode
        self.meshing = meshing
        self.user_options = {}
        for key in PALACE_Driven_Simulation.default_user_options:
            self.user_options[key] = user_options.get(key, PALACE_Driven_Simulation.default_user_options[key])
        self.view_design_gmsh_gui = view_design_gmsh_gui
        self.create_files = create_files
        self._ports = []
        super().__init__(meshing, mode, user_options)
        self.freqs = None


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
                ports[cur_port['port_name'] + 'a'] = list(comsol._model.java.component("comp1").physics("emw").feature(f"lport{m}").selection().entities())[1]
                ports[cur_port['port_name'] + 'b'] = list(comsol._model.java.component("comp1").physics("emw").feature(f"lport{m}").selection().entities())[0]

            #define length scale
            l0 = 1

            #file extension
            file_ext = '.mphbin'

        #get material parameters
        dielectric = Material(self.user_options["dielectric_material"])

        #Process Ports
        config_ports = []
        for m, cur_port in enumerate(self._ports):
            port_name, vec_perps = cur_port['port_name'], cur_port['vec_CPW2GND_1']
            leDict = {
                    "Index": m+1,
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
                }
            if 'impedance_R' in cur_port:
                leDict['R'] = cur_port['impedance_R']
            if 'impedance_L' in cur_port:
                leDict['L'] = cur_port['impedance_L']
            if 'impedance_C' in cur_port:
                leDict['C'] = cur_port['impedance_C']
            config_ports.append(leDict)
        config_ports[0]["Excitation"] = True

        if isinstance(self.freqs, tuple):
            fStart,fStop,fStep = self.freqs
        else:
            print("Warning: Frequencies have not been set properly...")
            fStart,fStop,fStep = (1,10,1)

        #Define python dictionary to convert to json file
        if self._output_subdir == "":
            self.set_local_output_subdir("", False)
        filePrefix = self._get_folder_prefix()
        config = {
            "Problem":
            {
                "Type": "Driven",
                "Verbose": 2,
                "Output": self._output_dir
            },
            "Model":
            {
                "Mesh":  filePrefix + self.name + file_ext,
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
                "Driven":
                {
                    "MinFreq": fStart,  # starting freequency
                    "MaxFreq": fStop,  # end frequency
                    "FreqStep": fStep,  # step size
                    "SaveStep": self.user_options["solns_to_save"] # Number of frequency steps for computed field modes to save to disk for visualization with ParaView
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
        self._sim_config = file
        self.set_local_output_subdir(self._output_subdir)
    
    def set_freq_values(self, startGHz, endGHz, stepGHz):
        self.freqs = (startGHz, endGHz, stepGHz)
        if self._sim_config != "":
            with open(self._sim_config, "r") as f:
                config_json = json.loads(f.read())
            config_json['Solver']['Driven']['MinFreq'] = startGHz
            config_json['Solver']['Driven']['MaxFreq'] = endGHz
            config_json['Solver']['Driven']['FreqStep'] = stepGHz
            with open(self._sim_config, "w") as f:
                json.dump(config_json, f, indent=2)

    def retrieve_data(self):
        raw_data = pd.read_csv(self._output_data_dir + '/port-S.csv')
        headers = raw_data.columns
        raw_data = raw_data.to_numpy().T

        ret_data = {'freqs' : raw_data[0]*1e9}
        for m in range(int((raw_data.shape[0]-1)/2)):
            #Header will be like: |S[1][1]| (dB)
            key = headers[2*m+1].strip().split('|')[1].replace('[','').replace(']','')
            amp = raw_data[2*m+1]
            phs = raw_data[2*m+2]
            sParam = 10**(amp/20)*np.exp(1j*phs/180*np.pi)
            ret_data[key] = sParam
        
        return ret_data
