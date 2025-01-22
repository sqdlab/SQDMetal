from SQDMetal.PALACE.Model import PALACE_Model_RF_Base
from SQDMetal.COMSOL.Model import COMSOL_Model
from SQDMetal.COMSOL.SimRFsParameter import COMSOL_Simulation_RFsParameters
from SQDMetal.Utilities.Materials import Material
from SQDMetal.Utilities.QUtilities import QUtilities
from SQDMetal.PALACE.PVDVTU_Viewer import PVDVTU_Viewer
import matplotlib.pyplot as plt
import numpy as np
import json
import os
import gmsh
import pandas as pd

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
                 "mesh_max": 100e-6,
                 "mesh_min": 10e-6,
                 "taper_dist_min": 30e-6,
                 "taper_dist_max": 200e-6,
                 "gmsh_dist_func_discretisation": 120,
                 "comsol_meshing": "Extremely fine",
                 "HPC_Parameters_JSON": "",
                 "fuse_threshold": 1e-9
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
        super().__init__(meshing, mode, user_options)

    def create_config_file(self, **kwargs):
        '''create the configuration file which specifies the simulation type and the parameters'''    

        if self.meshing == 'GMSH':

            #GMSH renderer object needed to get boundary conditions for config file
            gmsh_render_attrs = kwargs['gmsh_render_attrs']

            #GMSH config file variables
            material_air = gmsh_render_attrs['air_box']
            material_dielectric = gmsh_render_attrs['dielectric']
            dielectric_gaps = gmsh_render_attrs['dielectric_gaps']
            PEC_metals = gmsh_render_attrs['metals']
            far_field = gmsh_render_attrs['far_field']
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
                if cur_port['elem_type'] == 'cpw':
                    ports[cur_port['port_name'] + 'a'] = list(comsol._model.java.component("comp1").physics("emw").feature(f"lport{m}").feature('ue2').selection().entities())[0]
                    ports[cur_port['port_name'] + 'b'] = list(comsol._model.java.component("comp1").physics("emw").feature(f"lport{m}").feature('ue1').selection().entities())[0]
                else:
                    ports[cur_port['port_name']] = list(comsol._model.java.component("comp1").physics("emw").feature(f"lport{m}").selection().entities())[0]

            #define length scale
            l0 = 1

            #file extension
            file_ext = '.mphbin'

        #get material parameters
        dielectric = Material(self.user_options["dielectric_material"])

        #Process Ports
        config_ports = self._process_ports(ports)
        if self._rf_port_excitation > 0:
            config_ports[self._rf_port_excitation-1]["Excitation"] = True

        #Define python dictionary to convert to json file
        if self._output_subdir == "":
            self.set_local_output_subdir("", False)
        filePrefix = self.hpc_options["input_dir"]  + self.name + "/" if self.hpc_options["input_dir"] != "" else ""
        self._mesh_name = filePrefix + self.name + file_ext
        config = {
            "Problem":
            {
                "Type": "Eigenmode",
                "Verbose": 2,
                "Output": filePrefix  + "outputFiles"
            },
            "Model":
            {
                "Mesh":  self._mesh_name,
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
                        "LossTan": dielectric.loss_tangent
                    }
                ]
            },
            "Boundaries":
            {
                "PEC":
                {
                    "Attributes": PEC_metals,  # Metal trace
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
        if self.meshing == 'GMSH':
            self._setup_EPR_boundaries(config, dielectric_gaps, PEC_metals)
        if self._ff_type == 'absorbing':
            config['Boundaries']['Absorbing'] = {
                    "Attributes": far_field,
                    "Order": 1
                }
        else:
            config['Boundaries']['PEC']['Attributes'] += far_field
        # if self.meshing == 'GMSH':
        #     config['Solver']['Linear']['Type'] = "Default"
        #     config['Solver']['Linear']['KSPType'] = "GMRES"

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
    
    def set_freq_search(self, min_freq_Hz, num_freq):
        self.user_options["starting_freq"] = min_freq_Hz / 1e9
        self.user_options["number_of_freqs"] = num_freq
        if self._sim_config != "":
            with open(self._sim_config, "r") as f:
                config_json = json.loads(f.read())
            config_json['Solver']['Eigenmode']['Target'] = min_freq_Hz / 1e9
            config_json['Solver']['Eigenmode']['N'] = num_freq
            with open(self._sim_config, "w") as f:
                json.dump(config_json, f, indent=2)

    def retrieve_data(self):
        raw_data = pd.read_csv(self._output_data_dir + '/eig.csv')
        headers = raw_data.columns
        raw_data = raw_data.to_numpy()

        lePlots = self._output_data_dir + '/paraview/eigenmode/eigenmode.pvd'
        if os.path.exists(lePlots):
            leView = PVDVTU_Viewer(lePlots)
            for m in range(leView.num_datasets):
                try: 
                    leSlice = leView.get_data_slice(m)
                    fig = leSlice.plot(np.linalg.norm(leSlice.get_data('E_real'), axis=1), 'coolwarm', True)
                    fig.savefig(self._output_data_dir + f'/eig{m}_ErealMag.png')
                    plt.close(fig)
                except Exception as e:
                    print(f"Error in plotting: {e}")
            try:
                fig = leSlice.plot_mesh()
                fig.savefig(self._output_data_dir + '/mesh.png')
                plt.close(fig)
            except Exception as e:
                print(f"Error in plotting: {e}")

    def retrieve_field_plots(self):
        lePlots = self._output_data_dir + '/paraview/eigenmode/eigenmode.pvd'
        return PVDVTU_Viewer(lePlots)

    def retrieve_EPR_data(self):
        return PALACE_Eigenmode_Simulation.retrieve_EPR_data_from_file(self._sim_config, self._output_data_dir)

    @staticmethod
    def retrieve_EPR_data_from_file(json_sim_config, output_directory):
        if os.path.exists(output_directory + '/config.json'):
            #Newer version stores configuration in output directory... Use this as it is the exact one used in this simulation...
            json_sim_config = output_directory + '/config.json'
        with open(json_sim_config, "r") as f:
            config_json = json.loads(f.read())
        assert 'Boundaries' in config_json, "\"Boundaries\" not found in the JSON configuration."
        leDict = config_json['Boundaries']
        assert 'Postprocessing' in leDict, "\"Postprocessing\" not found in the JSON configuration."
        leDict = leDict['Postprocessing']
        assert 'Dielectric' in leDict, "\"Dielectric\" not found in the JSON configuration."
        leDict = leDict['Dielectric']
        epr_dict = {}
        for cur_dict in leDict:
            epr_dict[cur_dict['Type']] = cur_dict['Index']

        raw_data = pd.read_csv(output_directory + '/eig.csv')
        headers = raw_data.columns
        raw_data = raw_data.to_numpy()        
        col_Ref = [x for x in range(len(headers)) if headers[x].strip().startswith(r'Re{f}')][0]
        col_Ief = [x for x in range(len(headers)) if headers[x].strip().startswith(r'Im{f}')][0]
        col_Q = [x for x in range(len(headers)) if headers[x].strip().startswith(r'Q')][0]

        raw_dataEPR = pd.read_csv(output_directory + '/surface-Q.csv')
        headersEPR = raw_dataEPR.columns
        raw_dataEPR = raw_dataEPR.to_numpy()

        ret_list = []
        for m,cur_row in enumerate(raw_data):
            cur_mode = {}
            cur_mode['Frequency'] = (cur_row[col_Ref] + 1j*cur_row[col_Ief]) * 1e9
            cur_mode['Q'] = cur_row[col_Q]
            for cur_interface in epr_dict:
                cur_mode[cur_interface] = {}
                cur_mode[cur_interface]['p'] = raw_dataEPR[m,2*epr_dict[cur_interface]-1]
                cur_mode[cur_interface]['Q'] = raw_dataEPR[m,2*epr_dict[cur_interface]]
            ret_list.append(cur_mode)
        
        return ret_list
