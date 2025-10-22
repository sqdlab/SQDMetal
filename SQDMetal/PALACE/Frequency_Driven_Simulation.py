# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0

from SQDMetal.PALACE.Model import PALACE_Model_RF_Base
from SQDMetal.Utilities.Materials import Material
import matplotlib.pyplot as plt
import numpy as np
import json
import os
import pandas as pd
import re

class PALACE_Driven_Simulation(PALACE_Model_RF_Base):

    #Class Variables
    default_user_options = {
                 "fillet_resolution": 4,
                 "dielectric_material": "silicon",
                 "solns_to_save": 4,
                 "solver_order": 2,
                 "solver_tol": 1.0e-8,
                 "solver_maxits": 300,
                 "mesh_max": 100e-6,
                 "mesh_min": 10e-6,
                 "taper_dist_min": 30e-6,
                 "taper_dist_max": 200e-6,
                 "gmsh_dist_func_discretisation": 120,
                 "comsol_meshing": "Extremely fine",
                 "HPC_Parameters_JSON": "",
                 "fuse_threshold": 1e-9,
                 "gmsh_verbosity": 1,
                 "threshold": 1e-9,
                 "simplify_edge_min_angle_deg": -1,
                 'palace_mode': 'local',
                 'palace_wsl_spack_repo_directory': '~/repo'
                }

    #constructor
    def __init__(self, name, sim_parent_directory, mode, meshing, user_options = {}, 
                 view_design_gmsh_gui = False, create_files = False, **kwargs):
        self.name = name
        self.sim_parent_directory = sim_parent_directory
        self.mode = mode
        self.meshing = meshing
        self.user_options = {}
        for key in PALACE_Driven_Simulation.default_user_options:
            self.user_options[key] = user_options.get(key, PALACE_Driven_Simulation.default_user_options[key])
        self.view_design_gmsh_gui = view_design_gmsh_gui
        self.create_files = create_files
        self._ports = []
        super().__init__(meshing, mode, user_options, **kwargs)
        self.freqs = None
        self._current_sources = []


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
        config_ports, config_wports = self._process_ports(ports)
        if self._rf_port_excitation > 0:
            for cur_port in config_ports:
                if cur_port['Index'] == self._rf_port_excitation:
                    cur_port['Excitation'] = True
                    break
            for cur_port in config_wports:
                if cur_port['Index'] == self._rf_port_excitation:
                    cur_port['Excitation'] = True
                    break

        if isinstance(self.freqs, tuple):
            fStart,fStop,fStep = self.freqs[0]/1e9, self.freqs[1]/1e9, self.freqs[2]/1e9
        else:
            print("Warning: Frequencies have not been set properly...")
            fStart,fStop,fStep = (1,10,1)

        #Define python dictionary to convert to json file
        if self._output_subdir == "":
            self.set_local_output_subdir("", False)
        filePrefix = self._get_folder_prefix()
        self._mesh_name = filePrefix + self.name + file_ext
        config = {
            "Problem":
            {
                "Type": "Driven",
                "Verbose": 2,
                "Output": self._output_dir
            },
            "Model":
            {
                "Mesh":  self._mesh_name,
                "L0": l0,  
                "CrackDisplacementFactor":0,    #TODO: Remove if it is not required for both planar AND full-3D designs with CPW feeds. c.f. https://awslabs.github.io/palace/dev/config/model/
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
                "LumpedPort": config_ports,
                "WavePort": config_wports
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
        if self.meshing == 'GMSH':
            self._setup_EPR_boundaries(config, dielectric_gaps, PEC_metals)
        self._process_farfield(config, far_field)
        #Add in any surface current sources
        if len(self._current_sources) > 0:
            config['Boundaries']['SurfaceCurrent'] = []
            index = 1   #In case, the current source is rejected, keep index counting those actually added...
            for cur_cur_src in self._current_sources:
                if cur_cur_src[0] == 'region':
                    new_src = {'Index': index, "Direction": [0.0, 0.0, cur_cur_src[2]]}
                    if cur_cur_src[1] == 'metals':
                        new_src['Attributes'] = PEC_metals
                    elif cur_cur_src[1] == 'dielectric':
                        new_src['Attributes'] = dielectric_gaps
                    config['Boundaries']['SurfaceCurrent'].append(new_src)
                    index += 1

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
    
    def add_surface_current_source_region(self, region="", direction=1):
        assert region in ['metals', 'dielectric'], "Must set region to either 'metals' or 'dielectric'"
        assert direction == -1 or direction == 1, "Direction must be either 1 or -1 for +z or -z"
        self._current_sources.append(('region', region, direction))

    def set_freq_values(self, freq_start, freq_end, freq_step):
        #Frequencies in Hertz
        self.freqs = (freq_start, freq_end, freq_step)
        if self._sim_config != "":
            with open(self._sim_config, "r") as f:
                config_json = json.loads(f.read())
            config_json['Solver']['Driven']['MinFreq'] = freq_start / 1e9
            config_json['Solver']['Driven']['MaxFreq'] = freq_end / 1e9
            config_json['Solver']['Driven']['FreqStep'] = freq_step / 1e9
            with open(self._sim_config, "w") as f:
                json.dump(config_json, f, indent=2)

    @staticmethod
    def retrieve_data_from_file(file_path_port_S_csv):
        raw_data = pd.read_csv(file_path_port_S_csv)
        headers = raw_data.columns
        raw_data = raw_data.to_numpy().T

        num_s_params = int((raw_data.shape[0]-1)/2)

        freq_vals = raw_data[0]*1e9

        ret_data = {'freqs' : freq_vals}
        for m in range(num_s_params):
            #Header will be like: |S[1][1]| (dB)
            key = headers[2*m+1].strip().split('|')[1].replace('[','').replace(']','')
            amp = np.array(raw_data[2*m+1], dtype=np.float64)
            phs = np.array(raw_data[2*m+2], dtype=np.float64)
            sParam = 10**(amp/20)*np.exp(1j*phs/180*np.pi)
            ret_data[key] = sParam

        return ret_data

    def retrieve_data(self):
        if not os.path.exists(self._output_data_dir + '/port-S.csv'):
            return None

        raw_data = pd.read_csv(self._output_data_dir + '/port-S.csv')
        headers = raw_data.columns
        raw_data = raw_data.to_numpy().T

        num_s_params = int((raw_data.shape[0]-1)/2)

        fig, axs = plt.subplots(ncols=num_s_params)
        fig.set_figwidth(8*num_s_params)

        freq_vals = raw_data[0]*1e9

        ret_data = {'freqs' : freq_vals}
        for m in range(num_s_params):
            #Header will be like: |S[1][1]| (dB)
            key = headers[2*m+1].strip().split('|')[1].replace('[','').replace(']','')
            amp = np.array(raw_data[2*m+1], dtype=np.float64)
            phs = np.array(raw_data[2*m+2], dtype=np.float64)
            sParam = 10**(amp/20)*np.exp(1j*phs/180*np.pi)
            ret_data[key] = sParam
            #
            axs[m].grid()
            axs[m].set_title(f'|{key}| (dB)')
            axs[m].plot(freq_vals, amp, 'r-')
            axs[m].set_xlabel('Frequency (Hz)')
            axs[m].set_ylabel('Amplitude (dB)', color='red')
            axsPhs = axs[m].twinx()
            axsPhs.plot(freq_vals, phs, 'b-')
            axsPhs.set_ylabel('Phase (°)', color='blue')
            for label in axs[m].yaxis.get_ticklabels():
                label.set_color('red')
            for label in axsPhs.yaxis.get_ticklabels():
                label.set_color('blue')

        fig.tight_layout()
        fig.savefig(self._output_data_dir + '/s_params.png')
        plt.close(fig)

        return ret_data

    def get_waveport_modes(self):
        return PALACE_Driven_Simulation.get_waveport_modes_from_file(self.log_location, self._output_data_dir + '/port-S.csv')

    @staticmethod
    def get_waveport_modes_from_file(filename_outLog, filename_PortScsv):
        '''
        Returns the waveport modes by reading the out.log file and its associated port-S.csv file. 
        '''
        #The string is usually like:
        #   Port 1, mode 1: kₙ = 1.737e+02-5.076e-04i m⁻¹
        lines_with_k = []
        pattern = re.compile("Port (.*), mode (.*):")

        with open(filename_outLog, "r") as f:
            for line in f:
                cur_line = line.rstrip()
                cur_pattern = pattern.search(cur_line)
                if cur_pattern:
                    port_num = cur_pattern.group(1)
                    mode_num = cur_pattern.group(2)
                    cur_line = cur_line.split('=')[1].split('i')[0]
                    lines_with_k.append([int(port_num), int(mode_num), np.complex128(cur_line+'j')])
        with open(filename_PortScsv, "r") as f:
            raw_data = pd.read_csv(filename_PortScsv)
            raw_data = raw_data.to_numpy().T
            freq_vals = raw_data[0]*1e9
        num_data_per_freq = int(len(lines_with_k) / freq_vals.size)
        assert num_data_per_freq*freq_vals.size == len(lines_with_k), "For some reason the number of ports/modes found does not match the number of simulated frequencies. Check that the simulation completed properly."
        
        ret_data = []
        for m,cur_freq in enumerate(freq_vals):
            for n in range(num_data_per_freq):
                ret_data.append([cur_freq] + lines_with_k.pop(0))

        return ret_data
