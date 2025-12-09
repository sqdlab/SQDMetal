# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0

from SQDMetal.PALACE.Model import PALACE_Model_RF_Base
from SQDMetal.Utilities.Materials import Material
from SQDMetal.PALACE.Utilities.GMSH_Navigator import GMSH_Navigator
from SQDMetal.PALACE.PVDVTU_Viewer import PVDVTU_Viewer
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
import numpy as np
import json
import os
import re
import io
import sys
import pandas as pd

class PALACE_Eigenmode_Simulation(PALACE_Model_RF_Base):

    #Class Variables
    default_user_options = {
                 "starting_freq": 5.5e9,
                 "number_of_freqs": 5,
                 "solns_to_save": 5,
                 "solver_order": 2,
                 "solver_tol": 1.0e-8,
                 "solver_maxits": 100,
                 "HPC_Parameters_JSON": ""
                }

    #Parent Directory path
    simPC_parent_simulation_dir = "/home/experiment/PALACE/Simulations/input"

    #constructor
    def __init__(self, name, sim_parent_directory, mode, meshing, user_options = {}, 
                 view_design_gmsh_gui = False, create_files = False, **kwargs):
        self.name = name
        self.sim_parent_directory = sim_parent_directory
        self.mode = mode
        self.user_options = {}
        self._ff_type = {}
        self.set_farfield()
        for key in PALACE_Eigenmode_Simulation.default_user_options:
            self.user_options[key] = user_options.get(key, PALACE_Eigenmode_Simulation.default_user_options[key])
        self.view_design_gmsh_gui = view_design_gmsh_gui
        self.create_files = create_files
        self._ports = []
        super().__init__(meshing, mode, user_options, **kwargs)

    def _create_config_file(self, **kwargs):
        '''Create the configuration file for the simulation.
        
           This function creates the .json file which contains all the details of the simulation including simulation type,
           boundary conditions, postprocessing and solver conditions.

           Args:
                **kwargs: Arbitrary keyword arguments.
            
           Returns:
                None
        
        '''    

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
            #TODO: Add COMSOL Meshing support for separate far-field BCs...
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
                "LumpedPort": config_ports,
                "WavePort": config_wports
            },
            "Solver":
            {
                "Order": self.user_options["solver_order"],
                "Eigenmode":
                {
                    "N": self.user_options["number_of_freqs"],  # number of eigenfrequencies
                    "Tol": self.user_options["solver_tol"],  # solver tolerance
                    "Target": self.user_options["starting_freq"]/1e9,  # GHz - starting point
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
        
        #If kinetic inductance is incorporated change metals from PEC to Impedance boundary condition. This is done before  
        if self._use_KI:
            self._setup_kinetic_inductance(config, PEC_metals)

        if self.meshing == 'GMSH':
            self._setup_EPR_boundaries(config, dielectric_gaps, PEC_metals)
        
        self._process_farfield(config, far_field)
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
        '''Retrieve data from simulation output files.
        
           This function retrieves data from the files produced from the eigenmode simulation and creates plots of the electric fields.

           Args:
                None.
            
           Returns:
                None.
        
        '''   

        raw_data = pd.read_csv(self._output_data_dir + '/eig.csv')
        headers = raw_data.columns # noqa: F841 # abhishekchak52: headers is not used
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
            if self.meshing == 'GMSH':
                GMSH_Navigator(self._mesh_path).export_to_png(file_path=self._output_data_dir + '/mesh.png')
            else:
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

    def retrieve_mode_port_EPR(self):
        return PALACE_Eigenmode_Simulation.retrieve_mode_port_EPR_from_file(self._output_data_dir)

    def calculate_hamiltonian_parameters_EPR(self, modes_to_compare = [], print_output=True):
        return PALACE_Eigenmode_Simulation.calculate_hamiltonian_parameters_EPR_from_files(self._output_data_dir, self._sim_config, modes_to_compare=modes_to_compare, print_output=print_output)

    def plot_fields_with_data(self, save=True, columns=4):
        return PALACE_Eigenmode_Simulation.plot_fields_with_data_from_files(self._output_data_dir, save=save, columns=columns)

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
        headersEPR = raw_dataEPR.columns # noqa: F841 # abhishekchak52: headersEPR is not used
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

    @staticmethod
    def retrieve_mode_port_EPR_from_file(output_directory):
        #Returns matrix where modes on rows, ports on columns

        raw_data = pd.read_csv(output_directory + '/port-EPR.csv')
        headers = raw_data.columns
        raw_data = raw_data.to_numpy()

        raw_dataE = pd.read_csv(output_directory + '/eig.csv')
        headers = raw_dataE.columns
        raw_dataE = raw_dataE.to_numpy()
        col_Ref = [x for x in range(len(headers)) if headers[x].strip().startswith(r'Re{f}')][0]
        col_Ref_Q = [x for x in range(len(headers)) if headers[x].strip().startswith(r'Q')][0]

        return {'mat_mode_port': raw_data[:, 1:], 'eigenfrequencies': raw_dataE[:, col_Ref]*1e9, 'loaded_Q': raw_dataE[:, col_Ref_Q]}

    @staticmethod
    def calculate_hamiltonian_parameters_EPR_from_files(directory, config_json_path, modes_to_compare = [], print_output=True):
        '''Extracts Hamiltonian parameters from an eigenmode simulation using the EPR method'''
        
        #retrieve data from the simulation files
        mode_dict = PALACE_Eigenmode_Simulation.retrieve_mode_port_EPR_from_file(directory)
        
        #Constants used for calculations
        e = 1.60218e-19
        hbar = 1.05457e-34
        phi0 = hbar/(2*e)

        #get particpation ratios and frequencies
        participation_ratios = mode_dict['mat_mode_port']
        frequencies = mode_dict['eigenfrequencies']

        #get #modes and #junctions from arrays
        num_of_modes = np.shape(participation_ratios)[0]
        num_of_junctions = np.shape(participation_ratios)[1]

        #choose modes user is interested in
        modes = np.arange(num_of_modes) + 1 # array of modes
        modes_to_compare = np.array(modes_to_compare) #cast to array
        if modes_to_compare.size == 0:
            modes_to_compare = modes
        modes_to_delete = np.setdiff1d(modes, modes_to_compare)

        #update number of modes left for comparison
        num_of_modes = num_of_modes - modes_to_delete.size

        #create labels for the modes
        dataframe_label_mode = []
        for i, mode in enumerate(modes_to_compare):
            dataframe_label_mode.append("mode_" + str(mode))

        dataframe_label_junc = []
        for i in range(num_of_junctions):
            dataframe_label_junc.append("jj_" + str(i+1))

        #delete modes not wanted from participation ratio and frequency vectors
        if modes_to_delete.size > 0:
            for i,mode in enumerate(modes_to_delete):
                participation_ratios = np.delete(participation_ratios, (modes_to_delete[i] - (i+1)), axis=0)
                frequencies = np.delete(frequencies, (modes_to_delete[i] - (i+1)), axis=0)

        #participation ratio matrix - dimension (num_of_modes, num_of_junctions)
        P = np.matrix(np.abs(participation_ratios))

        #diagonal matrix for eigenfrequencies
        omega_vector = 2* np.pi * frequencies
        omega = np.diag(omega_vector)

        #Detuning matrix
        delta = np.ones((np.shape(omega)[0], np.shape(omega)[0])) #create delta/detuning matrix

        #Populate delta/detuning matrix
        for j in range(np.shape(omega)[0]):  
            delta[:,j] = omega_vector - omega_vector[j]

        #Extract port-inductances
        with open(config_json_path, "r") as f:
            config_json = json.loads(f.read())
        Lj = []
        for cur_port in config_json['Boundaries']['LumpedPort']:
            if 'L' in cur_port and cur_port['L'] > 0:
                Lj.append(cur_port['L'])
            #Note that Palace only stores the port-EPR values for inductive ports.
        #Get Ej for each junction
        Lj_matrix = np.matrix(Lj)
        Ej = np.diag((phi0**2 / Lj_matrix).A1)

        #Calcultate Kerr matrix: chi = hbar/4 * (omega*P) * inv_Ej * (omega*P)^T
        ###Please note all these values are angular frequencies radians/sec, to get frequencies in Hertz divide by 2*pi###
        chi = (hbar/4) * (omega*P) * np.linalg.inv(Ej) * np.transpose(omega*P)
        diag_matrix = (-1/2) * np.diag(np.ones(num_of_modes)) + np.ones(num_of_modes) #hacky way of creating ones matrix with values of 1/2 down the diagonal
        chi_anharm = np.multiply(diag_matrix,chi) #This is the kerr matrix but with the anharmonicity down the diagonal for the self-kerr term
        anharm = (1/2) * np.diagonal(chi)
        lamb_shift = (1/2) * chi * np.ones((np.shape(chi)[0],1))
        renormalised_freqs = np.transpose(np.matrix(frequencies * 2 * np.pi)) - lamb_shift #linear frequencies from eigenmode simulation are renormalised by lamb shift caused by non-linearity

        #Get values in megahertz (MHz) or gigahertz (GHz) to display to user
        chi_freq = chi / (2 * np.pi * 1e6)  # noqa: F841 # abhishekchak52: chi_freq is not used
        chi_anharm = chi_anharm / (2 * np.pi * 1e6) 
        anharm_freq = anharm / (2 * np.pi * 1e6)  # noqa: F841 # abhishekchak52: anharm_freq is not used
        lamb_shift_freq = lamb_shift / (2 * np.pi * 1e6)  
        delta_freq = delta / (2 * np.pi * 1e9) 
        renormalised_freqs = renormalised_freqs /  (2 * np.pi * 1e9)

        #wrapper to print to both notebook and buffer
        output_buffer = io.StringIO()
        original_stdout = sys.stdout
        def dual_print(*args, **kwargs):
            """Prints to both notebook and the output buffer."""
            print(*args, **kwargs)             # normal notebook output
            print(*args, **kwargs, file=output_buffer)  # captured copy

        freq_df = pd.DataFrame(frequencies / 1e9, dataframe_label_mode, ['Freq (GHz)'])
        freq_renorm_df = pd.DataFrame(renormalised_freqs, dataframe_label_mode, ['Freq (GHz)'])
        ratios_df = pd.DataFrame(np.abs(participation_ratios), dataframe_label_mode, dataframe_label_junc)
        chi_anharm_df = pd.DataFrame(chi_anharm, dataframe_label_mode, dataframe_label_mode)
        lamb_shift_df = pd.DataFrame(lamb_shift_freq, dataframe_label_mode, ['Lamb Shift (MHz)'])
        delta_df = pd.DataFrame(delta_freq, dataframe_label_mode, dataframe_label_mode)
        if print_output:
            dual_print('Simulation Mode Frequencies:')
            dual_print(freq_df)
            dual_print('______________________________\n')
            dual_print('Renormalised Frequencies:\n ***Freqs from simulation are adjusted for Lamb shift')
            dual_print(freq_renorm_df)
            dual_print('______________________________\n')
            dual_print('Participation Ratios:')
            dual_print(ratios_df)
            dual_print('______________________________\n')
            dual_print('Chi Matrix (MHz):\n ***Diag is Anharmonicity, Off-Diag is Cross-Kerr')
            dual_print(chi_anharm_df)
            dual_print('______________________________\n')
            dual_print('Lamb Shifts:')
            dual_print(lamb_shift_df)
            dual_print('______________________________\n')
            dual_print('Detuning (GHz):')
            dual_print(delta_df)
            dual_print('______________________________\n')

        #sSave output to file
        output_text = output_buffer.getvalue()
        file_path = directory + '/EPR_params.txt'
        with open(file_path, 'w') as f:
            f.write(output_text)

        #TODO: Strongly consider refactoring these keys/values...
        return {'f_modes_GHz': freq_df, 'f_norms_GHz': freq_renorm_df, 'EPR': ratios_df, 'Chi': chi_anharm_df, 'Lamb': lamb_shift_df, 'Detuning': delta_df}
    
    @staticmethod
    def plot_fields_with_data_from_files(directory, save=True, columns=4, skip_postprocessing=False):
        mode_dict = PALACE_Eigenmode_Simulation.retrieve_mode_port_EPR_from_file(
            directory)
        if not skip_postprocessing:
            try:
                participations = PALACE_Eigenmode_Simulation.retrieve_EPR_data_from_file(
                    json_sim_config=directory+r"/config.json", output_directory=directory)
            except:
                skip_postprocessing = True
        # collect png files
        r = re.compile(r"eig(\d+)_ErealMag\.png")
        png_files = []
        for filename in sorted(os.listdir(directory)):
            if r.match(filename):
                png_files.append(filename)
        n_modes = len(png_files)
        # plotting
        fig = plt.figure(figsize=(8, 3 * n_modes))
        gs = GridSpec(n_modes, 2, width_ratios=[2, 1], figure=fig)
        for i, filename in enumerate(png_files):
            eig_num = int(r.match(filename).group(1))
            img = mpimg.imread(os.path.join(directory, filename))
            # add image
            ax_img = fig.add_subplot(gs[i, 0])
            ax_img.imshow(img)
            ax_img.axis("off")
            # add text
            ax_txt = fig.add_subplot(gs[i, 1])
            ax_txt.axis("off")
            if skip_postprocessing:
                ax_txt.text(
                    0, 0.8,
                    f"Mode {eig_num}\n\n"
                    f"f = {mode_dict['eigenfrequencies'][eig_num].real * 1e-9:.3f} GHz\n"
                    f"Q = {mode_dict['loaded_Q'][eig_num]}"
                )
            else:
                ax_txt.text(
                    0, 0.8,
                    f"Mode {eig_num}\n\n"
                    f"f = {mode_dict['eigenfrequencies'][eig_num].real * 1e-9:.3f} GHz\n"
                    f"Q = {participations[eig_num]['Q']:.0f}\n"
                    f"kappa = {2*np.pi*mode_dict['eigenfrequencies'][eig_num].real/participations[eig_num]['Q'] * 1e-6:.3f} MHz\n"
                    f"p_MS = {participations[eig_num]['MS']['p']:.2e}\n"
                    f"p_MA = {participations[eig_num]['MA']['p']:.2e}\n"
                    f"p_SA = {participations[eig_num]['SA']['p']:.2e}\n",
                    fontsize=10,
                    va="top"
                )
        plt.tight_layout()
        if save:
            plt.savefig(directory + r"/eigenmodes.png")
            print(f"Plot saved: {directory}/eigenmodes.png")
        plt.show()
        if skip_postprocessing:
            return {'mode_dict': mode_dict}
        else:
            return {'mode_dict': mode_dict, 'participations': participations}

