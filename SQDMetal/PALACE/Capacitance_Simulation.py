# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0

from SQDMetal.PALACE.Model import PALACE_Model_Base
from SQDMetal.COMSOL.Model import COMSOL_Model
from SQDMetal.COMSOL.SimCapacitance import COMSOL_Simulation_CapMats
from SQDMetal.Utilities.Materials import Material
from SQDMetal.Utilities.QubitDesigner import FloatingTransmonDesigner
import matplotlib.pyplot as plt
import json
import os
import gmsh
import shapely
import numpy as np
import pandas as pd
import geopandas as gpd
from SQDMetal.PALACE.Utilities.GMSH_Geometry_Builder import GMSH_Geometry_Builder
from SQDMetal.PALACE.Utilities.GMSH_Mesh_Builder import GMSH_Mesh_Builder
from SQDMetal.PALACE.PVDVTU_Viewer import PVDVTU_Viewer

class PALACE_Capacitance_Simulation(PALACE_Model_Base):

    #Class Variables
    default_user_options = {
                 "solns_to_save": -1,
                 "solver_order": 2,
                 "solver_tol": 1.0e-8,
                 "solver_maxits": 100,
                 "HPC_Parameters_JSON": ""
    }

    # Parent Directory path
    HPC_parent_simulation_dir = "C:/PALACE_Simulations/"
    simPC_parent_simulation_dir = "/home/experiment/PALACE/Simulations/input"

    #constructor
    def __init__(self, name, sim_parent_directory, mode, meshing, user_options = {}, 
                                    view_design_gmsh_gui = False, create_files = False, **kwargs):
        self.name = name
        self.sim_parent_directory = sim_parent_directory
        self.mode = mode
        self.meshing = meshing
        self.user_options = {}
        for key in PALACE_Capacitance_Simulation.default_user_options:
            self.user_options[key] = user_options.get(key, PALACE_Capacitance_Simulation.default_user_options[key])
        self.view_design_gmsh_gui = view_design_gmsh_gui
        self.create_files = create_files
        self._cur_cap_terminals = []
        self.cap_matrix = None
        super().__init__(meshing, mode, user_options, **kwargs)

    def _create_directory(self, directory_name):
        '''create a directory to hold the simulation files'''

        if self.create_files:
            parent_simulation_dir = self._check_simulation_mode()

            # Directory
            directory = directory_name
    
            # Path
            path = os.path.join(parent_simulation_dir, directory)
    
            # Create the directory
            if not os.path.exists(path):
                os.mkdir(path)
                # print("Directory '% s' created" % directory)

    def _save_mesh_gmsh(self):
        parent_simulation_dir = self._check_simulation_mode()

        # file_name
        file_name = self.name + "/" + self.name + ".msh"
  
        # Path
        path = os.path.join(parent_simulation_dir, file_name)
        gmsh.write(path)
        
    def _save_mesh_comsol(self, comsol_obj):
        '''function used to save the comsol mesh file'''
        if self.create_files:
            parent_simulation_dir = self._check_simulation_mode()

            # file_name
            file_name = self.name + "/" + self.name + ".mphbin"
    
            # Path
            path = os.path.join(parent_simulation_dir, file_name)

            #COMSOL export commands
            comsol_obj._model.java.component("comp1").mesh("mesh1").export().set("filename", path)
            comsol_obj._model.java.component("comp1").mesh("mesh1").export(path)

    def _create_config_file(self, **kwargs):
        '''create the configuration file which specifies the simulation type and the parameters'''    

        if self.meshing == 'GMSH':

            #GMSH renderer object needed to get boundary conditions for config file
            gmsh_render_attrs = kwargs['gmsh_render_attrs']

            #GMSH config file variables
            material_air = gmsh_render_attrs['air_box']
            material_dielectric = gmsh_render_attrs['dielectric']
            far_field = [gmsh_render_attrs['far_field'][x] for x in gmsh_render_attrs['far_field']]
            
            self._cur_cap_terminals = gmsh_render_attrs['metalsShapely']

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

            self._cur_cap_terminals = [x[0] for x in sCapMats.model._cond_polys]

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
            self.set_local_output_subdir("")
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
                        "Save": len(Terminal) if self.user_options["solns_to_save"] == -1 else self.user_options["solns_to_save"]
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

    def prepare_simulation(self):
        '''Creates and saves the GMSH (.msh) and Palace configuration file (.json). 
        '''

        if self.meshing == 'GMSH':
            ggb = GMSH_Geometry_Builder(self._geom_processor, self.user_options['fillet_resolution'], self.user_options['gmsh_verbosity'])
            gmsh_render_attrs = ggb.construct_geometry_in_GMSH(self._metallic_layers, self._ground_plane, [],
                                                               self._fine_meshes, self.user_options["fuse_threshold"],
                                                               threshold=self.user_options["threshold"],
                                                               simplify_edge_min_angle_deg=self.user_options["simplify_edge_min_angle_deg"],
                                                               full_3D_params = self._full_3D_params,
                                                               boundary_distances = self._boundary_distances)
            #
            gmb = GMSH_Mesh_Builder(gmsh_render_attrs['fine_mesh_elems'], self.user_options)
            gmb.build_mesh()

            if self.create_files:
                #create directory to store simulation files
                self._create_directory(self.name)

                #create config file
                self._create_config_file(gmsh_render_attrs = gmsh_render_attrs)

                #create batch file
                # if self.mode == 'HPC':
                #     self.create_batch_file()

                self._save_mesh_gmsh()

            # abhishekchak52: commented out for now since pgr is not defined
            # if self.view_design_gmsh_gui:
            #     #plot design in gmsh gui
            #     pgr.view_design_components()
            
        if self.meshing == 'COMSOL':

            #initialise the COMSOL engine
            COMSOL_Model.init_engine()
            cmsl = COMSOL_Model('res')

            #Create COMSOL capacitance sim object
            simCapMats = COMSOL_Simulation_CapMats(cmsl)
            cmsl.initialize_model(self._geom_processor.design, [simCapMats], bottom_grounded = True, resolution = 10)     #TODO: Make COMSOL actually compatible rather than assuming Qiskit-Metal designs?

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

            if self.create_files:
                #create directory to store simulation files
                self._create_directory(self.name)

                #create config file
                self._create_config_file(comsol_obj = cmsl, simCap_object = simCapMats)

                #create batch file
                if self.mode == 'HPC':
                    self.create_batch_file()

            #save mesh
            self._save_mesh_comsol(comsol_obj = cmsl)

    def display_conductor_indices(self, save=False):
        '''Plots a coloured visualisation of the metallic conductors and their corresponding row/column indices of the capacitance matrix.
        
        Args:
            save (bool, optional): (Defaults to False) Choose whether or not to save the generated plot. 
        
        Returns:
            Matplotlib fig of the conductor visualisation.
        '''
        assert len(self._cur_cap_terminals) > 0, "There are no terminals. Ensure prepare_simulation() has been called."

        minX = (self._geom_processor.chip_centre[0] - self._geom_processor.chip_size_x*0.5)
        maxX = (self._geom_processor.chip_centre[0] + self._geom_processor.chip_size_x*0.5)
        minY = (self._geom_processor.chip_centre[1] - self._geom_processor.chip_size_y*0.5)
        maxY = (self._geom_processor.chip_centre[1] + self._geom_processor.chip_size_y*0.5)
        chip_bounding_poly = shapely.LineString([(minX,minY),(maxX,minY),(maxX,maxY),(minX,maxY),(minX,minY)])
        
        leGeoms = [chip_bounding_poly] + self._cur_cap_terminals
        leNames = ["Chip"] + [f"Cond{x+1}" for x in range(len(self._cur_cap_terminals))]
        gdf = gpd.GeoDataFrame({'names':leNames}, geometry=leGeoms)
        fig, ax = plt.subplots(1)
        gdf.plot(ax = ax, column='names', cmap='jet', alpha=0.5, categorical=True, legend=True)
        ax.set_xlabel('Position (mm)')
        ax.set_ylabel('Position (mm)')

        if save==True:
            fig.savefig("ConductorIndicides.png")

        return fig

    def retrieve_data(self):
        '''
        Retrieves output data from the capacitance simulation.
        
        Creates and saves plots of conductor indices (terminal_indices.png), field distribution results (cond1_V.png), and mesh (mesh.png)

        This function must be run after calling :func:`~SQDMetal.PALACE.Model.PALACE_Model_Base.run`.

        Returns:
            A NumPy array of the raw data.

        '''
        raw_data = pd.read_csv(self._output_data_dir + '/terminal-C.csv')
        self.cap_matrix = raw_data
        headers = raw_data.columns # noqa: F841 # abhishekchak52: headers is not used
        raw_data = raw_data.to_numpy()[:,1:]    #First column is just the indices...

        fig = self.display_conductor_indices()
        fig.savefig(self._output_data_dir + '/terminal_indices.png')
        plt.close(fig)

        lePlots = self._output_data_dir + '/paraview/electrostatic/electrostatic.pvd'
        if os.path.exists(lePlots):
            leView = PVDVTU_Viewer(lePlots)
            for m in range(leView.num_datasets):
                try: 
                    leSlice = leView.get_data_slice(m)
                    fig = leSlice.plot(leSlice.get_data('V'), 'coolwarm', True)
                    fig.savefig(self._output_data_dir + f'/cond{m+1}_V.png')
                    plt.close(fig)
                except Exception as e:
                    print(f"Error in plotting: {e}")
            try:
                fig = leSlice.plot_mesh()
                fig.savefig(self._output_data_dir + '/mesh.png')
                plt.close(fig)
            except Exception as e:
                print(f"Error in plotting: {e}")

        return raw_data
    
    def calc_params_floating_Transmon(self, **kwargs):
        """
        Class method mirroring :func:`~SQDMetal.PALACE.Capacitance_Simulation.PALACE_Capacitance_Simulation.calc_params_floating_Transmon_from_files` (without the ``directory`` parameter as it is inferred from the current simulation)

        This function must be run after calling :func:`~SQDMetal.PALACE.Model.PALACE_Model_Base.run`.
        """
        return PALACE_Capacitance_Simulation.calc_params_floating_Transmon_from_files(
            self._output_data_dir,
            capacitance_matrix=self.cap_matrix,
            **kwargs
        )

    @staticmethod
    def calc_params_floating_Transmon_from_files(directory:str, qubit_freq=None, res=None, capacitance_matrix=None, conductor_indices=None, print_all_capacitances=False, C_J=0, Z0_feedline=50):
        '''Static method to calculate key circuit parameters (chi, kappa, g, etc.) for a floating 
        Transmon qubit using the simulated capacitance matrix.

        The design can be either an isolated floating transmon or a floating transmon coupled to a 
        resonator, depending on the number of conductors in the capacitance matrix.

        There are three supported cases:
            3 conductors: 
                isolated floating transmon 
                (pad 1, pad 2, ground)
            4 conductors: 
                floating transmon coupled to resonator
                (ground, pad 1, pad 2, resonator)
            5 conductors: 
                floating transmon coupled to resonator-feedline
                (ground, pad 1, pad 2, resonator, feedline)

        Parameters
        ----------
        directory : str
            Directory containing Palace output files.
        qubit_freq : float
            (Defaults to None) Qubit freqeuncy in Hertz (Hz), used to calculate parameters.
        res : SQDMetal.Utilities.QubitDesigner.ResonatorBase
            (Defaults to None) A resonator object (:class:`~SQDMetal.Utilities.QubitDesigner.ResonatorBase`) is used 
            for the calculation of resonator-related parameters (resonator linewidth, dispersive shift, detuning etc.).
        capacitance_matrix : np.ndarray
            (Defaults to None) Simulated capacitance matrix. If None, fetches from the supplied directory.
        conductor_indeces : dict
            (Defaults to None) Dictionary containing indeces for the conductors corresponding to each index in the capacitance_matrix. 
                If None, the defaults are set to: ``{ 'ground': 0, 'pad1': 1, 'pad2': 2, 'res': 3, 'feed': 4 }``
        print_all_capacitances : bool
            (Defaults to False) If True, prints all capacitances (i.e. pad1-to-ground, pad2-to-ground, etc.).
        C_J : float
            (Defaults to 0) Optional capacitance of the Josephson junction.
        Z0_feedline : float
            (Defaults to 50) Impedence of the feedline.

        Returns
        -------
        params : dict
            Dictionary containing calculated parameters. Some entries may be "N/A" for the provided capacitance matrix.
            The dictionary is of the form:

            Energies and key circuit parameters:

            *   ``'E_C_GHz'`` (`float`):
                Charging energy in gigahertz (GHz)
            *   ``'C_sigma_fF'`` (`float`):
                Total capacitance in femtofarad (fF)
            *   ``'g_MHz'`` (`float`):
                Resonator-Qubit g-coupling in megahertz (MHz)
            *   ``'chi_MHz'`` (`float`):
                Resonator-Qubit :math:`\chi` in megahertz (MHz)
            *   ``'Delta_GHz'`` (`float`):
                Resonator-Qubit detuning in gigahertz (GHz)
            *   ``'anh_MHz'`` (`float`):
                Qubit anharmonicity in megahertz (MHz)
            *   ``'f_q_GHz'`` (`float`):
                Qubit frequency in gigahertz (GHz)
            *   ``'kappa_MHz'`` (`float`):
                Resonator-Qubit kappa (decay rate) in megahertz (MHz)
            *   ``'T1,p_ms'``" (`float`):
                Purcell-limited T1 in milliseconds (ms)

            Individual capacitances and inductances:

            *   ``'C1_ground_fF'`` (`float`):
                C1 to ground capacitance in femtofarad (fF)
            *   ``'C2_ground_fF'`` (`float`):
                C2_ground capacitance in femtofarad (fF)
            *   ``'C1_readout_fF'`` (`float`):
                C1_readout capacitance in femtofarad (fF)
            *   ``'C2_readout_fF'`` (`float`):
                C2_readout capacitance in femtofarad (fF)
            *   ``'C1_feed_fF'`` (`float`):
                C1_feed capacitance in femtofarad (fF)
            *   ``'C2_feed_fF'`` (`float`):
                C2_feed capacitance in femtofarad (fF)
            *   ``'C12_fF'`` (`float`):
                C12 capacitance in femtofarad (fF)
            *   ``'Cres_fF'`` (`float`):
                Cres capacitance in femtofarad (fF)
            *   ``'Lres_pH'`` (`float`):
                Resonator inductance in picohenries (pH)

            Junction parameters:

            *   ``'L_J_nH'`` (`float`):
                Inductance of Josephson junction in nanohenries (nH)
            *   ``'E_J_GHz'`` (`float`):
                Josephson energy in gigahertz (GHz)
            *   ``'I_C_nA'`` (`float`):
                Critical current of Josephson junction in nanoamperes (nA)
            *   ``'E_J/E_C'`` (`float`):
                E_J/E_C ratio required to test if qubit is in transmon regime.                
        '''
        #TODO: Have diagram of capacitances @6biscuits
        # Constants
        e = 1.602176634e-19  # Coulombs
        h = 6.62607015e-34   # J·s     
        # Fetch capacitance matrix (if needed)
        if capacitance_matrix is None:
            raw_data = pd.read_csv(directory + '/terminal-C.csv')
            capacitance_matrix = raw_data
        capMat = capacitance_matrix.values
        if np.allclose(capMat[:,0], np.arange(1, capMat.shape[0]+1)):
            capMat = capMat[:, 1:]
        # detect number of conductors
        num_conductors = len(capMat[0])
        #default indeces
        default_indeces = {
            'ground': 0,
            'pad1': 1,
            'pad2': 2,
            'res': 3,
            'feed': 4
        }
        #Update indices if given
        if conductor_indices:
            default_indeces.update(conductor_indices)
        idx_ground = default_indeces['ground']
        idx_pad1   = default_indeces['pad1']
        idx_pad2   = default_indeces['pad2']
        idx_res    = default_indeces['res']
        idx_feed   = default_indeces['feed']

        # define ground and pad capacitance (all cases)
        C12 = abs(capMat[idx_pad1, idx_pad2])
        C1_ground  = abs(capMat[idx_pad1, idx_ground])
        C2_ground  = abs(capMat[idx_pad2, idx_ground])
        
        # 3 conductor case: qubit only
        if num_conductors == 3:
            C1_readout = 0
            C2_readout = 0
            C1_feed    = 0
            C2_feed    = 0
        # 4 conductors: qubit, res
        elif num_conductors == 4:
            C1_readout = abs(capMat[idx_pad1, idx_res])
            C2_readout = abs(capMat[idx_pad2, idx_res])
            C1_feed    = 0
            C2_feed    = 0
        # 5 conductors: qubit, res, feedline
        elif num_conductors == 5:
            C1_readout = abs(capMat[idx_pad1, idx_res])
            C2_readout = abs(capMat[idx_pad2, idx_res])
            C1_feed    = abs(capMat[idx_pad1, idx_feed])
            C2_feed    = abs(capMat[idx_pad2, idx_feed])
            # Resonator couplings (for kappa)
            C_rf       = abs(capMat[idx_res, idx_feed])
            C_rg       = abs(capMat[idx_res, idx_ground])
            C_rp1      = abs(capMat[idx_res, idx_pad1])
            C_rp2      = abs(capMat[idx_res, idx_pad2])
        else:
            print("Error in capacitance matrix indexing. Exiting.")
            return 0
        
        ## CALCULATIONS
        # Calculate total C
        C_sigma = (((C1_ground+C1_readout)*(C2_ground+C2_readout)) \
                  / (C1_ground+C1_readout+C2_ground+C2_readout)) + C12 + C_J
        # Compute charging energy E_C
        E_C_J = e**2 / (2 * C_sigma)
        E_C_GHz = E_C_J / h / 1e9
        # Calculate qubit parameters (if possible)
        if (res is not None) and (num_conductors >= 4) and (qubit_freq is not None):
            # Use qubit calculator for chi, g, Delta etc.
            params = FloatingTransmonDesigner(res).optimise(
                {
                    "fQubit": qubit_freq,
                    "C_q1": C1_ground,
                    "C_q2": C2_ground,
                    "C_g1": C1_readout,
                    "C_g2": C2_readout,
                    "C_J": C12 + C_J,
                    "chi": (-1e6, -0.2e6),
                    "C_sigma": (1e-18, 1e-6),
                    "Ej/Ec": (50, 200),
                },
                print_results=False
            )
            g_MHz = -params['g_Hz'] * 1e-6
            chi_MHz = -params['chi_Hz'] * 1e-6
            Delta_GHz = params['Delta_Hz'] * 1e-9
            anh_MHz = params['anh_Hz'] * 1e-6
            f_r_Hz = res.get_res_frequency() # Hz
            Cres, Lres = res.get_res_capacitance(), res.get_res_inductance()
            # Calculate kappa, Purcell decay (if applicable)
            if num_conductors >= 5:
                C_r = C_rg + C_rp1 + C_rp2 # res capacitance (except feedline)
                C_c = C_rf # res-feedline mutual capacitance
                # External quality factor and linewidth (frequency units)
                Q_c = C_r / (2 * np.pi * f_r_Hz * Z0_feedline * C_c**2)
                kappa_MHz = f_r_Hz / Q_c * 1e-6
                T1p_ms = ((1 / (kappa_MHz*1e6)) * ((Delta_GHz * 1e9)**2 / (g_MHz*1e6)**2)) * 1e3
            else:
                kappa_MHz = None
                T1p_ms = None
        else:
            print("Could not calculate g, chi, kappa, or Purcell decay rate as there is no readout resonator present, or no defined qubit frequency.")
            g_MHz = 0
            chi_MHz = 0
            Delta_GHz = 0
            anh_MHz = 0
            Cres = 0
            Lres = 0
            kappa_MHz = 0
            T1p_ms = 0
        
        # Calculate junction parameters for target frequency (if applicable)
        if qubit_freq:
            phi_0 = h / (2*e) 
            E_J_GHz = (h * qubit_freq + E_C_J)**2 / (8 * E_C_J) / h / 1e9
            L_J_nH = phi_0**2 / (4 * np.pi**2 * h * E_J_GHz)
        else:
            E_J_GHz = 0
            L_J_nH = 0

        if print_all_capacitances:
            print("Capacitance Results")
            print("--------------------")
            print(f"{'C1-Ground':<16s} = {C1_ground * 1e15:>10.3f} fF")
            print(f"{'C2-Ground':<16s} = {C2_ground * 1e15:>10.3f} fF")
            print(f"{'C1-Readout':<16s} = {C1_readout * 1e15:>10.3f} fF")
            print(f"{'C2-Readout':<16s} = {C2_readout * 1e15:>10.3f} fF")
            print(f"{'C1-Feedline':<16s} = {C1_feed * 1e15:>10.3f} fF")
            print(f"{'C2-Feedline':<16s} = {C2_feed * 1e15:>10.3f} fF")
            print(f"{'Mutual C12':<16s} = {C12 * 1e15:>10.3f} fF")
            print(f"{'C_sigma':<16s} = {C_sigma * 1e15:>10.3f} fF\n")

        # Print Readout resonator parameters
        if (res is not None) and num_conductors >= 4:
            print("Readout Resonator")
            print("-------------------")
            print(f"{'f_res':<16s} = {f_r_Hz * 1e-9:>10.3f} GHz")
            print(f"{'L_res':<16s} = {Lres * 1e12:>10.3f} pH")
            print(f"{'C_res':<16s} = {Cres * 1e15:>10.3f} fF")
            if num_conductors >= 5:
                print(f"{'C_res-Feedline':<16s} = {C_c * 1e15:>10.3f} fF")
            print()

        # Print charging energy
        print("Circuit Parameters")
        print("-------------------")
        if qubit_freq is not None:
            print(f"{'f_qubit':<16s} = {qubit_freq * 1e-9:>10.3f} GHz")
        print(f"{'E_C':<16s} = {E_C_GHz * 1e3:>10.3f} MHz")
        print(f"{'g':<16s} = {g_MHz if isinstance(g_MHz, str) else f'{g_MHz:>10.3f}'} MHz")
        print(f"{'chi':<16s} = {chi_MHz if isinstance(chi_MHz, str) else f'{chi_MHz:>10.3f}'} MHz")
        print(f"{'Delta':<16s} = {Delta_GHz if isinstance(Delta_GHz, str) else f'{Delta_GHz:>10.3f}'} GHz")
        print(f"{'Anharmonicity':<16s} = {anh_MHz if isinstance(anh_MHz, str) else f'{anh_MHz:>10.3f}'} MHz")
        if kappa_MHz is not 'N/A':
            print(f"{'kappa':<16s} = {kappa_MHz:>10.3f} MHz")
            print(f"{'T_1,P':<16s} = {T1p_ms:>10.3f} ms")
        print()

        # Print required junction parameters for desired frequency
        if qubit_freq is not None:
            print(f"Target Junction Parameters\n (for f_qubit = {qubit_freq * 1e-9:.1f} GHz)")
            print("--------------------------")
            print(f"{'E_J':<16s} = {E_J_GHz:>10.3f} GHz")
            print(f"{'L_J':<16s} = {L_J_nH:>10.3f} nH\n")
            print(f"{'E_J/E_C':<16s} = {E_J_GHz/E_C_GHz:>10.3f}")

        return {
            # Energies and key circuit parameters
            "E_C_GHz": E_C_GHz,
            "C_sigma_fF": C_sigma * 1e15,
            "g_MHz": g_MHz,
            "chi_MHz": chi_MHz,
            "Delta_GHz": Delta_GHz,
            "anh_MHz": anh_MHz,
            "f_q_GHz": qubit_freq * 1e-9,
            "kappa_MHz": kappa_MHz,
            "T1,p_ms": T1p_ms,

            # Individual capacitances and inductances
            "C1_ground_fF": C1_ground * 1e15,
            "C2_ground_fF": C2_ground * 1e15,
            "C1_readout_fF": C1_readout * 1e15,
            "C2_readout_fF": C2_readout * 1e15,
            "C1_feed_fF": C1_feed * 1e15,
            "C2_feed_fF": C2_feed * 1e15,
            "C12_fF": C12 * 1e15,
            "Cres_fF": Cres * 1e15,
            "Lres_pH": Lres * 1e12,

            # Junction parameters
            'L_J_nH': L_J_nH,
            'E_J_GHz': E_J_GHz,
            'I_C_nA': E_J_GHz/0.495,
            'E_J/E_C': E_J_GHz/E_C_GHz
        }

    def calc_params_2_floating_transmons(self, **kwargs):
        """
        Class method mirroring :func:`~SQDMetal.PALACE.Capacitance_Simulation.PALACE_Capacitance_Simulation.calc_params_2_floating_transmons_from_files` (without the ``directory`` parameter as it is inferred from the current simulation)

        This function must be run after calling :func:`~SQDMetal.PALACE.Model.PALACE_Model_Base.run`.
        """
        return PALACE_Capacitance_Simulation.calc_params_2_floating_transmons_from_files(
            self._output_data_dir,
            capacitance_matrix=self.cap_matrix,
            **kwargs
        )

    @staticmethod
    def calc_params_2_floating_transmons_from_files(directory:str, qubit_freq=None, res=None, capacitance_matrix=None, conductor_indices=None, print_all_capacitances=False, C_J=0, Z0_feedline=50):
        '''Static method to calculate key circuit parameters (chi, kappa, g, etc.) for 2 Floating Transmons qubits using the simulated capacitance matrix.

        The design can be either 2 isolated floating transmons or 2 floating transmons coupled to a 
        resonator, depending on the number of conductors in the capacitance matrix.

        There are three supported cases:
            5 conductors: 
                2 floating transmons coupled through mutual capacitance
                (ground, pad 1, pad 2, pad 3, pad 4)
            6 conductors:
                2 floating transmons coupled through mutual capacitance, and coupled to resonator
                (ground, pad 1, pad 2, pad 3, pad 4, resonator)

        Parameters
        ----------
        directory : str
            Directory containing Palace output files.
        qubit_freq : List of floats
            (Defaults to None) Qubits freqeuncy in Hertz (Hz), used to calculate parameters.
        res : SQDMetal.Utilities.QubitDesigner.ResonatorBase
            (Defaults to None) A resonator object (:class:`~SQDMetal.Utilities.QubitDesigner.ResonatorBase`) is used 
            for the calculation of resonator-related parameters (resonator linewidth, dispersive shift, detuning etc.).
        capacitance_matrix : np.ndarray
            (Defaults to None) Simulated capacitance matrix. If None, fetches from the supplied directory.
        conductor_indeces : dict
            (Defaults to None) Dictionary containing indeces for the conductors corresponding to each index in the capacitance_matrix. 
                If None, the defaults are set to: ``{ 'ground': 0, 'pad1': 1, 'pad2': 2, 'pad 3': 3, 'pad 4': 4 }``
        print_all_capacitances : bool
            (Defaults to False) If True, prints all capacitances (i.e. pad1-to-ground, pad2-to-ground, etc.).
        C_J : float
            (Defaults to 0) Optional capacitance of the Josephson junction.

        Returns
        -------
        params : dict
            Dictionary containing calculated parameters. Some entries may be "N/A" for the provided capacitance matrix.
            The dictionary is of the form:

            Energies and key circuit parameters:

            *   ``'E_C1_GHz'`` (`float`):
                Charging energy in gigahertz (GHz) for qubit 1
            *   ``'E_C2_GHz'`` (`float`):
                Charging energy in gigahertz (GHz) for qubit 2
            *   ``'E_Cc_GHz'`` (`float`):
                Charging energy in gigahertz (GHz) for the coupled qubits
            *   ``'C1_sigma_fF'`` (`float`):
                Total capacitance in femtofarad (fF) for qubit 1
            *   ``'C2_sigma_fF'`` (`float`):
                Total capacitance in femtofarad (fF) for qubit 2
            *   ``'g1_MHz'`` (`float`):
                Resonator-Qubit1 g-coupling in megahertz (MHz)
            *   ``'g2_MHz'`` (`float`):
                Resonator-Qubit2 g-coupling in megahertz (MHz)
            *   ``'g12_MHz'`` (`float`):
                Qubit-Qubit g-coupling in megahertz (MHz)
            *   ``'chi_MHz'`` (`float`):
                Resonator-Qubit :math:`\chi` in megahertz (MHz)
            *   ``'Delta_GHz'`` (`float`):
                Resonator-Qubit detuning in gigahertz (GHz)
            *   ``'anh_MHz'`` (`float`):
                Qubit anharmonicity in megahertz (MHz)
            *   ``'f_q_GHz'`` (`float`):
                Qubit frequency in gigahertz (GHz)
            *   ``'kappa_MHz'`` (`float`):
                Resonator-Qubit kappa (decay rate) in megahertz (MHz)
            *   ``'T1,p_ms'``" (`float`):
                Purcell-limited T1 in milliseconds (ms)

            Individual capacitances and inductances:

            *   ``'C1_ground_fF'`` (`float`):
                C1 to ground capacitance in femtofarad (fF)
            *   ``'C2_ground_fF'`` (`float`):
                C2_ground capacitance in femtofarad (fF)
            *   ``'C3_ground_fF'`` (`float`):
                C3 to ground capacitance in femtofarad (fF)
            *   ``'C4_ground_fF'`` (`float`):
                C4 to ground capacitance in femtofarad (fF)
            *   ``'C1_res_fF'`` (`float`):
                C1_res capacitance in femtofarad (fF)
            *   ``'CQ1 fF'`` (`float`):
                Effective capacitance of qubit 1 (capacitance between pad 1 and pad 2, including their capacitance to ground) in femtofarad (fF)
            *   ``'CQ2 fF'`` (`float`):
                Effective capacitance of qubit 2 (capacitance between pad 3 and pad 4, including their capacitance to ground) in femtofarad (fF)
            *   ``'C2_res_fF'`` (`float`):
                C2_res capacitance in femtofarad (fF)
            *   ``'C3_res_fF'`` (`float`):
                C3_res capacitance in femtofarad (fF)
            *   ``'C4_res_fF'`` (`float`):
                C4_res capacitance in femtofarad (fF)
            *   ``'CQ1Q2 fF'`` (`float`):
                Coupling capacitance between qubit 1 and qubit 2 (capacitance between the two qubits excluding coupling through the resonator) in femtofarad (fF)
            *   ``'Cres_fF'`` (`float`):
                Cres capacitance in femtofarad (fF)
            *   ``'Lres_pH'`` (`float`):
                Resonator inductance in picohenries (pH)

            Junction parameters:

            *   ``'L1_J_nH'`` (`float`):
                Inductance of Josephson junction in nanohenries (nH)
            *  ``'L2_J_nH'`` (`float`):
                Inductance of Josephson junction in nanohenries (nH)
            *   ``'E1_J_GHz'`` (`float`):
                Josephson energy in gigahertz (GHz) for qubit 1
            *   ``'E2_J_GHz'`` (`float`):
                Josephson energy in gigahertz (GHz) for qubit 2
            *   ``'I_C_nA'`` (`float`):
                Critical current of Josephson junction in nanoamperes (nA)
            *   ``'E_J/E_C'`` (`float`):
                E_J/E_C ratio required to test if qubit is in transmon regime.                
        '''
        #TODO: Have diagram of capacitances @6biscuits
        # Constants
        e = 1.602176634e-19  # Coulombs
        h = 6.62607015e-34   # J·s     
        hbar = h / (2 * np.pi)
        # Fetch capacitance matrix (if needed)
        if capacitance_matrix is None:
            raw_data = pd.read_csv(directory + '/terminal-C.csv')
            capacitance_matrix = raw_data
        capMat = capacitance_matrix.values
        if np.allclose(capMat[:,0], np.arange(1, capMat.shape[0]+1)):
            capMat = capMat[:, 1:]
        # detect number of conductors
        num_conductors = len(capMat[0])
        #default indeces
        default_indeces = {
            'ground': 0,
            'pad1': 1,
            'pad2': 2,
            'res': 3,
            'pad3': 4,
            'pad4': 5,
            
         }
        idx_res = default_indeces['res']
        #Update indices if given
        if conductor_indices:
            default_indeces.update(conductor_indices)
        idx_ground = default_indeces['ground']
        idx_pad1   = default_indeces['pad1']
        idx_pad2   = default_indeces['pad2']
        idx_pad3    = default_indeces['pad3']
        idx_pad4   = default_indeces['pad4']

        # define ground and pad capacitance (all cases)
        C12 = abs(capMat[idx_pad1, idx_pad2])    # capacitance between pad 1 and pad 2 of qubit 1
        C34 = abs(capMat[idx_pad3, idx_pad4])    # capacitance between pad 3 and pad 4 of qubit 2
        C13 = abs(capMat[idx_pad1, idx_pad3])    # capacitance between pad 1 and pad 3 (qubit 1 to qubit 2)
        C14 = abs(capMat[idx_pad1, idx_pad4])    # capacitance between pad 1 and pad 4 (qubit 1 to qubit 2)
        C23 = abs(capMat[idx_pad2, idx_pad3])    # capacitance between pad 2 and pad 3 (qubit 1 to qubit 2)
        C24 = abs(capMat[idx_pad2, idx_pad4])    # capacitance between pad 2 and pad 4 (qubit 1 to qubit 2)
        C1_ground  = abs(capMat[idx_pad1, idx_ground])
        C2_ground  = abs(capMat[idx_pad2, idx_ground])
        C3_ground  = abs(capMat[idx_pad3, idx_ground])
        C4_ground  = abs(capMat[idx_pad4, idx_ground])
        # Coupling capacitance between qubit 1 and qubit 2 (capacitance between the two qubits excluding coupling through the resonator)
        CQ1Q2 = C13 + C14 + C23 + C24                            ####### TO BE VERIFIED, CHANGED

        C1_res = abs(capMat[idx_pad1, idx_res])
        C2_res = abs(capMat[idx_pad2, idx_res])
        C3_res = abs(capMat[idx_pad3, idx_res])
        C4_res = abs(capMat[idx_pad4, idx_res])
        #New CQ1 and CQ2 including coupling through the resonator
        CQ1 = C12 + ((C1_ground+C1_res)*(C2_ground+C2_res))/(C1_ground + C2_ground + C1_res + C2_res)
        CQ2 = C34 + ((C3_ground+C3_res)*(C4_ground+C4_res))/(C3_ground + C4_ground + C3_res + C4_res)
        # Resonator couplings (for kappa)
        C_r       = abs(capMat[idx_res, idx_ground])
        Cg1 = (C2_ground*C1_res-C1_ground*C2_res)/(C1_ground + C2_ground + C1_res + C2_res) 
        Cg2 = (C3_ground*C4_res-C4_ground*C3_res)/(C3_ground + C4_ground + C3_res + C4_res) 
        # Global capacitance calculation
        C_sigma = CQ1Q2*(Cg1+Cg2)*(CQ1+CQ2) + C_r*CQ1*CQ2 + CQ1Q2*C_r*(Cg1+Cg2+CQ1+CQ2) + Cg1*CQ2*(CQ1+C_r) + Cg2*CQ1*(CQ2+C_r) + Cg1*Cg2*(CQ1+CQ2+C_r)
        C_1_prime = (C_sigma)/(-(Cg2**2) + (CQ1Q2+Cg2+CQ2)*(Cg1+Cg2+C_r))
        C_2_prime = (C_sigma)/(-(Cg1**2) + (CQ1Q2+Cg1+CQ1)*(Cg1+Cg2+C_r))
        C_r_prime = (C_sigma)/((Cg1+CQ1)*(Cg2+CQ2) + CQ1Q2*(Cg1+Cg2+CQ1+CQ2))
        Cg1_prime = (C_sigma)/(2*(Cg1*(CQ2+Cg2)+CQ1Q2*(Cg1+Cg2)))
        Cg2_prime = (C_sigma)/(2*(Cg2*(CQ1+Cg1)+CQ1Q2*(Cg1+Cg2)))
        CQ1Q2_prime =  (C_sigma)/(2*(Cg1*Cg2 + CQ1Q2*(Cg1+Cg2+C_r)))
        
        ## CALCULATIONS
        # Compute charging energy E_C1 and E_C2
        E_C1_J = e**2 / (2 * C_1_prime)
        E_C1_GHz = E_C1_J / h / 1e9
        E_C2_J = e**2 / (2 * C_2_prime)
        E_C2_GHz = E_C2_J / h / 1e9
        E_Cc_J = e**2 / (2 * CQ1Q2_prime)
        E_Cc_GHz = E_Cc_J / h / 1e9
        E_Cr = e**2 / (2 * C_r_prime)
        E_Cr_GHz = E_Cr / h / 1e9
        # Calculate qubit parameters (if possible)
        if (res is not None) and (qubit_freq is not None):    ###### TO BE IMPLEMENTED
            # Use qubit calculator for chi, g, Delta etc.
            params = FloatingTransmonDesigner(res).optimise(
                {
                    "fQubit": qubit_freq[0],
                    "C_q1": C1_ground,
                    "C_q2": C2_ground,
                    "C_g1": C1_res,
                    "C_g2": C2_res,
                    "C_J": C12 ,
                    "chi": (-1e6, -0.2e6),
                    "C_sigma": (1e-18, 1e-6),
                    "Ej/Ec": (50, 200),
                },
                print_results=False
            )
            g1_MHz = -params['g_Hz'] * 1e-6
            chi_MHz = -params['chi_Hz'] * 1e-6
            Delta_GHz = params['Delta_Hz'] * 1e-9
            anh_MHz = params['anh_Hz'] * 1e-6
            f_r_Hz = res.get_res_frequency() # GHz
            Cres, Lres = res.get_res_capacitance(), res.get_res_inductance()
            kappa_MHz = None
            T1p_ms = None
            params = FloatingTransmonDesigner(res).optimise(
                {
                    "fQubit": qubit_freq[1],
                    "C_q1": C3_ground,
                    "C_q2": C4_ground,
                    "C_g1": C3_res,
                    "C_g2": C4_res,
                    "C_J": C34 ,
                    "chi": (-1e6, -0.2e6),
                    "C_sigma": (1e-18, 1e-6),
                    "Ej/Ec": (50, 200),
                },
                print_results=False
            )
            g2_MHz = -params['g_Hz'] * 1e-6

        else:
            print("Could not calculate g, chi, kappa, or Purcell decay rate as there is no readout resonator present, or no defined qubit frequency.")
            g1_MHz = 0
            g2_MHz = 0
            chi_MHz = 0
            Delta_GHz = 0
            anh_MHz = 0
            Cres = 0
            Lres = 0
            kappa_MHz = 0
            T1p_ms = 0
        
        # Calculate junction parameters for target frequency (if applicable)
        if qubit_freq:
            phi_0 = h / (2*e) 
            E1_J = (h * qubit_freq[0] + E_C1_J)**2 / (8 * E_C1_J)
            E1_J_GHz = (h * qubit_freq[0] + E_C1_J)**2 / (8 * E_C1_J) / h / 1e9
            L1_J_nH = phi_0**2 / (4 * np.pi**2 * h * E1_J_GHz)
            E2_J = (h * qubit_freq[1] + E_C2_J)**2 / (8 * E_C2_J)
            E2_J_GHz = (h * qubit_freq[1] + E_C2_J)**2 / (8 * E_C2_J) / h / 1e9
            L2_J_nH = phi_0**2 / (4 * np.pi**2 * h * E2_J_GHz)
        else:
            E1_J_GHz = 0
            L1_J_nH = 0
            E2_J_GHz = 0
            L2_J_nH = 0
        
        # Compute qubit-qubit coupling J12, omega_1, omega_2 (if possible) + Z_onoff_ratio

        if qubit_freq:
                # Calculate qubit-qubit coupling J12 (in Hz)
              g1 = ((e/np.sqrt(hbar))/(2*Cg1_prime))*(((E1_J_GHz*C_r_prime)/(8*E_C1_GHz*Lres))**(1/4))
              g2 = ((e/np.sqrt(hbar))/(2*Cg2_prime))*(((E2_J_GHz*C_r_prime)/(8*E_C2_GHz*Lres))**(1/4))
              g12_GHz = ((e**2)/(2*CQ1Q2_prime*hbar))*(((E1_J_GHz*E2_J_GHz)/(4*E_C1_GHz*E_C2_GHz))**(1/4))*1e-9
              g12_MHz = g12_GHz * 1e3
              #Calculate qubit and resonator frequencies including coupling (in Hz)
              omega_1_GHz = (np.sqrt(8*E_C1_GHz*E1_J_GHz) - E_C1_GHz)
              omega_2_GHz = (np.sqrt(8*E_C2_GHz*E2_J_GHz) - E_C2_GHz)
              omega_1_MHz = omega_1_GHz * 1e3
              omega_2_MHz = omega_2_GHz * 1e3
              omega_r_GHz = 1/(np.sqrt(Lres*C_r_prime)) * 1e-9
              omega_r_MHz = omega_r_GHz * 1e3
              

             # Calculate Z on/off ratio (assuming on state is when qubits are resonant, off state is when they are detuned by 500 MHz)
              Delta = abs(E1_J_GHz - E2_J_GHz) # detuning in Hz
              Z_onoff_ratio = Delta / (g12_GHz) if g12_GHz != 0 else "N/A"

        
        if print_all_capacitances:
            print("Capacitance Results")
            print("--------------------")
            print(f"{'C1-Ground':<16s} = {C1_ground * 1e15:>10.3f} fF")
            print(f"{'C2-Ground':<16s} = {C2_ground * 1e15:>10.3f} fF")
            print(f"{'C3-Ground':<16s} = {C3_ground * 1e15:>10.3f} fF")
            print(f"{'C4-Ground':<16s} = {C4_ground * 1e15:>10.3f} fF")
            print(f"{'CQ1':<16s} = {CQ1 * 1e15:>10.3f} fF")
            print(f"{'CQ2':<16s} = {CQ2 * 1e15:>10.3f} fF")
            print(f"{'Cg1':<16s} = {Cg1 * 1e15:>10.3f} fF")
            print(f"{'Cg2':<16s} = {Cg2 * 1e15:>10.3f} fF")
            print(f"{'C_r':<16s} = {C_r * 1e15:>10.3f} fF")
            print(f"{'CQ1Q2':<16s} = {CQ1Q2 * 1e15} fF")
            print(f"{'Mutual C12_prime':<16s} = {CQ1Q2_prime * 1e15:>10.3f} fF")
            print(f"{'C1_prime':<16s} = {C_1_prime * 1e15:>10.3f} fF")
            print(f"{'C2_prime':<16s} = {C_2_prime * 1e15:>10.3f} fF")
            print(f"{'C_r_prime':<16s} = {C_r_prime * 1e15:>10.3f} fF")
            print(f"{'Cg1_prime':<16s} = {Cg1_prime * 1e15:>10.3f} fF")
            print(f"{'Cg2_prime':<16s} = {Cg2_prime * 1e15:>10.3f} fF\n")

        # Print Readout resonator parameters
        if (res is not None) and num_conductors == 6:
            print("Readout Resonator")
            print("-------------------")
            print(f"{'f_res':<16s} = {f_r_Hz * 1e-9:>10.3f} GHz")
            print(f"{'L_res':<16s} = {Lres * 1e12:>10.3f} pH")
            print(f"{'C_res':<16s} = {Cres * 1e15:>10.3f} fF")

        # Print charging energy
        print("Circuit Parameters")
        print("-------------------")
        if qubit_freq is not None:
            print(f"{'f_qubit1':<16s} = {qubit_freq[0] * 1e-9:>10.3f} GHz")
            print(f"{'f_qubit2':<16s} = {qubit_freq[1] * 1e-9:>10.3f} GHz")
        print(f"{'E_C1':<16s} = {E_C1_GHz * 1e3:>10.3f} MHz")
        print(f"{'E_C2':<16s} = {E_C2_GHz * 1e3:>10.3f} MHz")
        print(f"{'E_Cc':<16s} = {E_Cc_GHz * 1e3:>10.3f} MHz")
        print(f"{'g1':<16s} = {g1_MHz if isinstance(g1_MHz, str) else f'{g1_MHz:>10.3f}'} MHz")
        print(f"{'g2':<16s} = {g2_MHz if isinstance(g2_MHz, str) else f'{g2_MHz:>10.3f}'} MHz")
        print(f"{'chi':<16s} = {chi_MHz if isinstance(chi_MHz, str) else f'{chi_MHz:>10.3f}'} MHz")
        print(f"{'Delta':<16s} = {Delta_GHz if isinstance(Delta_GHz, str) else f'{Delta_GHz:>10.3f}'} GHz")
        print(f"{'Anharmonicity':<16s} = {anh_MHz if isinstance(anh_MHz, str) else f'{anh_MHz:>10.3f}'} MHz")
        #if kappa_MHz is not 'N/A':
        #    print(f"{'kappa':<16s} = {kappa_MHz:>10.3f} MHz")
        #    print(f"{'T_1,P':<16s} = {T1p_ms:>10.3f} ms")
        print()

        # Print required junction parameters for desired frequency
        if qubit_freq is not None:
            print(f"Target Junction Parameters\n (for f_qubit1 = {qubit_freq[0] * 1e-9:.1f} GHz)")
            print("--------------------------")
            print(f"{'E_J1':<16s} = {E1_J_GHz:>10.3f} GHz")
            print(f"{'L_J1':<16s} = {L1_J_nH:>10.3f} nH")
            print(f"{'E_J1/E_C1':<16s} = {E1_J_GHz/E_C1_GHz:>10.3f}\n")
            print(f"Target Junction Parameters\n (for f_qubit2 = {qubit_freq[1] * 1e-9:.1f} GHz)")
            print("--------------------------")
            print(f"{'E_J2':<16s} = {E2_J_GHz:>10.3f} GHz")
            print(f"{'L_J2':<16s} = {L2_J_nH:>10.3f} nH")
            print(f"{'E_J2/E_C2':<16s} = {E2_J_GHz/E_C2_GHz:>10.3f}\n")
            print(f"{'Qubit-Qubit coupling g12':<16s} = {g12_GHz:>10.3f} GHz")
            print(f"{'Qubit 1 g1 coupling':<16s} = {g1*1e-6:>10.3f} MHz")
            print(f"{'Qubit 2 g2 coupling':<16s} = {g2*1e-6:>10.3f} MHz")
            print(f"{'Qubit 1 frequency (uncoupled)':<16s} = {omega_1_GHz:>10.3f} GHz")
            print(f"{'Qubit 2 frequency (uncoupled)':<16s} = {omega_2_GHz:>10.3f} GHz")
            print(f"{'Resonator frequency (uncoupled)':<16s} = {omega_r_GHz:>10.3f} GHz")
            print(f"{'Z on/off ratio':<16s} = {Z_onoff_ratio if isinstance(Z_onoff_ratio, str) else f'{Z_onoff_ratio:>10.3f}'}")

        return {
            # Energies and key circuit parameters
            "E_C1_GHz": E_C1_GHz,
            "E_C2_GHz": E_C2_GHz,
            "C_1_fF": C_1_prime * 1e15,
            "C_2_fF": C_2_prime * 1e15,
            "g1_MHz": g1_MHz,
            "chi_MHz": chi_MHz,
            "Delta_GHz": Delta_GHz,
            "anh_MHz": anh_MHz,
            "f_q1_GHz": qubit_freq[0] * 1e-9,
            "f_q2_GHz": qubit_freq[1] * 1e-9,
            "kappa_MHz": kappa_MHz,
            "T1,p_ms": T1p_ms,

            # Individual capacitances and inductances
            "C1_ground_fF": C1_ground * 1e15,
            "C2_ground_fF": C2_ground * 1e15,
            "C3_ground_fF": C3_ground * 1e15,
            "C4_ground_fF": C4_ground * 1e15,
            "CQ1_fF": CQ1 * 1e15,
            "CQ2_fF": CQ2 * 1e15,
            "C_res_ground": C_r * 1e15,
            "C12_fF": CQ1Q2_prime * 1e15,
            "Cres_fF": Cres * 1e15,
            "Lres_pH": Lres * 1e12,

            # Junction parameters
            'L_J1_nH': L1_J_nH,
            'E_J1_GHz': E1_J_GHz,
            'I_C1_nA': E1_J_GHz/0.495,
            'L_J2_nH': L2_J_nH,
            'E_J2_GHz': E2_J_GHz,
            'I_C2_nA': E2_J_GHz/0.495,
            'E_J1/E_C1': E1_J_GHz/E_C1_GHz,
            'E_J2/E_C2': E2_J_GHz/E_C2_GHz,
            'g12_MHz': g12_MHz,
            'g1_MHz': g1,
            'g2_MHz': g2,
            'omega_1_MHz': omega_1_MHz,
            'omega_2_MHz': omega_2_MHz,
            'omega_r_MHz': omega_r_MHz,
            'Z_onoff_ratio': Z_onoff_ratio
            
        }