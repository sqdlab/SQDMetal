import qiskit_metal as metal
from qiskit_metal import designs, draw
from qiskit_metal import MetalGUI, Dict, Headings
import pyEPR as epr
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
from qiskit_metal.qlibrary.qubits.transmon_pocket import TransmonPocket

from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
from SQDMetal.Utilities.Materials import MaterialInterface

import shutil
import unittest

class TestPalace(unittest.TestCase):
    ERR_TOL = 5e-13
    
    def initialise(self):
        self._folder_path = 'TestDesign'

    def cleanup(self):
        shutil.rmtree(self._folder_path)

    def test_TransmonResDesign(self):
        return
        self.initialise()

        design = designs.DesignPlanar({}, True)
        design.chips.main.size.center_x = '0.5mm'
        design.chips.main.size.center_y = '0.1mm'
        design.chips.main.size['size_x'] = '2.8mm'
        design.chips.main.size['size_y'] = '2mm'
        q1 = TransmonPocket(design, 'Q1', options = dict(
            pad_width = '425 um',
            pocket_height = '650um',
            connection_pads=dict(
                readout = dict(loc_W=+1,loc_H=+1, pad_width='200um')
            )))

        otg = OpenToGround(design, 'open_to_ground', options=dict(pos_x='1.75mm',  pos_y='0um', orientation='0'))
        RouteMeander(design, 'readout',  Dict(
                total_length='6 mm',
                hfss_wire_bonds = True,
                fillet='90 um',
                lead = dict(start_straight='100um'),
                pin_inputs=Dict(
                start_pin=Dict(component='Q1', pin='readout'),
                end_pin=Dict(component='open_to_ground', pin='open')), ))
        #
        #
        #Eigenmode Simulation Options
        user_defined_options = {
                        "mesh_refinement":  0,                             #refines mesh in PALACE - essetially divides every mesh element in half
                        "dielectric_material": "silicon",                  #choose dielectric material - 'silicon' or 'sapphire'
                        "starting_freq": 5.5e9,                              #starting frequency in Hz 
                        "number_of_freqs": 3,                              #number of eigenmodes to find
                        "solns_to_save": 3,                                #number of electromagnetic field visualizations to save
                        "solver_order": 2,                                 #increasing solver order increases accuracy of simulation, but significantly increases sim time
                        "solver_tol": 1.0e-8,                              #error residual tolerance for iterative solver
                        "solver_maxits": 200,                              #number of solver iterations
                        "fillet_resolution":12,                            #number of vertices per quarter turn on a filleted path
                        "palace_dir":"~/spack/opt/spack/linux-ubuntu24.04-zen2/gcc-13.3.0/palace-develop-36rxmgzatchgymg5tcbfz3qrmkf4jnmj/bin/palace",#"PATH/TO/PALACE/BINARY",
                        "num_cpus": 16
                        }
        #Create the Palace Eigenmode simulation
        eigen_sim = PALACE_Eigenmode_Simulation(name =self._folder_path,                                              #name of simulation
                                                metal_design = design,                                      #feed in qiskit metal design
                                                sim_parent_directory = "",                                  #choose directory where mesh file, config file and HPC batch file will be saved
                                                mode = 'simPC',                                             #choose simulation mode 'HPC' or 'simPC'                                          
                                                meshing = 'GMSH',                                           #choose meshing 'GMSH' or 'COMSOL'
                                                user_options = user_defined_options,                        #provide options chosen above
                                                create_files = True)                                        #create mesh and config files

        #add in metals from the first layer
        eigen_sim.add_metallic(1)
        #add ground plane into simulation
        eigen_sim.add_ground_plane()
        #Create a lumped element port for the Josephson junction and assign Jospehson inductance and junction capacitance
        eigen_sim.create_port_JosephsonJunction('Q1', L_J = 11e-9, C_J = 0e-15)
        #Fine mesh the qubit and resonator - min_size/max_size is the min/max mesh element size
        eigen_sim.fine_mesh_around_comp_boundaries(['Q1'], min_size=12e-6, max_size=100e-6, taper_dist_min=10e-6, metals_only=True)
        eigen_sim.fine_mesh_along_path(qObjName='readout', dist_resolution=10e-6, min_size=12e-6, max_size=150e-6, taper_dist_min=10e-6)
        #Lossy participatoin ratios calculated for MA, SA and MS
        eigen_sim.setup_EPR_interfaces(metal_air=MaterialInterface('Aluminium-Vacuum'), substrate_air=MaterialInterface('Silicon-Vacuum'), substrate_metal=MaterialInterface('Silicon-Aluminium'))
        #Prepare the simulation files - mesh file (.msh) and config file (.json)
        eigen_sim.prepare_simulation()

        self.cleanup()

if __name__ == '__main__':
    TestPalace().test_TransmonResDesign()
    unittest.main()
