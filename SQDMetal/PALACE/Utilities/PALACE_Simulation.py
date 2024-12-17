from GMSH_Geometrey_Builder import GMSH_Geometry_Builder
from GMSH_Mesh_Builder import GMSH_Mesh_Builder
from Eigenmode_Simulation import Eigenmode_Simulation
from Driven_Simulation import Driven_Simulation
from Capacitance_Simulation import Capacitance_Simulation
from Simulation_Files_Builder import Simulation_Files_Builder
import gmsh


class PALACE_Simulation:

    using_hpc = False

    def __init__(self, simulation_type, name, design, user_options, ports = [], hpc_options = {}):

        self.simulation_type = simulation_type
        self.design = design
        self.name = name
        self.user_options = user_options
        self.ports = ports
        self.hpc_options = hpc_options

    def run_simulation(self):
        
        #create gmsh geometry builder object and construct qiskit metal design in gmsh
        GGB = GMSH_Geometry_Builder(self.design, self.simulation_type, self.ports)
        _, _, _, dielectric_cutouts, ports_dict, metal_cap_physical_group, metal_cap_names, jj_dict = GGB.construct_geometry_in_GMSH()

        #create simulation object
        if self.simulation_type == 'Eigenmode':
            eigen_sim = Eigenmode_Simulation(self.name, ports_dict, self.user_options, jj_dict, self.hpc_options)
            physical_groups = eigen_sim.prepare_simulation()
            sim_config_file = eigen_sim.create_sim_config_file(physical_groups)

        elif self.simulation_type == 'Driven':
            driven_sim = Driven_Simulation(self.name, ports_dict, self.user_options, jj_dict, self.hpc_options)
            physical_groups = driven_sim.prepare_simulation()
            sim_config_file = driven_sim.create_sim_config_file(physical_groups)

        elif self.simulation_type == 'Capacitance':
            cap_sim = Capacitance_Simulation(self.name, metal_cap_physical_group, metal_cap_names, self.user_options, self.hpc_options)
            physical_groups, metals = cap_sim.prepare_simulation()
            sim_config_file = cap_sim.create_sim_config_file(physical_groups)

        else:
            raise Exception("Simulation type incorrectly specified. Simulation type must be either 'Eigenmode', 'Driven', or 'Capacitance'.")

        #create gmsh mesh builder object and build the mesh for the design - currently using dielectric cutouts for fine meshing, can change to metals
        GMB = GMSH_Mesh_Builder(dielectric_cutouts, self.user_options)
        GMB.build_mesh()

        #create Simulation Files Builder object to handle creation of config file and mesh file
        SFB = Simulation_Files_Builder(self.name, self.user_options, sim_config_file, self.hpc_options)
        SFB.create_simulation_files()

        #open gmsh
        gmsh.fltk.run()

        



    