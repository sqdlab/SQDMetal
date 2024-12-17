import gmsh
from Utilities.Materials import Material

class Capacitance_Simulation():

    def __init__(self, name, metal_cap_physical_group, metal_cap_names, user_options, hpc_options = []):
        self.name = name
        self.metal_cap_physical_group = metal_cap_physical_group
        self.metal_cap_names = metal_cap_names
        self.user_options = user_options
        self.hpc_options = hpc_options

    def prepare_simulation(self):
        
        #dictionary for physical groups/bundary conditions
        physical_groups = {}

        physical_group_ids = gmsh.model.get_physical_groups()
        for _,group in enumerate(physical_group_ids):
            group_name = gmsh.model.get_physical_name(group[0], group[1])
            if group_name == 'dielectric_gaps':
                physical_groups['dielectric_gaps'] = group[1]
            elif group_name == 'ground_plane':
                physical_groups['ground_plane'] = group[1]
            elif group_name == 'dielectric_substrate':
                physical_groups['dielectric_substrate'] = group[1]
            elif group_name == 'air_box':
                physical_groups['air_box'] = group[1]
            elif group_name == 'far_field':
                physical_groups['far_field']  = group[1]
        
        #add individual metals to metals dictionary
        metals = {}
        if len(self.metal_cap_names) != 0:
            for m,metal_name in enumerate(self.metal_cap_names):
                metals[metal_name] = self.metal_cap_physical_group[m]

        return physical_groups, metals
    

    def create_sim_config_file(self, physical_groups):
        
        #get the material propoerties for the dielectric specified by the user
        dielectric = Material(self.user_options['dielectric_material'])

        #prepare different metals for capacitance simulation
        metal_terminals = self._process_metals_for_config_file()

        #if hpc_options is not empty then supply correct location for mesh file and location to output simulation files
        if self.hpc_options:
            mesh_location = self.hpc_options['input_files_location'] + self.name + '/' + self.name + '.msh'
            output = self.hpc_options['output_files_location'] + self.name
        else:
            mesh_location = self.user_options['sim_directory'] + self.name + '.msh'
            output = self.hpc_options['output_files_location'] + self.name

        #config file for eigenmode smiulation which will be output as a Json file
        config = {
            "Problem":
            {
                "Type": "Electrostatic",
                "Verbose": 2,
                "Output": output
            },
            "Model":
            {
                "Mesh":  mesh_location,
                "L0": 1e-3,  
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
                        "Attributes": [physical_groups['air_box']],  # Air
                        "Permeability": 1.0,
                        "Permittivity": 1.0,
                        "LossTan": 0.0
                    },
                    {
                        "Attributes": [physical_groups['dielectric_substrate']],  # Dielectric
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
                    "Attributes": [physical_groups['far_field']]  # Metal trace
                },
                "Terminal": metal_terminals
            },
            "Solver":
            {
                "Order": self.user_options["solver_order"],
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

        return config


    def _process_metals_for_config_file(self):
        
        metal_terminals = []

        for m, metal in enumerate(self.metal_cap_physical_group):
            terminals = {}
            terminals['Index'] = m+1
            terminals['Attributes'] = [metal]
            metal_terminals.append(terminals)

        return metal_terminals