from RF_Simulation import RF_Simulation
from Utilities.Materials import Material

class Eigenmode_Simulation(RF_Simulation):

    def __init__(self, name, ports_dict, user_options, jj_dict, hpc_options = {}):
        super().__init__(name, ports_dict, user_options, jj_dict, hpc_options)

    
    def create_sim_config_file(self, physical_groups):
        
        #get the material propoerties for the dielectric specified by the user
        dielectric = Material(self.user_options['dielectric_material'])

        #prepare ports - inherited method from parent class RF_Simulation
        config_ports = self._process_ports_for_config_file()

        #if hpc_options is not empty then supply correct location for mesh file and where to output simulation results
        if self.hpc_options:
            mesh_location = self.hpc_options['input_files_location'] + self.name + '/' + self.name + '.msh'
            output = self.hpc_options['output_files_location'] + self.name
        else:
            mesh_location = self.user_options['sim_directory'] + self.name + '.msh'
            output = self.user_options['sim_directory'] + self.name

        if physical_groups.get('dielectric_gaps') is None:
            physical_groups['dielectric_gaps'] = []

        #config file for eigenmode smiulation which will be output as a Json file
        config = {
            "Problem":
            {
                "Type": "Eigenmode",
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
                    "Attributes": [physical_groups['metals'], physical_groups['far_field']]  # Metal trace
                },
                "LumpedPort": config_ports,
                "Postprocessing":
                {
                   "Dielectric":
                    [
                        {
                        "Index": 1,
                        "Attributes": [physical_groups['dielectric_gaps']],
                        "Type": "SA",
                        "Thickness": 2.0e-6,  
                        "Permittivity": 10.0,
                        "LossTan": 0.9e-3
                        },
                        {
                        "Index": 2,
                        "Attributes": [physical_groups['metals']],
                        "Type": "MS",
                        "Thickness": 2.0e-6,  
                        "Permittivity": 10.0,
                        "LossTan": 0.3e-3
                        },
                        {
                        "Index": 3,
                        "Attributes": [physical_groups['metals']],
                        "Type": "MA",
                        "Thickness": 2.0e-6, 
                        "Permittivity": 10.0,
                        "LossTan": 0.5e-3
                        }
                    ] 
                }
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

        return config