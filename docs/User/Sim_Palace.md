# Simulations using PALACE

Once the design has been completed using *Qiskit-Metal*, the design object can be used to run simulations in [AWS Palace](https://awslabs.github.io/palace/stable/). Currently, *SQDMetal* supports:

- Capacitance matrix simulations
- RF s-parameter simulations
- Eigenmode simulations

## Installation

Make sure to install the `gmsh` in the Python environment. Now to run simulations in Palace, the distribution must be installed either on a local computer or a cluster. We offer several implementations for a local PC:

- [Ubuntu](Sim_Palace_Ubuntu.md) (recommended)
- [Windows via WSL](WSL_Sim_palace.md)



### Singularity/Apptainer container

Instructions to build palace as a Singularity/Apptainer container are given [here](https://awslabs.github.io/palace/stable/install/#Build-using-Singularity/Apptainer). By default, this setup builds the latest commit of [awslabs/palace on main](https://github.com/awslabs/palace) with OpenMP support (no GPU). 

To use these containers for simulations, set the `palace-dir` user option to point to the `palace.sif` file (or whatever you named your container).

### HPC cluster

Working installation instructions for the Australian Bunya cluster are given [here](HPC_documentation.md).

## Usage (local PC)


```python
from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
from SQDMetal.Utilities.Materials import MaterialInterface

#Eigenmode Simulation Options
user_defined_options = {
                "mesh_refinement":  0,                             #refines mesh in PALACE - essetially divides every mesh element in half
                "dielectric_material": "silicon",                  #choose dielectric material - 'silicon' or 'sapphire'
                "starting_freq": 5e9,                              #starting frequency in Hz 
                "number_of_freqs": 1,                              #number of eigenmodes to find
                "solns_to_save": 1,                                #number of electromagnetic field visualizations to save
                "solver_order": 2,                                 #increasing solver order increases accuracy of simulation, but significantly increases sim time
                "solver_tol": 1.0e-8,                              #error residual tolerance foriterative solver
                "solver_maxits": 200,                              #number of solver iterations
                "fillet_resolution":12,                            #number of vertices per quarter turn on a filleted path
                "palace_dir":"~/spack/opt/spack/linux-ubuntu24.04-zen2/gcc-13.3.0/palace-develop-36rxmgzatchgymg5tcbfz3qrmkf4jnmj/bin/palace",#"PATH/TO/PALACE/BINARY",
                "num_cpus": 16                                     #number of cpus used in simulation
                }

#Create the Palace Eigenmode simulation
eigen_sim = PALACE_Eigenmode_Simulation(name ='x-mon_test',                                         #name of simulation
                                        metal_design = design,                                      #feed in qiskit metal design
                                        sim_parent_directory = "",                                  #choose directory where mesh file and config file are saved
                                        mode = 'simPC',                                             #choose simulation mode 'HPC' or 'simPC'                                  
                                        meshing = 'GMSH',                                           #choose meshing 'GMSH' or 'COMSOL'
                                        user_options = user_defined_options,                        #provide options chosen above
                                        create_files = True)                                        #create mesh, config and HPC batch files

#Add in metals from layer 1 of the design file
eigen_sim.add_metallic(1)

#Add in ground plane for simulation
eigen_sim.add_ground_plane()

#Add in the Josephson junction as a lumped port
eigen_sim.create_port_JosephsonJunction('junction', L_J=4.3e-9, C_J=10e-15)

#Fine-mesh x-mon
eigen_sim.fine_mesh_around_comp_boundaries(['x-mon'], min_size=8e-6, max_size=100e-6, taper_dist_min=10e-6, metals_only=False)

#Sets up the lossy interfaces for MA, SA and MS interfaces 
eigen_sim.setup_EPR_interfaces(metal_air=MaterialInterface('Aluminium-Vacuum'), substrate_air=MaterialInterface('Silicon-Vacuum'), substrate_metal=MaterialInterface('Silicon-Aluminium'))

#Prepares the mesh file and config file
eigen_sim.prepare_simulation()
```

Details on viewing the data in the end are shown [here](Sim_Palace_Data_viewer.md).

## Usage (HPC)

Here, some extra parameters need to be supplied. Specifically the path to a JSON file (supplied in the user option: `"HPC_Parameters_JSON"`) containing the parameters. Here is a template:

```json
{
    "HPC_nodes": 4,
    "sim_memory": "300G",
    "sim_time": "20:00:00",
    "account_name": "sam",
    "palace_location": "/scratch/project/palace-sqdlab/Palace-Project/palace/build/bin/palace-x86_64.bin",
    "input_dir": "/scratch/project/palace-sqdlab/inputFiles/"
}
```
