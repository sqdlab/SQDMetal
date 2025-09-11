# Simulations using PALACE

Once the design has been completed using *Qiskit-Metal*, the design object can be used to run simulations in [AWS Palace](https://awslabs.github.io/palace/stable/). Currently, *SQDMetal* supports:

- Capacitance matrix simulations
- RF s-parameter simulations

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

#Eigenmode Simulation Options
user_defined_options = {
                 "mesh_refinement":  0,                             #refines mesh in PALACE - essetially divides every mesh element in half
                 "dielectric_material": "silicon",                  #choose dielectric material - 'silicon' or 'sapphire'
                 "starting_freq": 7.5e9,                            #starting frequency in Hz 
                 "number_of_freqs": 4,                              #number of eigenmodes to find
                 "solns_to_save": 4,                                #number of electromagnetic field visualizations to save
                 "solver_order": 2,                                 #increasing solver order increases accuracy of simulation, but significantly increases sim time
                 "solver_tol": 1.0e-8,                              #error residual tolerance foriterative solver
                 "solver_maxits": 100,                              #number of solver iterations
                 "comsol_meshing": "Extremely fine",                #level of COMSOL meshing: 'Extremely fine', 'Extra fine', 'Finer', 'Fine', 'Normal'
                 "mesh_max": 120e-3,                                #maxiumum element size for the mesh in mm
                 "mesh_min": 10e-3,                                 #minimum element size for the mesh in mm
                 "mesh_sampling": 130,                              #number of points to mesh along a geometry
                 "sim_memory": '300G',                              #amount of memory for each HPC node i.e. 4 nodes x 300 GB = 1.2 TB
                 "sim_time": '20:00:00',                            #allocated time for simulation 
                 "HPC_nodes": '4',                                  #number of Bunya nodes. By default 20 cpus per node are selected, then total cores = 20 x HPC_nodes
                 "fillet_resolution":12                             #Number of vertices per quarter turn on a filleted path
                }

#Creat the Palace Eigenmode simulation
eigen_sim = PALACE_Eigenmode_Simulation(name ='single_resonator_example_eigen',                     #name of simulation
                                        metal_design = design,                                      #feed in qiskit metal design
                                        sim_parent_directory = "",            #choose directory where mesh file, config file and HPC batch file will be saved
                                        mode = 'HPC',                                               #choose simulation mode 'HPC' or 'simPC'                                          
                                        meshing = 'GMSH',                                           #choose meshing 'GMSH' or 'COMSOL'
                                        user_options = user_defined_options,                        #provide options chosen above
                                        view_design_gmsh_gui = False,                               #view design in GMSH gui 
                                        create_files = True)                                        #create mesh, config and HPC batch files
eigen_sim.add_metallic(1)
eigen_sim.add_ground_plane()

#Add in the RF ports
eigen_sim.create_port_CPW_on_Launcher('LP1', 20e-3)
eigen_sim.create_port_CPW_on_Launcher('LP2', 20e-3)
#Fine-mesh routed paths
eigen_sim.fine_mesh_along_path(100e-6, 'resonator1', mesh_sampling=130, mesh_min=5e-3, mesh_max=120e-3)
eigen_sim.fine_mesh_along_path(100e-6, 'TL', mesh_sampling=130, mesh_min=7e-3, mesh_max=120e-3)
#Fine-mesh a rectangular region
eigen_sim.fine_mesh_in_rectangle(-0.14e-3, -1.33e-3, 0.14e-3, -1.56e-3, mesh_sampling=130, mesh_min=5e-3, mesh_max=120e-3)

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
