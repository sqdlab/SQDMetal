# Simulations using PALACE

Once the design has been completed using *Qiskit-Metal*, the design object can be used to run simulations in [AWS Palace](https://awslabs.github.io/palace/stable/). Currently, *SQDMetal* supports:

- Capacitance matrix simulations
- RF s-parameter simulations

## Installation

Make sure to install the `gmsh` in the Python environment. Now to run simulations in Palace, the distribution must be installed either on a local computer or a cluster.

### Local PC

This method has been tested to work on Kubuntu. Ensure GIT has been installed via:

```
sudo apt-get install git
```

Now install the  [Spack Package Manager](https://spack-tutorial.readthedocs.io/en/latest/tutorial_basics.html) by running in console (in some directory of preference - e.g. somewhere in the home directory):

```
git clone --depth=100 --branch=releases/v0.21 https://github.com/spack/spack.git ~/spack
cd ~/spack
. share/spack/setup-env.sh
```

Now run:

```
spack install palace
find -name palace*
```

The second command will find the path of the Palace binary. Basically it's somewhere like:

```
./spack/opt/spack/linux-ubuntu23.04-zen2/gcc-12.3.0/palace-0.11.2-v3z2x65yy2n6p7ryhwnzzxhroauejven/bin/palace
```

Executing the above block should show the command-line switches required for Palace. Finally, install [Paraview](https://www.paraview.org/):

```
sudo apt install paraview
```

### HPC cluster

Working installation instructions for the Australian Bunya cluster are given [here](HPC_documentation.md).

## Usage


```python
from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
from SQDMetal.PALACE.SQDGmshRenderer import Palace_Gmsh_Renderer

#Eigenmode Simulation Options
user_defined_options = {
                 "mesh_refinement":  0,                             #refines mesh in PALACE - essetially divides every mesh element in half
                 "dielectric_material": "silicon",                  #choose dielectric material - 'silicon' or 'sapphire'
                 "starting_freq": 7.5,                              #starting frequency in GHz 
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
