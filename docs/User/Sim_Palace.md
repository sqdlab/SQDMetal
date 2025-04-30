# Simulations using PALACE

Once the design has been completed using *Qiskit-Metal*, the design object can be used to run simulations in [AWS Palace](https://awslabs.github.io/palace/stable/). Currently, *SQDMetal* supports:

- Capacitance matrix simulations
- RF s-parameter simulations

## Installation

Make sure to install the `gmsh` in the Python environment. Now to run simulations in Palace, the distribution must be installed either on a local computer or a cluster.

### Local PC

This method has been tested to work on Kubuntu. Ensure GIT has been installed via:

```bash
sudo apt-get install git
```

Now install the  [Spack Package Manager](https://spack-tutorial.readthedocs.io/en/latest/tutorial_basics.html) by running in console (in some directory of preference - e.g. somewhere in the home directory):

```bash
cd ~/
git clone -c feature.manyFiles=true --depth=2 https://github.com/spack/spack.git
cd ~/spack
. share/spack/setup-env.sh
```

On MacOS, you will have to install and setup [Sundials](https://github.com/LLNL/sundials) before installing Palace. This is done by running:

```bash
spack install sundials
export SUNDIALS_DIR=$(spack location -i sundials)
```

Now run (`@develop` tag is required to get the latest version, i.e, 0.13 or newer):

```bash
spack install palace@develop
find -name palace*
```

The second command will find the path of the Palace binary. Basically it's somewhere like (remember to add the `~/` for the full absolute path):

```bash
./spack/opt/spack/linux-ubuntu24.04-zen2/gcc-13.3.0/palace-develop-36rxmgzatchgymg5tcbfz3qrmkf4jnmj/bin/palace
```

Executing the above block should show the command-line switches required for Palace. 

Finally, install [Paraview](https://www.paraview.org/):

```bash
sudo apt install paraview
```

Finally, in the Python virtual environment, if it has not been already installed, run:

```bash
pip install pyvista
```

This is required to open the simulation field data via the API functions.

**Troubleshooting**

Errors can appear after updating the Linux system. Simply reinstall palace as above (starting from the `setup-env.sh` command). However, first run:
```bash
sudo rm -r .spack/
```
Alternatively, [edit the compilers.yaml file](https://stackoverflow.com/questions/67899951/change-version-of-gcc-which-does-not-support-compiling-c-programs-using-the-co) to the latest GCC version.

If OpenMPI processes fail to spawn (i.e. Palace does not start and/or throws the error `local rank failed   --> Returned value Not found (-13) instead of ORTE_SUCCESS`), try running the following [command in terminal](https://askubuntu.com/questions/730/how-do-i-set-environment-variables):
```bash
export PMIX_MCA_gds=hash
```

In the case of a Jupyter Notebook (i.e. running the commands in *SQDMetal*), add this to the initialisation cell:
```python
import os
os.environ["PMIX_MCA_gds"]="hash"
```

If it throws an error for `libCEED` (see [here](https://github.com/awslabs/palace/issues/257) for details) or an error like `SUNDIALS: Core: *** NOT FOUND ***`, first open the file:
```bash
nano ~/spack/var/spack/repos/builtin/packages/palace/package.py
```
Edit the line which says `depends_on("libxsmm@main")` to `depends_on("libxsmm@=main")` (note the `=` sign). Then add `depends_on("sundials")` after the line `depends_on("eigen")`. Now proceed with the installation by running the palace installation command once more.

**Upgrading**

Delete all traces of `spack` via:

```bash
sudo rm -r .spack/
sudo rm -r spack/
```
Then start from the beginning to reinstall to the latest version of `spack` and subsequently the latest version of `palace`.

**Installing Palace v0.13**
As of 18 Feb, 2024, Palace v0.13 is not available by default to install in `spack`. To install Palace v0.13.0, you need to edit the `spack` palace installation information contained in `package.py`:

```bash
spack edit palace
```

Add the following line (where the Palace versions are listed) to add v0.13.0 to available installs:

```bash
version("0.13.0", tag="v0.13.0", commit="a61c8cbe0cacf496cde3c62e93085fae0d6299ac")
```

Save the package.py file for the Palace install, and then run

```bash
spack install palace@0.13
```

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
