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

## Usage
