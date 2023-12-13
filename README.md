# SQDMetal

Tools to aid in simulating and fabricating superconducting quantum devices. The tools are an extension of [Qiskit-Metal](https://github.com/Qiskit/qiskit-metal) to provide additional support for:

- COMSOL modelling (including calculation of capacitance matrices, ~~inductance matrices~~ and RF s-parameters)
- Visualising and simulating effects of shadow evaporation techniques used to fabricate qubits

There are two classes of documentation provided for this stack:

- [User documentation](docs/User/Readme.md)
- [Developer documentation](docs/Developer/Readme.md)

## Installation instructions

The following installation instructions automatically installs Qiskit-Metal along with SQDMetal. First run Anaconda prompt and run the following command to create an environment (in this example, the name is sqdmetal_env):

```
conda create -n sqdmetal_env
```

Now activate the environment and install Qiskit-Metal:

```
activate sqdmetal_env
```

Ideally, the pip installation of qiskit-metal should work, but there are some dependencies that are better manually installed first before proceeding. So run the following commands (hitting `y` when required to proceed with the installation, this will take some time):

```
conda install -c conda-forge gdspy
conda install -c conda-forge pyside2
conda install -c conda-forge gmsh python-gmsh
```

The last one is optional and only required if using Palace. Now choose a folder to house SQDMetal (idea is to create an editable folder such that the code can be modified and pushed without upsetting the pip package manager). Once navigating to this folder, run the usual GIT clone:

```
cd C:/Users/....../myFolder/
git clone https://github.com/sqdlab/SQDMetal.git
```

Now there is a SQDMetal folder in the current directory. Do not enter this new folder. Simply run:

```
pip install -e SQDMetal
```

This should install Qiskit-Metal and SQDMetal.


