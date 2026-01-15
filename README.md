# SQDMetal

[![Documentation](https://img.shields.io/badge/Documentation-blue)](https://sqdlab.github.io/SQDMetal/)
[![UnitTests](https://github.com/sqdlab/SQDMetal/actions/workflows/Windows.yaml/badge.svg)](https://github.com/sqdlab/SQDMetal/actions/workflows/Windows.yaml)

Tools to aid in simulating and fabricating superconducting quantum devices. The tools are an extension of [~~Qiskit~~Quantum-Metal](https://github.com/Qiskit/qiskit-metal) to provide additional support for:

- Extra components with more flexible user-friendly options (see [gallery](https://nbviewer.org/github/sqdlab/SQDMetal/blob/main/docs/User/Comps_All.ipynb))
- Visualising and simulating effects of shadow evaporation techniques used to fabricate qubits
- RF and DC simulations using **COMSOL** (including calculation of capacitance matrices, inductance ~~matrices~~ and RF s-parameters)
- RF and DC simulations using cluster-friendly **AWS PALACE** (including **meshing via either COMSOL or Gmsh**)
- GDS export and manipulation techniques to help with fabrication setup for multi-die wafers, arrayed structures, and more

There are two classes of documentation provided for this stack:

- [User documentation](https://sqdlab.github.io/SQDMetal/)
- [Developer documentation](docs/Developer/Readme.md)

## Installation instructions

The following installation instructions automatically installs Qiskit-Metal along with SQDMetal. First choose a folder to house SQDMetal (idea is to create an editable folder such that the code can be modified and pushed without upsetting the pip package manager). Once navigating to this folder, run *Anaconda prompt* and run the following command:

```bash
cd C:/Users/....../myFolder/
git clone https://github.com/sqdlab/SQDMetal.git
```

Now run (changing `sqdmetal_env` to any other desired name for the virtual environment):

```bash
conda create -n sqdmetal_env python==3.12
```

Now activate the environment and install Qiskit-Metal:

```bash
activate sqdmetal_env
pip install -e SQDMetal
```

This should install Qiskit-Metal and SQDMetal.


