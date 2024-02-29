# PALACE
## Introduction
PALACE (Parallel Large Scale Computational Electromagnetics) is an open-source electromagnetic simulation platform developed by Amazon Web Services (AWS) for use on high-performance computing (HPC) architectures in order to improve simulation time and performance for simlulations relating to quantum computing hardware. Although there are many proprietary software platforms that support computational electromagnetic simulations (i.e. COMSOL, ANSYS) there exist very few open-source options for highly-parallel finite element based computational electromagnetics. Palace is a solution to this problem and offers a variety of simulations and features including:
- Electrostatic Simulations for Capacitance matrices
- Magnetostatic Simulations for Inductance matrices
- Eigenmode Simulations
- Driven Radio Frequency Simulations
- Device, Mesh and Field Visualisations through Paraview
Palace was originally developed for users to run simulations on AWS cloud-based HPC clusters, however we have been able to install and run Palace on `Bunya' the UQ HPC cluster.
## Installation
### Installation on UQ HPC Cluster (Bunya)
To log into Bunya you can use command prompt (for windows users) or terminal (for mac users) as both have in-built Secure Shell Protocol (SSH) functionality which allows secure remote access to Bunya along with command line execution. Upon opening command prompt/terminal on your PC enter the following to establish a ssh connection to Bunya:

```ssh <your_usernmame>@bunya.rcc.uq.edu.au```

Then follow the prompts, entering in you password and completing the multi-factor authentication (MFA) step to gain access to Bunya. Once logged-in to your account on Bunya, change directories to your scratch directory with the following command:

```cd /scratch/user/<your_username> ``` 

The scratch directory provides you with 150 GB of data and is the best place to store project files and data in the medium-term. In the long-term you should transfer files from your scratch directory to your PC using Secure File Transfer Protocol (SFTP) software such as winSCP as the administrators will routinely remove data that has been stored for too long on Bunya. Once in your scratch directory create a folder to house the Palace repository and then change into that directory:

```mkdir Palace-Project```

```cd Palace-Project```

We will be following the 'Build from Source' installation instructions detailed on the Palace Github documentation page located at https://awslabs.github.io/palace/stable/install/, so feel free to review this page, however I will provide all necessary commands in these install instructions.

Lets now open an interactive session in Bunya to build and install Palace. Run the following command to set-up an interactive session (note: this command should all be one line when you enter it into command prompt):

``` salloc --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --mem=50G --job-name=build_palace --time=01:00:00 --partition=debug--account=a_fedorov srun --export=PATH,TERM,HOME,LANG --pty /bin/bash -l ```

We need to load the following modules to complete the build and installation process, so run the following commands:

```module load foss/2021a```

```module load cmake/3.20.1-gcccore-10.3.0```

```module load pkgconfig/1.5.4-gcccore-10.3.0-python```

Now we can clone the Palace repository into our folder:

```git clone https://github.com/awslabs/palace.git --recurse-submodules```

Before we build Palace, we need to remove a few tests that will run and fail during the build process causing the installation to fail. These tests can be removed and it won't affect the functionality of Palace when it comes to running simulations. The first test we will remove is for an eigenmode solver called PetSc. Run the following commands from your Palace-Project folder to get to the directory where the file ExternalPETSc.cmake is stored:

```cd /palace/cmake/extern/```

We can open and then edit the 'ExternalPETSc.cmake' file with the following command:  

```nano ExternalPETSc.cmake```

Scroll down the page using the down arrow until you reach the following code block:

![Alt text](Petsc_before.png)

Here, we want to delete the line 'TEST_BEFORE_INSTALL' and change the option for 'TEST_COMMAND' to '' '', so the code block should now look like the following:

![Alt text](Petsc_after.png)

Save these changes to the file by holding 'Ctrl' and 'o', then hit 'Enter' to write the changes. You can then exit from the file by holding 'Ctrl' and 'x', this should take you back to the command prompt terminal.

The next test we want to remove is located in the file 'PkgConfigHelpers.cmake' and is located in the directory above our current working directory. If we run the following commands we can open this file:

```nano PkgConfigHelpers.cmake```

Scroll down to the following code block:

![Alt text](SLEPC_before.png)

Here we are removing the test code for SLEPc which is an eigenvalue/eigenvalue solver for large sparse matrices. Now, delete the 'try_run' block and the 'if else' block leaving only the following code:

![Alt text](SLEPC_after.png)

Again, save these changes to the file by holding 'Ctrl' and 'o', then hit 'Enter' to write the changes. Exit from the file by holding 'Ctrl' and 'x'.

We are now in a position to build Palace. Navigate to the directory where the repository was cloned, you can do this by ```cd ..``` and then opening the directory palace (lowercase) with ```cd palace```. Now, run the following commands to build Palace:

```mkdir build && cd build```

```cmake ..```

```make -j```

The build process generally takes between 10 - 15 minutes. After the build, you can check if the binary executable file has been created by navigating to the bin directory and listing the files:

```cd bin```

```ls```

A screenshot of running these commands and the output is shown below, the binary executable file is \url{palace-x86_64.bin}.

![Alt text](../bin_exec.png)