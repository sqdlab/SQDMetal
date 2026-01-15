Installing AWS Palace on Ubuntu and Mac
=======================================

In terminal, run:

.. code-block:: bash

    sudo apt update
    sudo apt install file bzip2 ca-certificates g++ gcc gfortran git gzip lsb-release patch python3 tar unzip xz-utils zstd
    sudo apt install libxft-dev libglu1-mesa libxi-dev libxmu-dev libglu1-mesa-dev build-essential

Make a new directory where we clone all repositories:

.. code-block:: bash

    cd ~
    mkdir repo
    cd repo

Now clone the spack source and enter its environment here:

.. code-block:: bash

    git clone --depth=2 https://github.com/spack/spack.git
    . spack/share/spack/setup-env.sh

Update spack compilers:


.. code-block:: bash

    spack compiler find


Ensure that it says gcc is at least above version 11.

Clone palace and `install using spack <https://github.com/awslabs/palace/issues/360#issuecomment-2874057931>`__:

.. code-block:: bash

    git clone https://github.com/awslabs/palace.git

From here, we go through the Palace installation process for :ref:`Ubuntu` and :ref:`Mac` separately, before covering the :ref:`installation of Paraview and Pyvista <Recommended add-ons>`.

Ubuntu
------

Now the idea is that we create a local virtual environment in spack called `spack-env`, add local palace install as a repo in it and concretize+install. Note that you can skip the creation and activation of the environment, but you may run into issues due to dependecies if you have used spack previously. This is done via:

.. code-block:: bash

    source ./spack/share/spack/setup-env.sh
    spack repo add ./palace/spack_repo/local
    spack env create spack-env
    spack env activate spack-env
    spack add local.palace@develop

Now concretize+install via:

.. code-block:: bash

    spack concretize -f
    spack install --only-concrete
    find -name palace*

The second command will find the path of the Palace binary. Basically it's somewhere like (remember to add the `~/` for the full absolute path):

.. code-block:: bash

    ./spack/opt/spack/linux-ubuntu24.04-zen2/gcc-13.3.0/palace-develop-36rxmgzatchgymg5tcbfz3qrmkf4jnmj/bin/palace

Executing the above block should show the command-line switches required for Palace. 

Mac
---

Note that the exact process for Mac installation changes depending on your chip and operating system. The method documented here has been tested on a M5 MacBook with Tahoe 26.1 installed.

To install Palace, ope the terminal from the directory containing ``/spack``, run:

.. code-block:: bash

    source ./spack/share/spack/setup-env.sh
    spack install veclibfort
    spack install palace +mumps +slepc ^blas=veclibfort ^lapack=veclibfort

The first line activates ``spack`` in your shell. The second line installs a wrapper for the builtin MacOS Accelerator (which can replace BLAS and LAPACK). The third line initiates a full Palace installation that is suitable for arm64 MacBooks.

After the installation is complete (roughly 30 mins) the Palace binary can be found in ``spack/opt/spack/darwin-m1/palace-0.14.0-aijcrhkhxpjki7lbvx336dbzvwh52enp/bin/palace``.


Recommended add-ons
-------------------

Optionally, install `Paraview <https://www.paraview.org/>`__:

.. code-block:: bash

    sudo apt install paraview

If this fails (e.g. some dependency clash), use its `flatpak`:

.. code-block:: bash

    sudo apt install flatpak
    flatpak remote-add --if-not-exists flathub https://dl.flathub.org/repo/flathub.flatpakrepo
    sudo flatpak install flathub org.paraview.ParaView

Finally, in the Python virtual environment, if it has not been already installed, run:

.. code-block:: bash

    pip install pyvista

This is required to open the simulation field data via the API functions.


OLD NOTES LEFT IN CASE OF ISSUES - TO BE REMOVED IN FUTURE
----------------------------------------------------------

This method has been tested to work on Kubuntu. Ensure GIT has been installed via:

```bash
sudo apt-get install git
```

Now install the `Spack Package Manager <https://spack-tutorial.readthedocs.io/en/latest/tutorial_basics.html>`__ by running in console (in some directory of preference - e.g. somewhere in the home directory):

.. code-block:: bash

    cd ~/
    git clone -c feature.manyFiles=true --depth=2 https://github.com/spack/spack.git
    cd ~/spack
    . share/spack/setup-env.sh

On MacOS, you will have to install and setup `Sundials <https://github.com/LLNL/sundials>`__ before installing Palace. This is done by running:

.. code-block:: bash

    spack install sundials
    export SUNDIALS_DIR=$(spack location -i sundials)

Now run (`@develop` tag is required to get the latest version, i.e, 0.13 or newer):

.. code-block:: bash

    spack install palace@develop
    find -name palace*

The second command will find the path of the Palace binary. Basically it's somewhere like (remember to add the `~/` for the full absolute path):

.. code-block:: bash

    ./spack/opt/spack/linux-ubuntu24.04-zen2/gcc-13.3.0/palace-develop-36rxmgzatchgymg5tcbfz3qrmkf4jnmj/bin/palace

Executing the above block should show the command-line switches required for Palace. 

Optionally, install `Paraview <https://www.paraview.org/>`__:

.. code-block:: bash

    sudo apt install paraview

If this fails (e.g. some dependency clash), use its `flatpak`:

.. code-block:: bash

    sudo apt install flatpak
    flatpak remote-add --if-not-exists flathub https://dl.flathub.org/repo/flathub.flatpakrepo
    sudo flatpak install flathub org.paraview.ParaView

Finally, in the Python virtual environment, if it has not been already installed, run:

.. code-block:: bash

    pip install pyvista

This is required to open the simulation field data via the API functions.

**Troubleshooting**

Errors can appear after updating the Linux system. Simply reinstall palace as above (starting from the ``setup-env.sh`` command). However, first run:

.. code-block:: bash

    sudo rm -r .spack/

Alternatively, `edit the compilers.yaml file <https://stackoverflow.com/questions/67899951/change-version-of-gcc-which-does-not-support-compiling-c-programs-using-the-co>`__ to the latest GCC version.

If OpenMPI processes fail to spawn (i.e. Palace does not start and/or throws the error ``local rank failed   --> Returned value Not found (-13) instead of ORTE_SUCCESS``), try running the following `command in terminal <https://askubuntu.com/questions/730/how-do-i-set-environment-variables>`__`:

.. code-block:: bash

    export PMIX_MCA_gds=hash

In the case of a Jupyter Notebook (i.e. running the commands in *SQDMetal*), add this to the initialisation cell:


.. code-block:: python

    import os
    os.environ["PMIX_MCA_gds"]="hash"

If it throws an error for `libCEED` (see `here <https://github.com/awslabs/palace/issues/257>`__ for details) or an error like `SUNDIALS: Core: *** NOT FOUND ***`, first open the file:

.. code-block:: bash

    nano ~/spack/var/spack/repos/builtin/packages/palace/package.py

Edit the line which says `depends_on("libxsmm@main")` to `depends_on("libxsmm@=main")` (note the `=` sign). Then add `depends_on("sundials")` after the line `depends_on("eigen")`. Now proceed with the installation by running the palace installation command once more.

**Upgrading**

Delete all traces of `spack` via:

.. code-block:: bash

    sudo rm -r .spack/
    sudo rm -r spack/

Then start from the beginning to reinstall to the latest version of `spack` and subsequently the latest version of `palace`.

**Installing Palace v0.13**
As of 18 Feb, 2024, Palace v0.13 is not available by default to install in `spack`. To install Palace v0.13.0, you need to edit the `spack` palace installation information contained in `package.py`:

.. code-block:: bash

    spack edit palace

Add the following line (where the Palace versions are listed) to add v0.13.0 to available installs:

.. code-block:: bash

    version("0.13.0", tag="v0.13.0", commit="a61c8cbe0cacf496cde3c62e93085fae0d6299ac")

Save the package.py file for the Palace install, and then run

.. code-block:: bash

    spack install palace@0.13
