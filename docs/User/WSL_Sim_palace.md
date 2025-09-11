# Install guide for Palace on WSL


# Install WSL

ref: https://learn.microsoft.com/en-us/windows/wsl/install

**NOTE: you mst be on windows10 version 2004 or Windows 11**

>check your windows version in settings, then system section, scroll down to  "about".

**Open powershell in admin mode**

>press start, type powershell, right click on it and select "Run as administrator".

```powershell
wsl --set-default-version 2
wsl --install Ubuntu-22.04
```


## Install spack basic on a new WSL system:

ref : https://spack.readthedocs.io/en/latest/installing_prerequisites.html#verify-spack-prerequisites

```bash
sudo apt update
sudo apt install file bzip2 ca-certificates g++ gcc gfortran git gzip lsb-release patch python3 tar unzip xz-utils zstd
```

## Install some more ubuntu build stuff:

```bash
sudo apt install libxft-dev libglu1-mesa libxi-dev libxmu-dev libglu1-mesa-dev build-essential
```

**CLOSE THE SHELL NOW AND REOPEN**.


# Prepare DIR
## Make a new directory where we clone all repos:

```bash
mkdir repo
cd repo
```


# Prepare SPACK

## Clone and source Spack here:

```bash
git clone --depth=2 https://github.com/spack/spack.git
. spack/share/spack/setup-env.sh
```

## Update spack compilers:

```bash
spack compiler find
```

# Clone palace and install using spack:

ref: https://github.com/awslabs/palace/issues/360#issuecomment-2874057931

## Clone palace :

and checkout cameron's branch, which has fixes that lets spack correctly install palace0.13.

```bash
git clone https://github.com/awslabs/palace.git

# These changes are now in master but need to be tested again.
#git checkout cameronrutherford/0.13.0-spack-cmake-constraints
```

## More steps:

* Create local env of spack called `spack-env`.
* add local palace install as a repo in it.
* concretize and install.

```bash
source ./spack/share/spack/setup-env.sh
spack env create -d ./spack-env
spack env activate ./spack-env
spack repo add ./palace/spack/local
spack add local.palace@develop
```

now concretize and install error free!!!

```bash
spack concretize -f
spack install --only-concrete
```

## No need for EXTERNAL MFEM, libCEED or sundial !! palace ships with its own !!!!
