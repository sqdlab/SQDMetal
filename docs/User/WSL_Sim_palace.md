# Prepare UBUNTU

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
