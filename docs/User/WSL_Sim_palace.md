# Install guide for Palace on WSL


## Install WSL

**NOTE: you must be on Windows 10 version 2004 or Windows 11**

>check your windows version in settings, then system section, scroll down to  "about".

Using the instructins from [here](https://learn.microsoft.com/en-us/windows/wsl/install), first **Open powershell in admin mode** (i.e. press start, type `powershell`, right click on it and select "Run as administrator") and type:

```powershell
wsl --set-default-version 2
wsl --install Ubuntu-22.04
```

If this fails (e.g. `WslRegisterDistribution failed`), then run `wsl --update --web-download` first. Restart windows if it has pending updates. Note the username/password supplied during installation.

From now onwards, when we state *enter WSL*, we mean:
- Enter powershell (administrative mode)
- Type `wsl` and hit `ENTER`


## Install spack basic on a new WSL system:

Enter WSL and verify that the nameservers are correct by first running:

```bash
sudo nano /etc/resolv.conf
```

Ensure it says `nameserver 8.8.8.8` (if not, edit it, hit `CTRL-X` and `Y`).

Now install spack as outlined [here](https://spack.readthedocs.io/en/latest/installing_prerequisites.html#verify-spack-prerequisites) by typing:

```bash
sudo apt update
sudo apt install file bzip2 ca-certificates g++ gcc gfortran git gzip lsb-release patch python3 tar unzip xz-utils zstd
sudo apt install libxft-dev libglu1-mesa libxi-dev libxmu-dev libglu1-mesa-dev build-essential
```

**CLOSE THE POWERSHELL NOW AND REOPEN**.

Make a new directory where we clone all repositories:

```bash
wsl
cd ~
mkdir repo
cd repo
```

Note that we use `~/repo` as the default palace repository directory. Since we will be using virtual environments if it is in a different directory, then make sure to set `'palace_wsl_spack_repo_directory'` in the `user_options` dictionary in the palace simulation object later.

Now clone the spack source and enter its environment here:

```bash
git clone --depth=2 https://github.com/spack/spack.git
. spack/share/spack/setup-env.sh
```

Update spack compilers:

```bash
spack compiler find
```

## Install Palace

Clone palace and [install using spack](https://github.com/awslabs/palace/issues/360#issuecomment-2874057931):

```bash
git clone https://github.com/awslabs/palace.git

#Checkout cameron's branch, which has fixes that lets spack correctly install palace0.13.
#But these changes are now in master but need to be tested again.
#git checkout cameronrutherford/0.13.0-spack-cmake-constraints
```

Now the idea is that we create a local virtual environment in spack called `spack-env`, add local palace install as a repo in it and concretize+install. This is done via:

```bash
source ./spack/share/spack/setup-env.sh
spack env create -d ./spack-env
spack env activate ./spack-env
spack repo add ./palace/spack_repo/local
spack add local.palace@develop
```

Now concretize+install via:

```bash
spack concretize -f
spack install --only-concrete
```

There is no need for EXTERNAL MFEM, libCEED or sundial as palace ships with its own versions. Now run:

```bash
find -name palace*
```

The command will find the path of the Palace binary. Basically it's somewhere like:

```bash
./spack/opt/spack/linux-ubuntu24.04-zen2/gcc-13.3.0/palace-develop-36rxmgzatchgymg5tcbfz3qrmkf4jnmj/bin/palace
```

Executing the above block should show the command-line switches required for Palace. 

## When running the simulations

When running Palace simulations, the simulation object receives a dictionary via the argument `user_options`. Make sure to set:

- `'palace_mode': 'wsl'`,
- `'palace_wsl_spack_repo_directory': '~/repo'`
- `'palace_dir':"~/repo/spack/opt/spack/linux-zen2/palace-develop-3ofp7n4fjqj5i6slvei3w6nptzdiwdma/bin/palace"`

where the second key points to a different directory if the repository directory was chosen to be different during the installation procedure. Note that the final key must be the absolute path. Thus, one adds `~/repo` to the path found using `find -name palace*` earlier.

