Installing AWS Palace in Singularity/Apptainer container
========================================================

Instructions to build palace as a Singularity/Apptainer container are given `here <https://awslabs.github.io/palace/stable/install/#Build-using-Singularity/Apptainer>`__. By default, this setup builds the latest commit of `awslabs/palace on main <https://github.com/awslabs/palace>`__ with *OpenMP* support (no GPU).

To use these containers for simulations, set the ``palace-dir`` user option to point to the ``palace.sif`` file (or whatever you named your container).

TODO: Properly verify this...

