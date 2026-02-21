Simulations using PALACE
========================

Once the design has been completed using *Qiskit-Metal*, the design object can be used to run simulations in `AWS Palace <https://awslabs.github.io/palace/stable/>`__. Currently, *SQDMetal* supports:

- Capacitance matrix simulations
- RF s-parameter simulations
- Eigenmode simulations

The raw data can be viewed via *ParaView* or our :doc:`in-built visualisation class <simpalacedv>`.

Installation
------------

Make sure to install the `gmsh` in the Python environment. Now to run simulations in Palace, the distribution must be installed either on a local computer or a cluster. We offer several implementations for a local PC:

.. toctree::
   :maxdepth: 1
   :caption: Helper utility modules

   Ubuntu/Mac <instpalaceunix>
   Windows via WSL <isntpalacewsl>
   Apptainer container <instpalaceapp>

Usage (local PC)
----------------

.. note::
    We recommend browsing our :doc:`worked examples <../examples/index>` for further details beyond the example usage shown below.

.. code-block:: python

    from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
    from SQDMetal.Utilities.Materials import MaterialInterface

    #Eigenmode Simulation Options
    user_defined_options = {
                    "mesh_refinement":  0,                             #refines mesh in PALACE - essetially divides every mesh element in half
                    "dielectric_material": "silicon",                  #choose dielectric material - 'silicon' or 'sapphire'
                    "starting_freq": 5e9,                              #starting frequency in Hz 
                    "number_of_freqs": 1,                              #number of eigenmodes to find
                    "solns_to_save": 1,                                #number of electromagnetic field visualizations to save
                    "solver_order": 2,                                 #increasing solver order increases accuracy of simulation, but significantly increases sim time
                    "solver_tol": 1.0e-8,                              #error residual tolerance foriterative solver
                    "solver_maxits": 200,                              #number of solver iterations
                    "fillet_resolution":12,                            #number of vertices per quarter turn on a filleted path
                    "palace_dir":"~/spack/opt/spack/linux-ubuntu24.04-zen2/gcc-13.3.0/palace-develop-36rxmgzatchgymg5tcbfz3qrmkf4jnmj/bin/palace",#"PATH/TO/PALACE/BINARY",
                    "num_cpus": 16                                     #number of cpus used in simulation
                    }

    #Create the Palace Eigenmode simulation
    eigen_sim = PALACE_Eigenmode_Simulation(name ='x-mon_test',                                         #name of simulation
                                            metal_design = design,                                      #feed in qiskit metal design
                                            sim_parent_directory = "",                                  #choose directory where mesh file and config file are saved
                                            mode = 'simPC',                                             #choose simulation mode 'HPC' or 'simPC'                                  
                                            meshing = 'GMSH',                                           #choose meshing 'GMSH' or 'COMSOL'
                                            user_options = user_defined_options,                        #provide options chosen above
                                            create_files = True)                                        #create mesh, config and HPC batch files

    #Add in metals from layer 1 of the design file
    eigen_sim.add_metallic(1)

    #Add in ground plane for simulation
    eigen_sim.add_ground_plane()

    #Add in the Josephson junction as a lumped port
    eigen_sim.create_port_JosephsonJunction('junction', L_J=4.3e-9, C_J=10e-15)

    #Fine-mesh x-mon
    eigen_sim.fine_mesh_components(['x-mon'], min_size=8e-6, max_size=100e-6, taper_dist_min=10e-6, metals_only=False)

    #Sets up the lossy interfaces for MA, SA and MS interfaces 
    eigen_sim.setup_EPR_interfaces(metal_air=MaterialInterface('Aluminium-Vacuum'), substrate_air=MaterialInterface('Silicon-Vacuum'), substrate_metal=MaterialInterface('Silicon-Aluminium'))

    #Prepares the mesh file and config file
    eigen_sim.prepare_simulation()

Details on viewing the data in the end are shown :doc:`here <simpalacedv>`.

Usage (HPC)
-----------

Here, some extra parameters need to be supplied. Specifically the path to a JSON file (supplied in the user option: `"HPC_Parameters_JSON"`) containing the parameters. Here is a template:

.. code-block:: json

    {
        "HPC_nodes": 4,
        "sim_memory": "300G",
        "sim_time": "20:00:00",
        "account_name": "sam",
        "palace_location": "/scratch/project/palace-sqdlab/Palace-Project/palace/build/bin/palace-x86_64.bin",
        "input_dir": "/scratch/project/palace-sqdlab/inputFiles/"
    }

Functions
---------

General functions
~~~~~~~~~~~~~~~~~

- `add_metallic`
- `add_ground_plane`

- `enforce_full_3D_simulation`

- `set_xBoundary_as_proportion` (note it uses that percentage on both sides. So 0.1 on a 1mm chip means 100um on negative side and 100um on positive side)
- `set_yBoundary_as_proportion`
- `set_zBoundary_as_proportion`
- `set_xBoundary_as_absolute`
- `set_yBoundary_as_absolute`
- `set_zBoundary_as_absolute`


- `fine_mesh_features`
- `fine_mesh_along_path` (QM)
- `fine_mesh_in_rectangle`
- `fine_mesh_components` (QM)

- `enable_mesh_refinement`

- `prepare_simulation`
- `set_local_output_subdir`
- `run`

- `retrieve_data`
- `retrieve_simulation_sizes`
- `retrieve_simulation_sizes_from_file` (Static)
- `open_mesh_gmsh`

RF (Eigenmode and Frequency-Driven)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `create_port_2_conds` (QM)
- `create_port_JosephsonJunction` (QM)
- `create_port_CPW_on_Launcher` (QM)
- `create_port_CPW_on_Route` (QM)
- `create_port_CPW_via_edge_point`

- `create_CPW_feed_Uclip_on_Launcher` (QM)
- `create_CPW_feed_Uclip_on_Route` (QM)

- `set_port_impedance`
- `create_waveport_on_boundary`

- `setup_EPR_interfaces`
- `add_kinetic_inductance`

- `set_farfield`

Eigenmode
~~~~~~~~~

- `set_freq_search`
- `retrieve_field_plots`
- `retrieve_EPR_data`
- `retrieve_EPR_data_from_file`
- `retrieve_mode_port_EPR`
- `retrieve_mode_port_EPR_from_file`

- `calculate_hamiltonian_parameters_EPR`
- `calculate_hamiltonian_parameters_EPR_from_files`

Frequency-Driven
~~~~~~~~~~~~~~~~

- `set_port_excitation`
- `set_freq_values`

- `add_surface_current_source_region`

- `retrieve_data_from_file`

- `get_waveport_modes`
- `get_waveport_modes_from_file`


Capacitance Matrix
~~~~~~~~~~~~~~~~~~

- `display_conductor_indices`

- `calc_params_floating_Transmon`
- `calc_params_floating_Transmon_from_files`

