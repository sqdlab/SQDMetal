# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0

from pathlib import Path
from SQDMetal.Utilities.QUtilities import QUtilities
from SQDMetal.Utilities.ShapelyEx import ShapelyEx
from SQDMetal.Utilities.Materials import MaterialInterface, MaterialConductor
from SQDMetal.COMSOL.Model import COMSOL_Model
from SQDMetal.COMSOL.SimRFsParameter import COMSOL_Simulation_RFsParameters
from SQDMetal.PALACE.Utilities.GMSH_Geometry_Builder import GMSH_Geometry_Builder
from SQDMetal.PALACE.Utilities.GMSH_Mesh_Builder import GMSH_Mesh_Builder
from SQDMetal.PALACE.Utilities.GMSH_Navigator import GMSH_Navigator

from SQDMetal.Utilities.GeometryProcessors.GeomQiskitMetal import GeomQiskitMetal
from SQDMetal.Utilities.GeometryProcessors.GeomGDS import GeomGDS

import numpy as np
import os
import subprocess
import sys
import gmsh
import json
import platform
import shapely
import shutil

class PALACE_Model:
    #TODO: Refactor with suffix "Base"

    default_user_options_parent = {
        "dielectric_material": "silicon",
        "fillet_resolution": 12,
        "mesh_min": 10e-6,
        "mesh_max": 100e-6,
        "taper_dist_min": 30e-6,
        "taper_dist_max": 200e-6,
        "threshold": 1e-9,
        "fuse_threshold": 1e-9,
        "simplify_edge_min_angle_deg": -1,
        "gmsh_verbosity": 1,
        "gmsh_dist_func_discretisation": 120,
        'palace_mode': 'local',
        'palace_dir': '',
        'palace_wsl_spack_repo_directory': '~/repo',
        "comsol_meshing": "Extremely fine"
    }

    def __init__(self, meshing:str, mode:str, options: dict, **kwargs):
        """
        Base class for PALACE simulation classes.

        Parameters
        ----------
        meshing : str 
            The meshing engine to use. It can be either 'GMSH' or 'COMSOL' with the latter
            requiring a local COMSOL installation.
        mode : str 
            Either 'PC' to run locally or 'HPC' to just generate the scripts to be run on
            a cluster/computer later.
        options : dict
            One must supply **one** of these keyword arguments depending on the input format of
            the design to simulate:

            *   ``'metal_design'`` (QDesign):
                    The Qiskit-Metal design object to simulate
            *   ``'gds_design'`` (`str`):
                    Path to the GDS file with the design to simulate

            Additional options which include those of the daughter classes and:

            *   ``'dielectric_material'`` (`str`):
                    Name of the dielectric to use on the substrate
            *   ``'fillet_resolution'`` (`float`):
                    All curved sections of paths (e.g. a CPW line) will be discretised into line
                    segments. This value is the number of such segments observed over a 90° arc.
            *   ``'mesh_min'`` (`float`):
                    The default minimum mesh element size to use when using the fine-meshing functions.
            *   ``'mesh_max'`` (`float`):
                    The default maximum mesh element size to use when using the fine-meshing functions.
            *   ``'taper_dist_min'`` (`float`):
                    The default distance to which the size of the mesh elements remain at mesh_min.
            *   ``'taper_dist_max'`` (`float`):
                    The default distance to which the size of the mesh elements begin to be set at mesh_max.
            *   ``'threshold'`` (`float`):
                    The default threshold (in metres) to which elements in a given layer are simplified. That
                    is, adjacent vertices that are closer than this value are simplified into a single vertex.
            *   ``'fuse_threshold'`` (`float`):
                    The default threshold (in metres) to which cracks are bridged together. Gaps that are twice
                    this value are effectively bridged together so that there are no seams etc.
            *   ``'simplify_edge_min_angle_deg'`` (`float`):
                    The default simplification used to delete unnecessary vertices when they are adjacent
                    and collinear. It weights them based on the distance from the central vertex. For
                    example, consider three vertices ABC. If AB and BC are equidistant and the angle ABC
                    is less than the supplied value, then B is removed as it is considered collinear. If
                    the length of AB is greater than BC, then the ratio of the lengths R is multiplied by
                    the angle ABC before checking for this threshold. That is, if there is a genuine sharp
                    feature, it is not simplified. This is mostly designed for curved edges that may have
                    an unnecessary number of vertices.
            *   ``'gmsh_verbosity'`` (`int`):
                    The default level of messages to print from GMSH. The levels symbolised by the integers
                    are: Nothing (0), Errors (1), Warnings (2), Direct (3), Information (4), Status (5),
                    Debug (99). The verbosity of the messages increase with the listed levels.
            *   ``'gmsh_dist_func_discretisation'`` (`int`):
                    When a path or polygon is given for a fine-mesh field, GMSH needs to compute distances
                    from these 1D or 2D regions in order to compute the mesh element size at a given location.
                    Thus, it must be able to numerically approximate these lines/polygons. It does so by chopping
                    every line-segment or polygon into smaller portions given by this supplied value. Note
                    that making this too large will make the generation of the mesh slower as it will take
                    longer to compute all the mesh field distances.
            *   ``'palace_mode'`` (`str`):
                    When running Palace locally, one may run it via a 'local' installation (e.g. Mac-OS,
                    Linux etc.) or 'wsl' (e.g. Windows).
            *   ``'palace_dir'`` (`str`):
                    When running Palace locally, the location of the Palace binary must be supplied so that
                    SQDMetal can call it directly.
            *   ``'palace_wsl_spack_repo_directory'`` (`str`):
                    When running Palace in WSL (i.e. Windows), the location of the cloned repository (inside
                    the WSL VM) needs to be supplied unless it is in the default '~/repo' folder.
            *   ``'comsol_meshing'`` (`str`):
                    The default meshing size to use if using COMSOL to mesh the simulation. Acceptable values
                    are: Can be: 'Extremely fine', 'Extra fine', 'Finer', 'Fine', 'Normal', 'Coarse', 'Coarser',
                    'Extra coarse' or 'Extremely coarse'.

        Returns
        -------
            The constructed PALACE_Model object.
        """

        #TODO: Look to refactor/remove mode?
        #TODO: Add link to WSL installation details when mentioning `palace_wsl_spack_repo_directory`?

        self.meshing = meshing
        self._metallic_layers = []
        self._ground_plane = {'omit': True}
        self._fine_meshes = []
        self._sim_config = ""
        self._input_dir = ""
        self._output_subdir = ""
        self._EPR_setup = False
        self._use_KI = False
        self._KI = 0
        #
        self._mesh_refinement = {"UniformLevels": 0}

        self._process_geometry_type(**kwargs)
        self._boundary_distances = {}
        self.set_xBoundary_as_proportion(0.1)
        self.set_yBoundary_as_proportion(0.1)
        self.set_zBoundary_as_proportion(1.0,0.0)

        for key in PALACE_Model.default_user_options_parent:
            self.user_options[key] = options.get(key, PALACE_Model.default_user_options_parent[key])

        if mode == 'HPC':
            with open(options["HPC_Parameters_JSON"], "r") as f:
                self.hpc_options = json.loads(f.read())
        else:
            self.palace_dir = options.get('palace_dir', 'palace')
            self.palace_mode = options.get('palace_mode', 'local')  #Can be local, wsl, docker
            self._palace_wsl_repo_location = options.get('palace_wsl_spack_repo_directory', '~/repo')
            self._palace_wsl_spack_env_name = options.get('palace_wsl_spack_env_name', 'spack-env')
            self.hpc_options = {"input_dir":""}
        self._num_cpus = options.get('num_cpus', 16)
        self._num_threads = options.get('num_threads', 1)
        self._full_3D_params = {'metal_thickness':0, 'substrate_trenching':0}


    def _process_geometry_type(self, **kwargs):
        if 'metal_design' in kwargs:
            self._geom_processor = GeomQiskitMetal(kwargs['metal_design'], **kwargs)
        elif 'gds_design' in kwargs:
            self._geom_processor = GeomGDS(kwargs['gds_design'], **kwargs)
    
    
    def create_batch_file(self):
        
        #note: I have disabled naming the output file by setting '# SBATCH' instead of '#SBATCH' 
        #so I can get the slurm job number to use for testing

        sbatch = {
                "header": "#!/bin/bash --login",
                "job_name": "#SBATCH --job-name=" + self.name,
                "output_loc": "# SBATCH --output=" + self.name + ".out",
                "error_out": "#SBATCH --error=" + self.name + ".err",
                "partition": "#SBATCH --partition=general",
                "nodes": "#SBATCH --nodes=" + self.hpc_options["HPC_nodes"],
                "tasks": "#SBATCH --ntasks-per-node=20",
                "cpus": "#SBATCH --cpus-per-task=1",
                "memory": "#SBATCH --mem=" + self.hpc_options["sim_memory"],
                "time": "#SBATCH --time=" + self.hpc_options['sim_time'],
                "account": "#SBATCH --account=" + self.hpc_options['account_name'],
                "foss": "module load foss/2021a",
                "cmake": "module load cmake/3.20.1-gcccore-10.3.0",
                "pkgconfig": "module load pkgconfig/1.5.4-gcccore-10.3.0-python",
                "run_command": f"srun {self.hpc_options['palace_location']} " + self.hpc_options["input_dir"] + self.name + "/" + self.name + ".json"
        }
    
        #check simulation mode and return appropriate parent directory 
        parent_simulation_dir = self._check_simulation_mode()

        #create sbatch file name
        sim_file_name = self.name + '.sbatch'

        #destination for config file
        simulation_dir = parent_simulation_dir + str(self.name)

        #save to created directory
        file = os.path.join(simulation_dir, sim_file_name)

        #write sbatch dictionary to file
        with open(file, "w+", newline = '\n') as f:
            for value in sbatch.values():
                f.write('{}\n'.format(value))

    def _create_config_file(self):
        pass
    def _prepare_simulation(self, metallic_layers, ground_plane):
        raise NotImplementedError()

    def prepare_simulation(self):
        """
        Creates the JSON and mesh files as relevant to the defined Palace simulation. The model will
        have everything required to run the simulation either externally via command line or internally
        using :func:`~SQDMetal.PALACE.Model.PALACE_Model.run`.
        """
        self._prepare_simulation(self._metallic_layers, self._ground_plane)

    def add_metallic(self, layer_id, **kwargs):
        """
        Adds metallic conductors from the design input onto the surface layer of the chip simulation.
        If the particular layer has fancy PVD evaporation steps, the added metallic layer will account
        for said steps and merge the final result. In addition, all metallic contiguous elements are
        merged into single blobs.

        This function must be run before calling :func:`~SQDMetal.PALACE.Model.PALACE_Model.prepare_simulation`.

        Parameters
        ----------
        layer_id : int
            The layer ID containing the metallic layer to add.
        kwargs : dict
            Optional arguments:

            *   ``'threshold'`` (`float`):
                    Defaults to -1 (i.e. no threshold). This is the threshold in metres, below which
                    consecutive vertices along a given polygon are combined into a single vertex. This
                    simplification helps with meshing as the engine will not overdo the meshing. If
                    this argument is negative, the argument is ignored.
            *   ``'fuse_threshold'`` (`float`):
                    Defaults to 1e-12. This is the minimum distance between metallic elements, below
                    which they are considered to be a single polygon and thus, the polygons are merged
                    with the gap filled. This accounts for floating-point errors that make adjacent
                    elements fail to merge as a single element, due to infinitesimal gaps between them.
            *   ``'evap_mode'`` (`str`):
                    Defaults to 'separate_delete_below'. These are the methods upon which to separate
                    or merge overlapping elements across multiple evaporation steps. See documentation
                    on PVD_Shadows for more details on the available options.
            *   ``'group_by_evaporations'`` (`bool`):
                    Defaults to False. If set to True, if elements on a particular evaporation step are
                    separated due to the given evap_mode, they will still be selected as a part of the
                    same conductor (useful for example, in capacitance matrix simulations).
            *   ``'evap_trim'`` (`float`):
                    Defaults to 20e-9. This is the trimming distance used in certain evap_mode profiles.
                    See documentation on PVD_Shadows for more details on its definition.
            
            If using GDS inputs:
            
            *   ``'cell_index'`` (`int`):
                    Defaults to 0. This specifies the cell from which the layer is to be extracted.
        """
        #TODO: Link PVD_Shadows to documentation in above docstring
        new_metallic_layer = {
                'type': 'design_layer',
                'layer_id': layer_id,
                'threshold': kwargs.pop('threshold', -1),
                'fuse_threshold': kwargs.pop('fuse_threshold', 1e-12),
                'evap_mode': kwargs.pop('evap_mode', 'separate_delete_below'),
                'group_by_evaporations': kwargs.pop('group_by_evaporations', False),
                'evap_trim': kwargs.pop('evap_trim', 20e-9),
                'other' : kwargs
            }
        self._metallic_layers.append(new_metallic_layer)

    def add_ground_plane(self, **kwargs):
        """
        *Applies Qiskit-Metal designs only*
        
        Adds metallic ground-plane from the Qiskit-Metal design object onto the surface layer of the
        chip simulation.

        This function must be run before calling :func:`~SQDMetal.PALACE.Model.PALACE_Model.prepare_simulation`.

        Parameters
        ----------
        kwargs : dict
            Optional argument:

            *   ``'threshold'`` (`float`):
                    Defaults to -1 (i.e. no threshold). This is the threshold in metres, below which
                    consecutive vertices along a given polygon are combined into a single vertex. This
                    simplification helps with meshing as the engine will not overdo the meshing. If
                    this argument is negative, the argument is ignored.
        """
        self._ground_plane = {'omit':False, 'threshold':kwargs.get('threshold', -1)}


    def _is_native_arm64(self):
        '''
        Helper method written by Sadman Ahmed Shanto @shanto268 in this issue: https://github.com/sqdlab/SQDMetal/issues/9
        '''
        try:
            # Run the sysctl command to check the native architecture
            result = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True, check=True)
            return "M" in result.stdout  # M is for Apple M1/M2 chips
        except subprocess.CalledProcessError:
            # print(f"Error checking native architecture: {e}")
            return False

    def _run_local(self):
        """
        Runs the PALACE simulation locally.
        """

        config_file = self._sim_config
        leFile = os.path.basename(os.path.realpath(config_file))
        leDir = os.path.dirname(os.path.realpath(config_file))

        if not os.path.exists(self._output_data_dir):
            os.makedirs(self._output_data_dir)
        self.log_location = f"{self._output_data_dir}/out.log"
        with open("temp.sh", "w+") as f:
            f.write(f"cd \"{leDir}\"\n")
            f.write(f"{self.palace_dir} -np {self._num_cpus} -nt {self._num_threads} {leFile} | tee \"{self.log_location}\"\n")

        # Set execute permission on temp.sh
        os.chmod("temp.sh", 0o755)

        #Just create the output log file in case the terminal shell requires it to exist in order to write into it...
        with open(self.log_location, 'w'):
            pass

        # Get the current running architecture
        running_arch = platform.machine()

        # If the machine is native arm64 but running under x86_64, run temp.sh under arm64
        if self._is_native_arm64() and running_arch == "x86_64":
            # print("Running arm64 process from x86_64 environment...")
            self.cur_process = subprocess.Popen(["arch", "-arm64", "bash", "./temp.sh"], shell=False)
        else:
            # Run the temp.sh script as usual
            # print("Running temp.sh script under current architecture...")
            self.cur_process = subprocess.Popen("./temp.sh", shell=True)
        
        try:
            self.cur_process.wait()
        except KeyboardInterrupt:
            self.cur_process.kill()
        self.cur_process = None

    def _run_local_wsl(self):
        """
        Runs the PALACE simulation locally.
        """

        config_file = self._sim_config
        leFile = os.path.basename(os.path.realpath(config_file))
        leDir = os.path.dirname(os.path.realpath(config_file))
        leCWD = os.getcwd()

        def conv_winpath_to_WSLmount(lePath):
            lePathWSL = lePath.replace('\\','/').replace(':','')  #Turn \ into / and remove the : after the drive letter...
            lePathWSL = lePathWSL[0].lower() + lePathWSL[1:]        #Turn the C or D into c or d...
            return '/mnt/' + lePathWSL

        leDirWSL = conv_winpath_to_WSLmount(leDir)
        leCWDWSL = conv_winpath_to_WSLmount(leCWD)

        if not os.path.exists(self._output_data_dir):
            os.makedirs(self._output_data_dir)
        self.log_location = f"{self._output_data_dir}/out.log"
        log_locationWSL = conv_winpath_to_WSLmount(self.log_location)

        with open("temp.sh", "w+", newline='\n') as f:
            f.write(f"cd {self._palace_wsl_repo_location}\n")
            f.write(f"source ./spack/share/spack/setup-env.sh\n")
            f.write(f"spack env activate ./{self._palace_wsl_spack_env_name}\n")
            f.write(f"cd \"{leDirWSL}\"\n")
            f.write(f"{self.palace_dir} -np {self._num_cpus} -nt {self._num_threads} {leFile} | tee \"{log_locationWSL}\"\n")

        # Set execute permission on temp.sh
        os.chmod("temp.sh", 0o755)

        self.cur_process = subprocess.Popen(f"wsl -e bash -ic \"{leCWDWSL}/temp.sh\"",
                                                stdout=subprocess.PIPE,
                                                stderr=subprocess.STDOUT,
                                                bufsize=1,
                                                text=True
                                            )
        while True:
            output_line = self.cur_process.stdout.readline()
            if output_line == '' and self.cur_process.poll() is not None:
                break  # Process finished and no more output
            if output_line:
                print(output_line.strip()) # Print and remove trailing newline
                sys.stdout.flush() # Ensure immediate display in the notebook
        
        try:
            self.cur_process.wait()
        except KeyboardInterrupt:
            self.cur_process.kill()
        self.cur_process = None

    def _run_local_container(self, container_executable_name: str = "apptainer", **kwargs):
        """
        Runs the PALACE simulation using a container.

        This method executes the PALACE simulation inside an Apptainer or Singularity
        container image specified by `self.palace_dir`. It uses the simulation
        configuration file prepared by `prepare_simulation` and the number of
        CPUs specified by `self._num_cpus`. It also sets the number of threads per CPU
        by 'self._num_threads'.

        Palace installation instructions for Apptainer can be found here: https://awslabs.github.io/palace/dev/install/#Build-using-Singularity/Apptainer

        Args:
            container_executable_name (str): The name of the container executable
                                             to use ('apptainer' or 'singularity').
                                             Default is 'apptainer'.

        Raises:
            AssertionError: If `prepare_simulation` has not been run.
            ValueError: If `container_executable_name` is not 'apptainer' or 'singularity'.
            AssertionError: If the container image specified by `self.palace_dir` does not exist.
            Exception: If there is an error resolving the container path.
        """

        #TODO: Make container_executable_name a global option if even required...abs
        #Then delete all kwargs in run

        if container_executable_name not in ["apptainer", "singularity"]:
            raise ValueError(
                f"Invalid container executable name: {container_executable_name}. Must be 'apptainer' or 'singularity'."
            )

        container_executable_path = shutil.which(container_executable_name)
        if container_executable_path is None:
            raise FileNotFoundError(
                f"Container executable {container_executable_name} not found in PATH."
            )

        if "~" in self.palace_dir:
            container_path = Path(self.palace_dir).expanduser()
        else:
            container_path = Path(self.palace_dir)

        config_file = Path(self._sim_config).resolve()
        output_data_dir = Path(self._output_data_dir).resolve()
        if not output_data_dir.exists():
            output_data_dir.mkdir(parents=True)

        # Prepare the command to run inside the container
        # Assumes 'palace' executable is in the container's PATH
        # The config file (config_file_name) will be accessible inside the container
        # because subprocess.run's cwd is set to config_file_dir, which Apptainer mounts by default.

        # Prepare the full command to execute using apptainer
        # Use the resolved container path and pass the config file name as an argument
        full_command = [
            container_executable_path,
            "run",
            str(container_path),
            "-np",
            str(self._num_cpus),
            "-nt",
            str(self._num_threads),
            config_file.name,
        ]

        # Open log file and run subprocess using nested context managers
        log_output = []  # List to buffer output

        try:
            # Execute the command using subprocess.Popen for real-time output
            # text=True decodes stdout/stderr as text
            # cwd sets the working directory
            with subprocess.Popen(
                full_command,
                cwd=config_file.parent,
                stdout=subprocess.PIPE,
                text=True,  # Decode stdout/stderr as text
            ) as process:
                # Read output in real-time and print to console and buffer
                for line in iter(process.stdout.readline, ""):
                    print(line, end="")
                    log_output.append(line)  # Append to list

                # The context manager waits for the process to terminate
                # and checks the return code (raises CalledProcessError for non-zero exit)
                # Any remaining output will have been read by the iter loops above

            log_location = output_data_dir / "out.log"
            with open(log_location, "w") as log_file:
                log_file.writelines(log_output)
            print(f"Output saved to {log_location}")

            shutil.copy(self._sim_config, self._output_data_dir + "/config.json")

            return self.retrieve_data()

        except Exception as e:
            print(f"Error: {e}")
            raise

    def run(self, **kwargs):
        """
        Runs the PALACE simulation locally either via a local installation, WSL or a container VM depending on the configuration.
        """
        #TODO: Remove the kwargs - the Apptainer should be specified by the user options anyway...
        assert self._sim_config != "", "Must run prepare_simulation at least once."

        if self.palace_mode == 'local':
            self._run_local(**kwargs)
        elif self.palace_mode == 'wsl':
            self._run_local_wsl()
        elif self.palace_dir.endswith(".sif"):  # If palace is run using a container    #TODO: Need to test this properly...
            self._run_local_container(**kwargs)

        shutil.copy(self._sim_config, self._output_data_dir + '/config.json')
        return self.retrieve_data()

    def retrieve_data(self):
        """
        Retrieves relevant data after running the current simulation. Note that this is only relevant when using
        :func:`~SQDMetal.PALACE.Model.PALACE_Model.run` to run the simulation locally. In addition, the daughter class
        must override this function to actually retrieve and return data relevant to the simulation.

        Returns
        -------
            Returns data relevant to the simulation; if this function is not overridden, nothing if returned.
        """
        pass
    

    def _check_simulation_mode(self):
        """method to check the type of simualtion being run and return
        the parent directory to store the simulation files"""

        parent_simulation_dir = None

        # Check the simulation mode is valid.
        if self.mode in ["HPC", "simPC", "PC"]:
            parent_simulation_dir = self.sim_parent_directory
        else:
            raise ValueError("Invalid simulation mode entered.")

        return parent_simulation_dir



    def _create_directory(self, directory_name):
        '''create a directory to hold the simulation files'''

        parent_simulation_dir = self._check_simulation_mode()

        # Directory
        directory = directory_name
  
        # Path
        path = os.path.join(parent_simulation_dir, directory)
  
        # Create the directory
        if not os.path.exists(path):
            os.makedirs(path)
            print("Directory '% s' created" % directory)


    def _save_mesh_gmsh(self):
        '''function used to save the gmsh mesh file'''
        if self.create_files:
            parent_simulation_dir = self._check_simulation_mode()

            # file_name
            file_name = self.name + "/" + self.name + ".msh"
    
            # Path
            path = os.path.join(parent_simulation_dir, file_name)
            if os.path.exists(path):
                os.remove(path)
            gmsh.write(path)
            self._mesh_path = path
            GMSH_Navigator(self._mesh_path).export_to_png()
    
    def open_mesh_gmsh(self):
        '''
        Note that if this is used in a Notebook environment, it'll be a blocking function until the Gmsh GUI is closed...
        '''
        gmsh_nav = GMSH_Navigator(self.path_mesh)
        gmsh_nav.open_GUI()

    def _save_mesh_comsol(self, comsol_obj):
        '''function used to save the comsol mesh file'''
        if self.create_files:
            parent_simulation_dir = self._check_simulation_mode()

            # file_name
            comsol_obj.save(parent_simulation_dir + self.name + "/" + self.name)
    
            # Path
            path = os.path.abspath(parent_simulation_dir) + "/" + self.name + "/" + self.name + ".mphbin"

            #COMSOL export commands
            comsol_obj._model.java.component("comp1").mesh("mesh1").export().set("filename", path)
            comsol_obj._model.java.component("comp1").mesh("mesh1").export(path)

    def _get_folder_prefix(self):
        #TODO: Consider removing?
        return self.hpc_options["input_dir"]  + self.name + "/" if self.hpc_options["input_dir"] != "" else ""

    def set_local_output_subdir(self, name):
        """
        Sets a subdirectory within the simulation directory so that when :func:`~SQDMetal.PALACE.Model.PALACE_Model.run`
        is called, the simulation output is placed within this subdirectory. For example, this is useful when running the
        same model, but varying the frequency range through which it is swept (i.e. no changes to geometry, mesh etc.).
        This function sets the subdirectory by modifying the simulation JSON file that is used by the Palace simulations
        when calling :func:`~SQDMetal.PALACE.Model.PALACE_Model.run`

        Parameters
        ----------
        name : str
            Name of the subdirectory.
        """
        self._output_subdir = str(name)
        self._output_dir = self._get_folder_prefix()  + "outputFiles"
        if self._output_subdir != "":
            self._output_dir += "/" + self._output_subdir
        if self._sim_config != "":
            with open(self._sim_config, "r") as f:
                config_json = json.loads(f.read())
            config_json['Problem']['Output'] = self._output_dir
            with open(self._sim_config, "w") as f:
                json.dump(config_json, f, indent=2)
            self._output_data_dir = os.path.dirname(os.path.realpath(self._sim_config)) + "/" + self._output_dir

    def _check_if_QiskitMetalDesign(self):
        assert isinstance(self._geom_processor, GeomQiskitMetal), "This function can only be used on QiskitMetal designs."

    def fine_mesh_features(self, max_feature_size:float, **kwargs):
        """
        Automatically detect the finer features on the device and select them to be meshed finely. Fine features
        are defined as fingers/slivers that have a particular feature width. For example, the feature size of a
        CPW line would be its conductor width.

        This function must be run before calling :func:`~SQDMetal.PALACE.Model.PALACE_Model.prepare_simulation`.

        Parameters
        ----------
        max_feature_size : float
            Maximum feature size in metres; any feature below this size will be selected automatically for fine
            meshing.
        kwargs : dict
            There are additional options that may be supplied:

            *   ``'mesh_min'`` (`float`):
                    The default minimum mesh element size to use when using the fine-meshing functions. If it
                    is not supplied, the default value from the user options will be taken.
            *   ``'mesh_max'`` (`float`):
                    The default maximum mesh element size to use when using the fine-meshing functions. If it
                    is not supplied, the default value from the user options will be taken.
            *   ``'taper_dist_min'`` (`float`):
                    The default distance to which the size of the mesh elements remain at mesh_min. If it is
                    not supplied, the default value from the user options will be taken.
            *   ``'taper_dist_max'`` (`float`):
                    The default distance to which the size of the mesh elements begin to be set at mesh_max. If
                    it is not supplied, the default value from the user options will be taken.

        Returns
        -------
            Returns nothing.
        """
        polys = self._geom_processor.get_polys_of_fine_features(max_feature_size, self._metallic_layers, self._ground_plane)
        self._fine_meshes.append({
            'type': 'arb_polys',
            'polys': polys,
            'min_size': kwargs.get('min_size', self.user_options['mesh_min']),
            'max_size': kwargs.get('max_size', self.user_options['mesh_max']),
            'taper_dist_min': kwargs.get('taper_dist_min', self.user_options['taper_dist_min']),
            'taper_dist_max': kwargs.get('taper_dist_max', self.user_options['taper_dist_max'])
        })

    def fine_mesh_along_path(self, dist_resolution:float, qObjName:str, trace_name:str='', **kwargs):
        """
        *Applies Qiskit-Metal designs only*

        Mesh finely around the areas surrounding a Qiskit-Metal path object. This is useful with features like
        a CPW.

        This function must be run before calling :func:`~SQDMetal.PALACE.Model.PALACE_Model.prepare_simulation`.

        Parameters
        ----------
        dist_resolution : float
            Distance in metres to which the curve is numerically approximated before passing into Gmsh. That is,
            the curves will be approximated by this segmentation distance.
        qObjName : str
            The name of the Qiskit-Metal QComponent object that contains said path.
        trace_name : str
            The trace name within the QComponent object. If left at the default blank string '', the path will
            be automatically inferred.
        kwargs : dict
            There are additional options that may be supplied:

            *   ``'mesh_min'`` (`float`):
                    The default minimum mesh element size to use when using the fine-meshing functions. If it
                    is not supplied, the default value from the user options will be taken.
            *   ``'mesh_max'`` (`float`):
                    The default maximum mesh element size to use when using the fine-meshing functions. If it
                    is not supplied, the default value from the user options will be taken.
            *   ``'taper_dist_min'`` (`float`):
                    The default distance to which the size of the mesh elements remain at mesh_min. If it is
                    not supplied, the default value from the user options will be taken.
            *   ``'taper_dist_max'`` (`float`):
                    The default distance to which the size of the mesh elements begin to be set at mesh_max. If
                    it is not supplied, the default value from the user options will be taken.

        Returns
        -------
            Returns nothing.
        """
        #TODO: Look into whether dist_resolution should be inferred? Perhaps from the filletting? Make sure to check/test with GMSH GUI
        self._check_if_QiskitMetalDesign()
        leUnits = QUtilities.get_units(self._geom_processor.design)
        lePath = QUtilities.calc_points_on_path(dist_resolution/leUnits, self._geom_processor.design, qObjName, trace_name)[0] * leUnits
        self._fine_meshes.append({
            'type': 'path',
            'path': lePath,
            'min_size': kwargs.get('min_size', self.user_options['mesh_min']),
            'max_size': kwargs.get('max_size', self.user_options['mesh_max']),
            'taper_dist_min': kwargs.get('taper_dist_min', self.user_options['taper_dist_min']),
            'taper_dist_max': kwargs.get('taper_dist_max', self.user_options['taper_dist_max'])
        })

    def fine_mesh_in_rectangle(self, x1:float, y1:float, x2:float, y2:float, **kwargs):
        """
        Fine mesh a rectangular region.

        This function must be run before calling :func:`~SQDMetal.PALACE.Model.PALACE_Model.prepare_simulation`.

        Parameters
        ----------
        x1 : float
            Minimum x-boundary of the rectangle (in metres).
        x2 : float
            Maximum x-boundary of the rectangle (in metres).
        y1 : float
            Minimum y-boundary of the rectangle (in metres).
        y2 : float
            Maximum y-boundary of the rectangle (in metres).
        kwargs : dict
            There are additional options that may be supplied:

            *   ``'mesh_min'`` (`float`):
                    The default minimum mesh element size to use when using the fine-meshing functions. If it
                    is not supplied, the default value from the user options will be taken.
            *   ``'mesh_max'`` (`float`):
                    The default maximum mesh element size to use when using the fine-meshing functions. If it
                    is not supplied, the default value from the user options will be taken.
            *   ``'taper_dist_min'`` (`float`):
                    The default distance to which the size of the mesh elements remain at mesh_min. If it is
                    not supplied, the default value from the user options will be taken.
            *   ``'taper_dist_max'`` (`float`):
                    The default distance to which the size of the mesh elements begin to be set at mesh_max. If
                    it is not supplied, the default value from the user options will be taken.

        Returns
        -------
            Returns nothing.
        """
        self._fine_meshes.append({
            'type': 'box',
            'x_bnds': (x1, x2),
            'y_bnds': (y1, y2),
            'min_size': kwargs.get('min_size', self.user_options['mesh_min']),
            'max_size': kwargs.get('max_size', self.user_options['mesh_max']),
            'taper_dist_min': kwargs.get('taper_dist_min', self.user_options['taper_dist_min']),
            'taper_dist_max': kwargs.get('taper_dist_max', self.user_options['taper_dist_max'])
        })

    def fine_mesh_around_comp_boundaries(self, list_comp_names, **kwargs):
        #TODO: Deprecate in March 2026
        print("fine_mesh_around_comp_boundaries has been refactored! It's called fine_mesh_components now; please make changes before it gets deprecated")
        self.fine_mesh_components(list_comp_names, **kwargs)

    def fine_mesh_components(self, list_comp_names:list[str], **kwargs):
        """
        *Applies Qiskit-Metal designs only*

        Mesh finely within the boundaries of a given list of components.

        This function must be run before calling :func:`~SQDMetal.PALACE.Model.PALACE_Model.prepare_simulation`.

        Parameters
        ----------
        list_comp_names : list[str]
            Names of the components (over which the regions are to be finely meshed) given as a list of strings.
        kwargs : dict
            There are additional options that may be supplied:

            *   ``'metals_only'`` (`bool`):
                    If set to True (defaults to False), only the metallic regions are selected for fine-meshing.
                    Otherwise, the entire components (i.e. both dielectric gaps as well as the metallic regions)
                    are selected.
            *   ``'mesh_min'`` (`float`):
                    The default minimum mesh element size to use when using the fine-meshing functions. If it
                    is not supplied, the default value from the user options will be taken.
            *   ``'mesh_max'`` (`float`):
                    The default maximum mesh element size to use when using the fine-meshing functions. If it
                    is not supplied, the default value from the user options will be taken.
            *   ``'taper_dist_min'`` (`float`):
                    The default distance to which the size of the mesh elements remain at mesh_min. If it is
                    not supplied, the default value from the user options will be taken.
            *   ``'taper_dist_max'`` (`float`):
                    The default distance to which the size of the mesh elements begin to be set at mesh_max. If
                    it is not supplied, the default value from the user options will be taken.

        Returns
        -------
            Returns nothing.
        """
        self._check_if_QiskitMetalDesign()
        self._fine_meshes.append({
            'type': 'comp_bounds',
            'list_comp_names': list_comp_names,
            'min_size': kwargs.get('min_size', self.user_options['mesh_min']),
            'max_size': kwargs.get('max_size', self.user_options['mesh_max']),
            'taper_dist_min': kwargs.get('taper_dist_min', self.user_options['taper_dist_min']),
            'taper_dist_max': kwargs.get('taper_dist_max', self.user_options['taper_dist_max']),
            'metals_only': kwargs.get('metals_only', False)
        })

    def retrieve_simulation_sizes(self):
        """
        Retrieves simulation sizes after running the current simulation. Note that this is only relevant when using
        :func:`~SQDMetal.PALACE.Model.PALACE_Model.run` to run the simulation locally.

        Returns
        -------
        sizes : dict
            The simulation sizes are given via the following keys:

            *   ``'DoF'`` (`int`):
                    Degrees of Freedom used in this simulation run.
            *   ``'MeshElements'`` (`int`):
                    Number of mesh elements used in this simulation run.
        """
        return self.retrieve_simulation_sizes_from_file(self._output_data_dir + '/palace.json')

    def enforce_full_3D_simulation(self, metal_thickness:float, substrate_trenching:float=0.0):
        """
        Enable a full 3D simulation where the metals have a vertical thickness and dielectric gaps can be trenched.
        If this function is not called, then the simulations will be quasi-3D in which the metals are deemed to
        be infinitesimally thin and thus, approximated as planar sheets on the surface of the 3D substrate.

        This function must be run before calling :func:`~SQDMetal.PALACE.Model.PALACE_Model.prepare_simulation`.

        Parameters
        ----------
        metal_thickness : float
            Thickness of the metals in metres. Can be zero.
        substrate_trenching : float
            Depth of the trenches (of the dielectric gaps relative to the bottom of the metals) in metres.
            Default value is 0.0.

        Returns
        -------
            Returns nothing.
        """
        self._full_3D_params['metal_thickness'] = metal_thickness
        self._full_3D_params['substrate_trenching'] = substrate_trenching

    @staticmethod
    def retrieve_simulation_sizes_from_file(path_palace_json:str):
        """
        Retrieves simulation sizes by scanning the stored output of a previou run simulation.

        Parameters
        ----------
        path_palace_json : str
            File path (either relative or absolute) to the *palace.json* file generated when completing a simulation.

        Returns
        -------
        sizes : dict
            The simulation sizes are given via the following keys:

            *   ``'DoF'`` (`int`):
                    Degrees of Freedom used in this simulation run.
            *   ``'MeshElements'`` (`int`):
                    Number of mesh elements used in this simulation run.
        """
        #Returns dictionary of DoF and Mesh size...
        with open(path_palace_json, "r") as f:
            config_json = json.loads(f.read())
        return {
            'DoF': config_json['Problem']['DegreesOfFreedom'],
            'MeshElements': config_json['Problem']['MeshElements']
        }

    def set_xBoundary_as_proportion(self, x_prop:float):
        """
        Set the distance of the chip to the air-box on both its boundaries along the x-axis. The distance is
        given as a proportion of the chip size along the x-axis.

        This function must be run before calling :func:`~SQDMetal.PALACE.Model.PALACE_Model.prepare_simulation`.

        Parameters
        ----------
        x_prop : float
            Proportion of the chip size along the x-axis to use when separating the chip from the air-box on both
            sides along the x-axis.
        """
        assert x_prop >= 0, "Boundary distance proportion along the x-axis must be a non-negative number."
        self._boundary_distances.pop('x_pos',None)
        self._boundary_distances.pop('x_neg',None)
        self._boundary_distances['x_prop'] = x_prop

    def set_yBoundary_as_proportion(self, y_prop:float):
        """
        Set the distance of the chip to the air-box on both its boundaries along the y-axis. The distance is
        given as a proportion of the chip size along the y-axis.

        This function must be run before calling :func:`~SQDMetal.PALACE.Model.PALACE_Model.prepare_simulation`.

        Parameters
        ----------
        y_prop : float
            Proportion of the chip size along the y-axis to use when separating the chip from the air-box on both
            sides along the y-axis.
        """
        assert y_prop >= 0, "Boundary distance proportion along the y-axis must be a non-negative number."
        self._boundary_distances.pop('y_pos',None)
        self._boundary_distances.pop('y_neg',None)
        self._boundary_distances['y_prop'] = y_prop

    def set_zBoundary_as_proportion(self, z_prop_top:float, z_prop_bottom:float):
        """
        Set the distance of the chip to the air-box on both its boundaries along the z-axis. The distances are
        given as a proportion of the chip size along the z-axis.

        This function must be run before calling :func:`~SQDMetal.PALACE.Model.PALACE_Model.prepare_simulation`.

        Parameters
        ----------
        z_prop_top : float
            Proportion of the chip thickness (i.e. along z-axis) to use when separating the top of the chip to the
            top ceiling of the air-box.
        z_prop_bottom : float
            Proportion of the chip thickness (i.e. along z-axis) to use when separating the bottom of the chip to the
            bottom floor of the air-box.
        """
        assert z_prop_top >= 0, "Boundary distance proportion above the chip must be a non-negative number."
        assert z_prop_bottom >= 0, "Boundary distance proportion below the chip must be a non-negative number."
        self._boundary_distances.pop('z_pos',None)
        self._boundary_distances.pop('z_neg',None)
        self._boundary_distances['z_prop_top'] = z_prop_top
        self._boundary_distances['z_prop_bottom'] = z_prop_bottom

    def set_xBoundary_as_absolute(self, x_neg, x_pos):
        """
        Set the distances of the chip to the air-box on both its boundaries along the x-axis.

        This function must be run before calling :func:`~SQDMetal.PALACE.Model.PALACE_Model.prepare_simulation`.

        Parameters
        ----------
        x_neg : float
            Distance of the chip to the air-box (in metres) on the negative x-axis side.
        x_pos : float
            Distance of the chip to the air-box (in metres) on the positive x-axis side.
        """
        assert x_neg >= 0, "Boundary on the negative x-axis of the chip must be a non-negative number."
        assert x_pos >= 0, "Boundary on the positive x-axis of the chip must be a non-negative number."
        self._boundary_distances.pop('x_prop',None)
        self._boundary_distances['x_neg'] = x_neg
        self._boundary_distances['x_pos'] = x_pos

    def set_yBoundary_as_absolute(self, y_neg, y_pos):
        """
        Set the distances of the chip to the air-box on both its boundaries along the y-axis.

        This function must be run before calling :func:`~SQDMetal.PALACE.Model.PALACE_Model.prepare_simulation`.

        Parameters
        ----------
        y_neg : float
            Distance of the chip to the air-box (in metres) on the negative y-axis side.
        y_pos : float
            Distance of the chip to the air-box (in metres) on the positive y-axis side.
        """
        assert y_neg >= 0, "Boundary on the negative y-axis of the chip must be a non-negative number."
        assert y_pos >= 0, "Boundary on the positive y-axis of the chip must be a non-negative number."
        self._boundary_distances.pop('y_prop',None)
        self._boundary_distances['y_neg'] = y_neg
        self._boundary_distances['y_pos'] = y_pos

    def set_zBoundary_as_absolute(self, z_neg, z_pos):
        """
        Set the distances of the chip to the air-box on both its boundaries along the z-axis.

        This function must be run before calling :func:`~SQDMetal.PALACE.Model.PALACE_Model.prepare_simulation`.

        Parameters
        ----------
        z_neg : float
            Distance of the chip to the air-box (in metres) on the negative z-axis side.
        z_pos : float
            Distance of the chip to the air-box (in metres) on the positive z-axis side.
        """
        assert z_neg >= 0, "Boundary on the negative z-axis of the chip must be a non-negative number."
        assert z_pos >= 0, "Boundary on the positive z-axis of the chip must be a non-negative number."
        self._boundary_distances.pop('z_prop_top',None)
        self._boundary_distances.pop('z_prop_bottom',None)
        self._boundary_distances['z_neg'] = z_neg
        self._boundary_distances['z_pos'] = z_pos

    def enable_mesh_refinement(self, num_iterations:int, max_DoFs:int=0, tolerance:float=1e-2, Dorfler_marking_fraction:float = 0.7, save_iterations_data:bool=True, save_iterations_mesh:bool=False, nonconformal:bool=True):
        """
        Enable automatic mesh refinement (AMR) when simulating the model. See further details in the
        `AWS Palace reference <https://awslabs.github.io/palace/dev/config/model/#model[%22Refinement%22]>`__.
        Some of the parameter definitions have been copied from this reference.

        This function must be run before calling :func:`~SQDMetal.PALACE.Model.PALACE_Model.prepare_simulation`.

        Parameters
        ----------
        num_iterations : int
            Number of AMR iterations to run. Each simulation is run and stored in a new sub-folder as the
            mesh is refined during each iteration.
        max_DoFs : int 
            If set to a value above zero, the AMR iterations will cease to refine the mesh upon reaching
            this value for DoF (degrees of freedom). If set to zero (default value), then there is no limit
            to the maximum DoF (ultimately limited by RAM).
        tolerance : float
            Relative error convergence tolerance for adaptive mesh refinement (AMR). Defaults to 1e-2.
        Dorfler_marking_fraction : float
            Dörfler marking fraction used to specify which elements to refine. It defines smallest number
            of elements that make up this fractional value to the total error in the mesh. A larger value
            will refine more elements per iteration, at the cost of the final mesh being less efficient.
            Default value is 0.7.
        save_iterations_data : bool
            If True, the simulation output is stored for every AMR iteration. Defaults to True.
        save_iterations_mesh : bool
            If True, the mesh is stored for every AMR iteration. Defaults to True.
        nonconformal : bool
            Selects whether AMR should use nonconformal refinement. Nonconformal refinement is required
            for "non-simplex meshes". Somewhat required even if the seed mesh is conformal. Thus, the
            default value is True.

        Returns
        -------
            Returns nothing.
        """
        #TODO: look into how the mesh can be viewed/converted in Gmsh?!
        self._mesh_refinement['Tol'] = tolerance
        self._mesh_refinement['MaxIts'] = num_iterations
        self._mesh_refinement['MaxSize'] =  max_DoFs
        self._mesh_refinement['UpdateFraction'] = Dorfler_marking_fraction
        self._mesh_refinement['UniformLevels'] = 0
        self._mesh_refinement['Nonconformal'] = nonconformal
        self._mesh_refinement['SaveAdaptIterations'] = save_iterations_data
        self._mesh_refinement['SaveAdaptMesh'] = save_iterations_mesh

class PALACE_Model_RF_Base(PALACE_Model):
    def _prepare_simulation(self, metallic_layers, ground_plane):
        '''set-up the simulation'''
        

        if self.meshing == 'GMSH':
            #Prepare the ports...
            #assert len(self._ports) > 0, "There must be at least one port in the RF simulation - do so via the create_port_CPW_on_Launcher or create_port_CPW_on_Route function."
            lePorts = []
            #TODO: Move the coordinate creation code to here - i.e. don't hard-code it on calling the port functions; they store enough metadata anyway...
            for cur_port in self._ports:
                if cur_port['type'] == 'waveport':
                    lePorts.append({'type':'wave', 'metadata':cur_port})
                elif cur_port['elem_type'] == 'single':
                    lePorts.append({'type':'lumped', 'name':cur_port['port_name'], 'coords':cur_port['portCoords']})
                else:
                    lePorts.append({'type':'lumped', 'name':cur_port['port_name']+'a', 'coords':cur_port['portAcoords']})
                    lePorts.append({'type':'lumped', 'name':cur_port['port_name']+'b', 'coords':cur_port['portBcoords']})

            ggb = GMSH_Geometry_Builder(self._geom_processor, self.user_options['fillet_resolution'], self.user_options['gmsh_verbosity'])
            gmsh_render_attrs = ggb.construct_geometry_in_GMSH(self._metallic_layers, self._ground_plane, lePorts,
                                                               self._fine_meshes, self.user_options["fuse_threshold"],
                                                               threshold=self.user_options["threshold"],
                                                               simplify_edge_min_angle_deg=self.user_options["simplify_edge_min_angle_deg"],
                                                               full_3D_params = self._full_3D_params,
                                                               boundary_distances = self._boundary_distances)
            
            gmb = GMSH_Mesh_Builder(gmsh_render_attrs['fine_mesh_elems'], self.user_options)
            gmb.build_mesh()

            if self.create_files:
                #create directory to store simulation files
                self._create_directory(self.name)

                #create config file
                self._create_config_file(gmsh_render_attrs = gmsh_render_attrs)

                # #create batch file
                # if self.mode == 'HPC':
                #     self.create_batch_file()

                self._save_mesh_gmsh()

            # abhishekchak52: commented out for now since pgr is not defined
            # if self.view_design_gmsh_gui == True:
            #     #plot design in gmsh gui
            #     pgr.view_design_components() # TODO: error here: pgr not defined
                    
        if self.meshing == 'COMSOL':

            #initialise the COMSOL engine
            COMSOL_Model.init_engine()
            cmsl = COMSOL_Model('res')

            try:
                #Create COMSOL RF sim object
                sim_sParams = COMSOL_Simulation_RFsParameters(cmsl, adaptive = 'None')
                cmsl.initialize_model(self._geom_processor.design, [sim_sParams], bottom_grounded = True, resolution = 10)     #TODO: Make COMSOL actually compatible rather than assuming Qiskit-Metal designs?

                #Add metallic layers
                for cur_layer in metallic_layers:
                    if cur_layer.get('type') == 'design_layer':
                        cmsl.add_metallic(**cur_layer)
                    elif cur_layer.get('type') == 'Uclip':
                        if cur_layer['clip_type'] == 'inplaneLauncher':
                            sim_sParams.create_CPW_feed_Uclip_on_Launcher(cur_layer['qObjName'], cur_layer['thickness_side'], cur_layer['thickness_back'], cur_layer['separation_gap'])
                        elif cur_layer['clip_type'] == 'inplaneRoute':
                            sim_sParams.create_CPW_feed_Uclip_on_Route(cur_layer['route_name'], cur_layer['pin_name'], cur_layer['thickness_side'], cur_layer['thickness_back'], cur_layer['separation_gap'])
                if not ground_plane['omit']:
                    cmsl.add_ground_plane(threshold=ground_plane['threshold'])
                #cmsl.fuse_all_metals()

                #assign ports
                for cur_port in self._ports:
                    if cur_port['type'] == 'launcher':
                        sim_sParams.create_port_CPW_on_Launcher(cur_port['qObjName'], cur_port['len_launch'])
                    elif cur_port['type'] == 'route':
                        sim_sParams.create_port_CPW_on_Route(cur_port['qObjName'], cur_port['pin_name'], cur_port['len_launch'])
                    elif cur_port['type'] == 'single_rect':
                        sim_sParams.create_port_2_conds_by_position(cur_port['pos1'], cur_port['pos2'], cur_port['rect_width']) #TODO: Make COMSOL support cpw_point CPW feed...
                
                for fine_mesh in self._fine_meshes:
                    if fine_mesh['type'] == 'box':
                        x_bnds = fine_mesh['x_bnds']
                        y_bnds = fine_mesh['y_bnds']
                        cmsl.fine_mesh_in_rectangle(x_bnds[0], y_bnds[0], x_bnds[1], y_bnds[1], fine_mesh['min_size'], fine_mesh['max_size'])
                    elif fine_mesh['type'] == 'comp_bounds':
                        cmsl.fine_mesh_components(fine_mesh['list_comp_names'], fine_mesh['min_size'], fine_mesh['max_size'])

                #build model
                cmsl.build_geom_mater_elec_mesh(mesh_structure = self.user_options["comsol_meshing"])

                #plot model
                # cmsl.plot()

                #save comsol file
                #cmsl.save(self.name)

                if self.create_files:
                    #create directory to store simulation files
                    self._create_directory(self.name)

                    #create config file
                    self._create_config_file(comsol_obj = cmsl, simRF_object = sim_sParams)

                    #create batch file
                    if self.mode == 'HPC':
                        self.create_batch_file()

                #save mesh
                self._save_mesh_comsol(comsol_obj = cmsl)

            except Exception as error:
                cmsl.save("ERROR_" + self.name)
                assert False, f"COMSOL threw an error (file has been saved): {error}"


    def create_port_2_conds(self, qObjName1, pin1, qObjName2, pin2, rect_width=20e-6, impedance_R=50, impedance_L=0, impedance_C=0):
        self._check_if_QiskitMetalDesign()
        
        port_name = "rf_port_" + str(len(self._ports))

        unit_conv = QUtilities.get_units(self._geom_processor.design)
        pos1 = self._geom_processor.design.components[qObjName1].pins[pin1]['middle'] * unit_conv
        pos2 = self._geom_processor.design.components[qObjName2].pins[pin2]['middle'] * unit_conv
        v_parl = pos2-pos1
        v_parl /= np.linalg.norm(v_parl)

        portCoords = ShapelyEx.rectangle_from_line(pos1, pos2, rect_width, False)
        portCoords = [x for x in portCoords]
        portCoords = portCoords + [portCoords[0]]   #Close loop...

        self._ports += [{'port_name':port_name, 'type':'single_rect', 'elem_type':'single', 'pos1':pos1, 'pin2':pos2, 'rect_width': rect_width,
                         'vec_field': v_parl.tolist(),
                         'portCoords': portCoords,
                         'impedance_R':impedance_R, 'impedance_L':impedance_L, 'impedance_C':impedance_C}]

    def create_port_JosephsonJunction(self, qObjName, **kwargs):
        self._check_if_QiskitMetalDesign()

        junction_index = kwargs.get('junction_index', 0)

        comp_id = self._geom_processor.design.components[qObjName].id
        gsdf = self._geom_processor.design.qgeometry.tables['junction']
        gsdf = gsdf.loc[gsdf["component"] == comp_id]
        ls = gsdf['geometry'].iloc[junction_index]
        assert isinstance(ls, shapely.geometry.linestring.LineString), "The junction must be defined as a LineString object. Check the code for the part to verify this..."

        coords = [x for x in ls.coords]
        assert len(coords) == 2, "Currently only junctions with two points (i.e. a rectangle) are supported."
        rect_width = gsdf['width'].iloc[junction_index]

        if 'E_J_Hertz' in kwargs:
            elem_charge = 1.60218e-19
            hbar = 1.05457e-34
            h = hbar*2*np.pi
            phi0 = hbar/(2*elem_charge)
            E_J = kwargs.pop('E_J_Hertz') * h
            L_ind = phi0**2 / E_J
        elif 'L_J' in kwargs:
            L_ind = kwargs.pop('L_J')
        else:
            assert False, "Must supply either E_J_Hertz or L_J to define a Josephson Junction port."

        C_J = kwargs.get('C_J', 0)

        port_name = "rf_port_" + str(len(self._ports))

        unit_conv = QUtilities.get_units(self._geom_processor.design)
        pos1 = np.array(coords[0]) * unit_conv
        pos2 = np.array(coords[1]) * unit_conv
        v_parl = pos2-pos1
        v_parl /= np.linalg.norm(v_parl)

        portCoords = ShapelyEx.rectangle_from_line(pos1, pos2, rect_width * unit_conv, False)
        portCoords = [x for x in portCoords]
        portCoords = portCoords + [portCoords[0]]   #Close loop...

        self._ports += [{'port_name':port_name, 'type':'single_rect', 'elem_type':'single', 'pos1':pos1, 'pin2':pos2, 'rect_width': rect_width * unit_conv,
                         'vec_field': v_parl.tolist(),
                         'portCoords': portCoords,
                         'impedance_R':0, 'impedance_L':L_ind, 'impedance_C':C_J}]

    def create_port_CPW_on_Launcher(self, qObjName, len_launch = 20e-6, impedance_R=50, impedance_L=0, impedance_C=0):
        '''
        Creates an RF port on a CPW inlet. The elements form fins to ground from the central CPW stripline. Note that the first port will automatically be
        the 50Ohm excitation port in the E-field plots while the second port is a 50Ohm ground. The s-parameters will calculate S11 and S21 if 2 such ports
        are defined.
        Inputs:
            - CPW_obj - A CPW object that has the attributes: start, end, width, gap
            - is_start - If True, then the port is attached to the start of the CPW, while False attaches the port to the end of the CPW
            - len_launch - (Default: 20e-6) Length of the inlet port fins along the the CPW. It is a good idea to keep it thin w.r.t. CPW gap 
        '''
        self._check_if_QiskitMetalDesign()
        port_name = "rf_port_" + str(len(self._ports))
        
        launchesA, launchesB, vec_perp = QUtilities.get_RFport_CPW_coords_Launcher(self._geom_processor.design, qObjName, len_launch, 1)  #Units of m...

        #See here for details: https://awslabs.github.io/palace/stable/config/boundaries/#boundaries[%22LumpedPort%22]
        self._ports += [{'port_name':port_name, 'type':'launcher', 'elem_type':'cpw', 'qObjName':qObjName, 'len_launch': len_launch,
                         'portAcoords': launchesA + [launchesA[0]],
                         'portBcoords': launchesB + [launchesB[0]],
                         'vec_field': vec_perp.tolist(),
                         'impedance_R':impedance_R, 'impedance_L':impedance_L, 'impedance_C':impedance_C}]

    def create_port_CPW_on_Route(self, qObjName, pin_name='end', len_launch = 20e-6, impedance_R=50, impedance_L=0, impedance_C=0):
        self._check_if_QiskitMetalDesign()
        port_name = "rf_port_" + str(len(self._ports))
        
        launchesA, launchesB, vec_perp = QUtilities.get_RFport_CPW_coords_Route(self._geom_processor.design, qObjName, pin_name, len_launch, 1)  #Units of m...

        #See here for details: https://awslabs.github.io/palace/stable/config/boundaries/#boundaries[%22LumpedPort%22]
        self._ports += [{'port_name':port_name, 'type':'route', 'elem_type':'cpw', 'qObjName':qObjName, 'pin_name':pin_name, 'len_launch': len_launch,
                         'portAcoords': launchesA + [launchesA[0]],
                         'portBcoords': launchesB + [launchesB[0]],
                         'vec_field': vec_perp.tolist(),
                         'impedance_R':impedance_R, 'impedance_L':impedance_L, 'impedance_C':impedance_C}]

    def create_port_CPW_via_edge_point(self, pt_near_centre_end, len_launch, impedance_R=50, impedance_L=0, impedance_C=0, **kwargs):
        port_name = "rf_port_" + str(len(self._ports))

        kwargs["fuse_threshold"] = self.user_options["fuse_threshold"]
        kwargs["threshold"] = self.user_options["threshold"]
        launchesA, launchesB, vec_perp = self._geom_processor.create_CPW_feed_via_point(pt_near_centre_end, len_launch, self._metallic_layers, self._ground_plane, **kwargs)
        
        #See here for details: https://awslabs.github.io/palace/stable/config/boundaries/#boundaries[%22LumpedPort%22]
        self._ports += [{'port_name':port_name, 'type':'point', 'elem_type':'cpw', 'len_launch': len_launch,
                         'portAcoords': launchesA + [launchesA[0]],
                         'portBcoords': launchesB + [launchesB[0]],
                         'vec_field': vec_perp,
                         'impedance_R':impedance_R, 'impedance_L':impedance_L, 'impedance_C':impedance_C}]


    def set_port_impedance(self, port_ind, impedance_R=50, impedance_L=0, impedance_C=0):
        assert self._ports[port_ind-1]['type'] != 'waveport', "Cannot set impedances to a waveport."
        #Enumerate port_ind from 1...
        #TODO: Override if different for Eigenmode...
        self._ports[port_ind-1]['impedance_R'] = impedance_R
        self._ports[port_ind-1]['impedance_L'] = impedance_L
        self._ports[port_ind-1]['impedance_C'] = impedance_C
        if self._sim_config != "":
            with open(self._sim_config, "r") as f:
                config_json = json.loads(f.read())
            config_json['Boundaries']['LumpedPort'][port_ind-1]['R'] = impedance_R
            config_json['Boundaries']['LumpedPort'][port_ind-1]['L'] = impedance_L
            config_json['Boundaries']['LumpedPort'][port_ind-1]['C'] = impedance_C
            with open(self._sim_config, "w") as f:
                json.dump(config_json, f, indent=2)

    def create_waveport_on_boundary(self, plane:str, **kwargs):
        port_name = "rf_wport_" + str(len(self._ports))

        shape = kwargs.get('shape', 'rectangle')
        supported_shapes = ['rectangle']
        assert shape in supported_shapes, f"The shape '{shape}' is unsupported."

        if plane == 'x_pos' or plane == 'x_neg':
            assert 'y' in kwargs and 'z' in kwargs, "Must supply y and z arguments when placing waveport on positive x planar boundary."
            assert 'size_y' in kwargs and 'size_z' in kwargs, "Must supply size_y and size_z arguments when placing waveport on a x planar boundary."
            new_wvprt = {'y': kwargs['y'], 'z': kwargs['z'], 'size_y':kwargs['size_y'], 'size_z':kwargs['size_z']}
        elif plane == 'y_pos' or plane == 'y_neg':
            assert 'x' in kwargs and 'z' in kwargs, "Must supply x and z arguments when placing waveport on positive x planar boundary."
            assert 'size_x' in kwargs and 'size_z' in kwargs, "Must supply size_x and size_z arguments when placing waveport on a x planar boundary."
            new_wvprt = {'x': kwargs['x'], 'z': kwargs['z'], 'size_x':kwargs['size_x'], 'size_z':kwargs['size_z']}
        elif plane == 'z_pos' or plane == 'z_neg':
            assert 'x' in kwargs and 'y' in kwargs, "Must supply x and y arguments when placing waveport on positive x planar boundary."
            assert 'size_x' in kwargs and 'size_y' in kwargs, "Must supply size_y and size_z arguments when placing waveport on a z planar boundary."
            new_wvprt = {'x': kwargs['x'], 'y': kwargs['y'], 'size_x':kwargs['size_x'], 'size_y':kwargs['size_y']}
        else:
            assert False, "The argument 'plane' must be 'x_pos', 'y_pos', 'z_pos', 'x_neg', 'y_neg' or 'z_neg'."

        new_wvprt['shape'] = shape
        new_wvprt['plane'] = plane
        new_wvprt['port_name'] = port_name
        new_wvprt['type'] = 'waveport'
        self._ports.append(new_wvprt)

    def create_CPW_feed_Uclip_on_Launcher(self, qObjName, thickness_side=20e-6, thickness_back=20e-6, separation_gap=0e-6):
        self._metallic_layers += [{
            'type': 'Uclip',
            'clip_type':'inplaneLauncher',
            'qObjName':qObjName,
            'thickness_side':thickness_side,
            'thickness_back':thickness_back,
            'separation_gap':separation_gap
        }]

    def create_CPW_feed_Uclip_on_Route(self, route_name, pin_name, thickness_side=20e-6, thickness_back=20e-6, separation_gap=0e-6):
        self._metallic_layers += [{
            'type': 'Uclip',
            'clip_type':'inplaneRoute',
            'route_name':route_name,
            'pin_name':pin_name,
            'thickness_side':thickness_side,
            'thickness_back':thickness_back,
            'separation_gap':separation_gap
        }]
        #TODO: Do this in general for GDS as well?
        #TODO: Also consider a single-conductor feed onto the U-Clip?

    def _process_ports(self, ports):
        #Assumes that ports is a dictionary that contains the port names (with separate keys with suffixes a and b for multi-element ports)
        #where each value is a list of element IDs corresponding to the particular port...
        config_ports = []
        config_wports = []
        for m, cur_port in enumerate(self._ports):
            if cur_port['type'] == 'waveport':
                leDict = {
                        "Index": m+1,                    
                    }

                leDict['Attributes'] = [ports[cur_port['port_name']]]
                config_wports.append(leDict)
            else:
                port_name, vec_field = cur_port['port_name'], cur_port['vec_field']
                leDict = {
                        "Index": m+1,                    
                    }
                if port_name + 'a' in ports:
                    leDict['Elements'] = [
                            {
                            "Attributes": [ports[port_name + 'a']],
                            "Direction": vec_field + [0]
                            },
                            {
                            "Attributes": [ports[port_name + 'b']],
                            "Direction": [-x for x in vec_field] + [0]
                            }
                        ]
                else:
                    leDict['Attributes'] = [ports[port_name]]
                    leDict['Direction'] = vec_field + [0]

                if 'impedance_R' in cur_port:
                    leDict['R'] = cur_port['impedance_R']
                if 'impedance_L' in cur_port:
                    leDict['L'] = cur_port['impedance_L']
                if 'impedance_C' in cur_port:
                    leDict['C'] = cur_port['impedance_C']
                config_ports.append(leDict)
        return config_ports, config_wports

    def setup_EPR_interfaces(self, substrate_air : MaterialInterface, substrate_metal : MaterialInterface, metal_air : MaterialInterface, **kwargs):
        self.substrate_air = substrate_air
        self.substrate_metal = substrate_metal
        self.metal_air = metal_air
        self._EPR_setup = True
        self.substrate_air_thickness = kwargs.get('substrate_air_thickness', 2e-9)
        self.substrate_metal_thickness = kwargs.get('substrate_air_thickness', 2e-9)
        self.metal_air_thickness = kwargs.get('substrate_air_thickness', 2e-9)

    def _setup_EPR_boundaries(self, dict_json, id_dielectric_gaps, id_metals, **kwargs):
        if not self._EPR_setup:
            return
        len_scale = kwargs.get('len_scale', 1e-3)    #Only supports GMSH
        dict_json['Boundaries']['Postprocessing'] = {
                   "Dielectric":
                    [
                        {
                        "Index": 1,
                        "Attributes": id_dielectric_gaps,
                        "Type": "SA",
                        "Thickness": self.substrate_air_thickness/len_scale,
                        "Permittivity": self.substrate_air.permittivity,
                        "LossTan": self.substrate_air.loss_tangent
                        },
                        {
                        "Index": 2,
                        "Attributes": id_metals,
                        "Type": "MS",
                        "Thickness": self.substrate_metal_thickness/len_scale,  
                        "Permittivity": self.substrate_metal.permittivity,
                        "LossTan": self.substrate_metal.loss_tangent
                        },
                        {
                        "Index": 3,
                        "Attributes": id_metals,
                        "Type": "MA",
                        "Thickness": self.metal_air_thickness/len_scale, 
                        "Permittivity": self.metal_air.permittivity,
                        "LossTan": self.metal_air.loss_tangent
                        }
                    ] 
                }
        
    def add_kinetic_inductance(self, kinetic_inductance):
        '''Method for user to add kinetic inductance into simulations'''
        self._use_KI = True
        self._KI = kinetic_inductance
        
        
    def _setup_kinetic_inductance(self, dict_json, id_metals):
        '''Internal method to adjust Palace config file to change boundary condition from PEC to Surface Impedance to
            incorporate kinetic inductance.'''

        dict_json['Boundaries']['Impedance'] = [{
                    "Attributes": id_metals,
                    "Ls": self._KI                   
                }]
        
        dict_json['Boundaries']['PEC'] = {
                    "Attributes": [],
                }

    def set_farfield(self, ff_type='pec', ff_plane='', ff_material: MaterialConductor=None):
        #ff_type can be: 'absorbing', 'pec' or 'conductor'
        ff_type = ff_type.lower()
        assert ff_type == 'pec' or ff_type == 'absorbing' or ff_type == 'conductor', "ff_type must be: 'absorbing', 'pec' or 'conductor'"
        ff_planes = ['x_pos', 'x_neg', 'y_pos', 'y_neg', 'z_pos', 'z_neg']
        if ff_plane != '':
            assert ff_plane in ff_planes, "Farfield plane must be 'x_pos', 'x_neg', 'y_pos', 'y_neg', 'z_pos', 'z_neg'."
            ff_planes = [ff_plane]
        for cur_ff_plane in ff_planes:
            if ff_type == 'conductor':
                assert ff_material != None, "Must supply a MaterialConductor object when specifying farfields as a conductor."
                cur_ff_type = (ff_type, ff_material)
            else:
                cur_ff_type = (ff_type,None)
            self._ff_type[cur_ff_plane] = cur_ff_type

    def _process_farfield(self, dict_config, ff_physgrps):
        abs_list = []
        pec_list = []
        con_list = []
        for ff_plane in self._ff_type:
            if self._ff_type[ff_plane][0] == 'absorbing':
                abs_list.append(ff_physgrps[ff_plane])
            elif self._ff_type[ff_plane][0] == 'pec':
                pec_list.append(ff_physgrps[ff_plane])
            elif self._ff_type[ff_plane][0] == 'conductor':
                con_list.append((ff_physgrps[ff_plane], self._ff_type[ff_plane][1]))
        #
        if len(abs_list) > 0:
            dict_config['Boundaries']['Absorbing'] = {"Attributes": abs_list, "Order": 1}
        if len(pec_list) > 0:
            dict_config['Boundaries']['PEC']['Attributes'] += pec_list
        for cur_cond_plane in con_list:
            if not 'Conductivity' in dict_config['Boundaries']:
                dict_config['Boundaries']['Conductivity'] = []
            attrs, material_cond = cur_cond_plane
            dict_config['Boundaries']['Conductivity'].append(
                {
                    'Attributes': [attrs],
                    'Conductivity': material_cond.conductivity,
                    'Permeability': material_cond.permeability
                    # 'Thickness': <float>
                })
