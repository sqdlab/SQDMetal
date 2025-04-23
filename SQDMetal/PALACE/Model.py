from SQDMetal.Utilities.QUtilities import QUtilities
from SQDMetal.Utilities.ShapelyEx import ShapelyEx
from SQDMetal.Utilities.Materials import MaterialInterface
from SQDMetal.COMSOL.Model import COMSOL_Model
from SQDMetal.COMSOL.SimRFsParameter import COMSOL_Simulation_RFsParameters
from SQDMetal.PALACE.Utilities.GMSH_Geometry_Builder import GMSH_Geometry_Builder
from SQDMetal.PALACE.Utilities.GMSH_Mesh_Builder import GMSH_Mesh_Builder
from SQDMetal.PALACE.Utilities.GMSH_Navigator import GMSH_Navigator
import numpy as np
import os, subprocess
import gmsh
import json
import platform
import shapely
import shutil

class PALACE_Model:
    def __init__(self, meshing, mode, options):
        self.meshing = meshing
        self._metallic_layers = []
        self._ground_plane = {'omit': True}
        self._fine_meshes = []
        self._sim_config = ""
        self._input_dir = ""
        self._output_subdir = ""
        self.set_farfield()
        self._EPR_setup = False
        self._rf_port_excitation = -1

        if mode == 'HPC':
            with open(options["HPC_Parameters_JSON"], "r") as f:
                self.hpc_options = json.loads(f.read())
        else:
            self.palace_dir = options.get('palace_dir', 'palace')                
            self.hpc_options = {"input_dir":""}
        self._num_cpus = options.get('num_cpus', 16)

        #Assuming that metal_design attribute already exists

        restrict_rect = options.get('segment_rectangle', None)   #Given as [x1,y1,x2,y2]
        if isinstance(restrict_rect, list) or isinstance(restrict_rect, np.ndarray) or isinstance(restrict_rect, tuple):
            self.restrict_rect = [min([restrict_rect[0], restrict_rect[2]]), min([restrict_rect[1], restrict_rect[3]]),
                                  max([restrict_rect[0], restrict_rect[2]]), max([restrict_rect[1], restrict_rect[3]])]
            self.chip_len = self.restrict_rect[2]-self.restrict_rect[0]
            self.chip_wid = self.restrict_rect[3]-self.restrict_rect[1]
            self.chip_centre = [(self.restrict_rect[0]+self.restrict_rect[2])*0.5,
                                (self.restrict_rect[1]+self.restrict_rect[3])*0.5,
                                QUtilities.parse_value_length(self.metal_design.chips['main']['size']['center_z'])]
        else:
            self.chip_len = QUtilities.parse_value_length(self.metal_design.chips['main']['size']['size_x'])
            self.chip_wid = QUtilities.parse_value_length(self.metal_design.chips['main']['size']['size_y'])
            self.chip_centre = [QUtilities.parse_value_length(self.metal_design.chips['main']['size'][x]) for x in ['center_x', 'center_y', 'center_z']]
            self.restrict_rect = [self.chip_centre[0]-0.5*self.chip_len, self.chip_centre[1]-0.5*self.chip_wid,
                                  self.chip_centre[0]+0.5*self.chip_len, self.chip_centre[1]+0.5*self.chip_wid]

    def create_batch_file(self):
        pass
    def create_config_file(self):
        pass
    def _prepare_simulation(self, metallic_layers, ground_plane):
        raise NotImplementedError()

    def prepare_simulation(self):
        self._prepare_simulation(self._metallic_layers, self._ground_plane)

    def add_metallic(self, layer_id, **kwargs):
        '''
        Adds metallic conductors from the Qiskit-Metal design object onto the surface layer of the chip simulation. If the particular layer has
        fancy PVD evaporation steps, the added metallic layer will account for said steps and merge the final result. In addition, all metallic
        elements that are contiguous are merged into single blobs.

        Inputs:
            - layer_id - The index of the layer from which to take the metallic polygons
            - threshold - (Optional) Defaults to -1. This is the threshold in metres, below which consecutive vertices along a given polygon are
                          combined into a single vertex. This simplification helps with meshing as COMSOL will not overdo the meshing. If this
                          argument is negative, the argument is ignored.
            - fuse_threshold - (Optional) Defaults to 1e-12. This is the minimum distance between metallic elements, below which they are considered
                               to be a single polygon and thus, the polygons are merged with the gap filled. This accounts for floating-point errors
                               that make adjacent elements fail to merge as a single element, due to infinitesimal gaps between them.
            - evap_mode - (Optional) Defaults to 'separate_delete_below'. These are the methods upon which to separate or merge overlapping elements
                          across multiple evaporation steps. See documentation on PVD_Shadows for more details on the available options.
            - group_by_evaporations - (Optional) Defaults to False. If set to True, if elements on a particular evaporation step are separated due
                                      to the given evap_mode, they will still be selected as a part of the same conductor (useful for example, in
                                      capacitance matrix simulations).
            - evap_trim - (Optional) Defaults to 20e-9. This is the trimming distance used in certain evap_mode profiles. See documentation on
                          PVD_Shadows for more details on its definition.
        '''
        self._metallic_layers += [{
            'type': 'design_layer',
            'layer_id': layer_id,
            'threshold': kwargs.get('threshold', -1),
            'fuse_threshold': kwargs.get('fuse_threshold', 1e-12),
            'evap_mode': kwargs.get('evap_mode', 'separate_delete_below'),
            'group_by_evaporations': kwargs.get('group_by_evaporations', False),
            'evap_trim': kwargs.get('evap_trim', 20e-9),
        }]

    def add_ground_plane(self, **kwargs):
        '''
        Adds metallic ground-plane from the Qiskit-Metal design object onto the surface layer of the chip simulation.

        Inputs:
            - threshold - (Optional) Defaults to -1. This is the threshold in metres, below which consecutive vertices along a given polygon are
                          combined into a single vertex. This simplification helps with meshing as COMSOL will not overdo the meshing. If this
                          argument is negative, the argument is ignored.
        '''
        self._ground_plane = {'omit':False, 'threshold':kwargs.get('threshold', -1)}

    def set_farfield(self, ff_type='absorbing'):
        #ff_type can be: 'absorbing' or 'pec'
        ff_type = ff_type.lower()
        assert ff_type == 'pec' or ff_type == 'absorbing', "ff_type must be: 'absorbing' or 'pec'"
        self._ff_type = ff_type

    def _is_native_arm64(self):
        '''
        Helper method written by Sadman Ahmed Shanto @shanto268 in this issue: https://github.com/sqdlab/SQDMetal/issues/9
        '''
        try:
            # Run the sysctl command to check the native architecture
            result = subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True, check=True)
            return "M" in result.stdout  # M is for Apple M1/M2 chips
        except subprocess.CalledProcessError as e:
            # print(f"Error checking native architecture: {e}")
            return False

    def run(self):
        assert self._sim_config != "", "Must run prepare_simulation at least once."

        config_file = self._sim_config
        leFile = os.path.basename(os.path.realpath(config_file))
        leDir = os.path.dirname(os.path.realpath(config_file))

        if not os.path.exists(self._output_data_dir):
            os.makedirs(self._output_data_dir)
        log_location = f"{self._output_data_dir}/out.log"
        with open("temp.sh", "w+") as f:
            f.write(f"cd \"{leDir}\"\n")
            f.write(f"{self.palace_dir} -np {self._num_cpus} {leFile} | tee \"{log_location}\"\n")

        # Set execute permission on temp.sh
        os.chmod("temp.sh", 0o755)

        with open(log_location, 'w') as fp:
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

        shutil.copy(self._sim_config, self._output_data_dir + f'/config.json')

        return self.retrieve_data()

    def retrieve_data(self):
        pass

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

    def _check_simulation_mode(self):
        '''method to check the type of simualtion being run and return
            the parent directory to store the simulation files'''

        parent_simulation_dir = None

        if self.mode == "HPC":
            parent_simulation_dir = self.sim_parent_directory
        elif self.mode == "simPC" or self.mode == "PC":
            parent_simulation_dir = self.sim_parent_directory
        else:
            Exception('Invalid simulation mode entered.')
        
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
        if self.create_files == True:
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
            
    def _save_mesh_comsol(self, comsol_obj):
        '''function used to save the comsol mesh file'''
        if self.create_files == True:
            parent_simulation_dir = self._check_simulation_mode()

            # file_name
            comsol_obj.save(parent_simulation_dir + self.name + "/" + self.name)
    
            # Path
            path = os.path.abspath(parent_simulation_dir) + "/" + self.name + "/" + self.name + ".mphbin"

            #COMSOL export commands
            comsol_obj._model.java.component("comp1").mesh("mesh1").export().set("filename", path)
            comsol_obj._model.java.component("comp1").mesh("mesh1").export(path)

    def _get_folder_prefix(self):
        return self.hpc_options["input_dir"]  + self.name + "/" if self.hpc_options["input_dir"] != "" else ""

    def set_local_output_subdir(self, name, update_config_file=True):
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

    def fine_mesh_along_path(self, dist_resolution, qObjName, trace_name='', **kwargs):
        print("Try not to use fine_mesh_along_path for now. A bug makes multiple calls of this mesh around the edges...")
        leUnits = QUtilities.get_units(self.metal_design)
        lePath = QUtilities.calc_points_on_path(dist_resolution/leUnits, self.metal_design, qObjName, trace_name)[0] * leUnits
        self._fine_meshes.append({
            'type': 'path',
            'path': lePath,
            'min_size': kwargs.get('min_size', self.user_options['mesh_min']),
            'max_size': kwargs.get('max_size', self.user_options['mesh_max']),
            'taper_dist_min': kwargs.get('taper_dist_min', self.user_options['taper_dist_min']),
            'taper_dist_max': kwargs.get('taper_dist_max', self.user_options['taper_dist_max'])
        })

    def fine_mesh_in_rectangle(self, x1, y1, x2, y2, **kwargs):
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
        self._fine_meshes.append({
            'type': 'comp_bounds',
            'list_comp_names': list_comp_names,
            'min_size': kwargs.get('min_size', self.user_options['mesh_min']),
            'max_size': kwargs.get('max_size', self.user_options['mesh_max']),
            'taper_dist_min': kwargs.get('taper_dist_min', self.user_options['taper_dist_min']),
            'taper_dist_max': kwargs.get('taper_dist_max', self.user_options['taper_dist_max']),
            'metals_only': kwargs.get('metals_only', False)
        })

class PALACE_Model_RF_Base(PALACE_Model):
    def _prepare_simulation(self, metallic_layers, ground_plane):
        '''set-up the simulation'''
        
        if self.meshing == 'GMSH':
            #Prepare the ports...
            assert len(self._ports) > 0, "There must be at least one port in the RF simulation - do so via the create_port_CPW_on_Launcher or create_port_CPW_on_Route function."
            lePorts = []
            for cur_port in self._ports:
                if cur_port['elem_type'] == 'single':
                    lePorts += [(cur_port['port_name'], cur_port['portCoords'])]
                else:
                    lePorts += [(cur_port['port_name'] + 'a', cur_port['portAcoords'])]
                    lePorts += [(cur_port['port_name'] + 'b', cur_port['portBcoords'])]

            ggb = GMSH_Geometry_Builder(self.metal_design, self.user_options['fillet_resolution'], self.user_options['gmsh_verbosity'])
            gmsh_render_attrs = ggb.construct_geometry_in_GMSH(self._metallic_layers, self._ground_plane, lePorts, self._fine_meshes, self.user_options["fuse_threshold"], threshold=self.user_options["threshold"])
            #
            gmb = GMSH_Mesh_Builder(gmsh_render_attrs['fine_mesh_elems'], self.user_options)
            gmb.build_mesh()

            if self.create_files == True:
                #create directory to store simulation files
                self._create_directory(self.name)

                #create config file
                self.create_config_file(gmsh_render_attrs = gmsh_render_attrs)

                # #create batch file
                # if self.mode == 'HPC':
                #     self.create_batch_file()

                self._save_mesh_gmsh()

            if self.view_design_gmsh_gui == True:
                #plot design in gmsh gui
                pgr.view_design_components() # TODO: error here: pgr not defined
                    
        if self.meshing == 'COMSOL':

            #initialise the COMSOL engine
            COMSOL_Model.init_engine()
            cmsl = COMSOL_Model('res')

            try:
                #Create COMSOL RF sim object
                sim_sParams = COMSOL_Simulation_RFsParameters(cmsl, adaptive = 'None')
                cmsl.initialize_model(self.metal_design, [sim_sParams], bottom_grounded = True, resolution = 10)

                #Add metallic layers
                for cur_layer in metallic_layers:
                    if cur_layer.get('type') == 'design_layer':
                        cmsl.add_metallic(**cur_layer)
                    elif cur_layer.get('type') == 'Uclip':
                        if cur_layer['clip_type'] == 'inplaneLauncher':
                            sim_sParams.create_RFport_CPW_groundU_Launcher_inplane(cur_layer['qObjName'], cur_layer['thickness_side'], cur_layer['thickness_back'], cur_layer['separation_gap'], cur_layer['unit_conv_extra'])
                        elif cur_layer['clip_type'] == 'inplaneRoute':
                            sim_sParams.create_RFport_CPW_groundU_Route_inplane(cur_layer['route_name'], cur_layer['pin_name'], cur_layer['thickness_side'], cur_layer['thickness_back'], cur_layer['separation_gap'], cur_layer['unit_conv_extra'])
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
                        sim_sParams.create_port_2_conds_by_position(cur_port['pos1'], cur_port['pos2'], cur_port['rect_width']) 
                
                for fine_mesh in self._fine_meshes:
                    if fine_mesh['type'] == 'box':
                        x_bnds = fine_mesh['x_bnds']
                        y_bnds = fine_mesh['y_bnds']
                        cmsl.fine_mesh_in_rectangle(x_bnds[0], y_bnds[0], x_bnds[1], y_bnds[1], fine_mesh['min_size'], fine_mesh['max_size'])
                    elif fine_mesh['type'] == 'comp_bounds':
                        cmsl.fine_mesh_around_comp_boundaries(fine_mesh['list_comp_names'], fine_mesh['min_size'], fine_mesh['max_size'])

                #build model
                cmsl.build_geom_mater_elec_mesh(mesh_structure = self.user_options["comsol_meshing"])

                #plot model
                # cmsl.plot()

                #save comsol file
                #cmsl.save(self.name)

                if self.create_files == True:
                    #create directory to store simulation files
                    self._create_directory(self.name)

                    #create config file
                    self.create_config_file(comsol_obj = cmsl, simRF_object = sim_sParams)

                    #create batch file
                    if self.mode == 'HPC':
                        self.create_batch_file()

                #save mesh
                self._save_mesh_comsol(comsol_obj = cmsl)

            except Exception as error:
                cmsl.save("ERROR_" + self.name)
                assert False, f"COMSOL threw an error (file has been saved): {error}"

    def set_port_excitation(self, port_index):
        assert port_index > 0 and port_index <= len(self._ports), "Invalid index of port for excitation. Check if ports with R>0 have been correctly defined."
        assert self._ports[port_index-1]['impedance_R'] > 0, f"Port {port_index} does not have a non-zero resistive part in its impedance."
        self._rf_port_excitation = port_index

    def create_port_2_conds(self, qObjName1, pin1, qObjName2, pin2, rect_width=20e-6, impedance_R=50, impedance_L=0, impedance_C=0):
        port_name = "rf_port_" + str(len(self._ports))

        unit_conv = QUtilities.get_units(self.metal_design)
        pos1 = self.metal_design.components[qObjName1].pins[pin1]['middle'] * unit_conv
        pos2 = self.metal_design.components[qObjName2].pins[pin2]['middle'] * unit_conv
        v_parl = pos2-pos1
        v_parl /= np.linalg.norm(v_parl)

        portCoords = ShapelyEx.rectangle_from_line(pos1, pos2, rect_width, False)

        self._ports += [{'port_name':port_name, 'type':'single_rect', 'elem_type':'single', 'pos1':pos1, 'pin2':pos2, 'rect_width': rect_width,
                         'vec_field': v_parl.tolist(),
                         'portCoords': portCoords,
                         'impedance_R':impedance_R, 'impedance_L':impedance_L, 'impedance_C':impedance_C}]
        if self._rf_port_excitation == -1 and impedance_R > 0:
            self._rf_port_excitation = len(self._ports)

    def create_port_JosephsonJunction(self, qObjName, **kwargs):
        junction_index = kwargs.get('junction_index', 0)

        comp_id = self.metal_design.components[qObjName].id
        gsdf = self.metal_design.qgeometry.tables['junction']
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

        unit_conv = QUtilities.get_units(self.metal_design)
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
        port_name = "rf_port_" + str(len(self._ports))
        
        launchesA, launchesB, vec_perp = QUtilities.get_RFport_CPW_coords_Launcher(self.metal_design, qObjName, len_launch, 1)  #Units of m...

        #See here for details: https://awslabs.github.io/palace/stable/config/boundaries/#boundaries[%22LumpedPort%22]
        self._ports += [{'port_name':port_name, 'type':'launcher', 'elem_type':'cpw', 'qObjName':qObjName, 'len_launch': len_launch,
                         'portAcoords': launchesA + [launchesA[0]],
                         'portBcoords': launchesB + [launchesB[0]],
                         'vec_field': vec_perp.tolist(),
                         'impedance_R':impedance_R, 'impedance_L':impedance_L, 'impedance_C':impedance_C}]
        if self._rf_port_excitation == -1 and impedance_R > 0:
            self._rf_port_excitation = len(self._ports)

    def create_port_CPW_on_Route(self, qObjName, pin_name='end', len_launch = 20e-6, impedance_R=50, impedance_L=0, impedance_C=0):
        port_name = "rf_port_" + str(len(self._ports))
        
        launchesA, launchesB, vec_perp = QUtilities.get_RFport_CPW_coords_Route(self.metal_design, qObjName, pin_name, len_launch, 1)  #Units of m...

        #See here for details: https://awslabs.github.io/palace/stable/config/boundaries/#boundaries[%22LumpedPort%22]
        self._ports += [{'port_name':port_name, 'type':'route', 'elem_type':'cpw', 'qObjName':qObjName, 'pin_name':pin_name, 'len_launch': len_launch,
                         'portAcoords': launchesA + [launchesA[0]],
                         'portBcoords': launchesB + [launchesB[0]],
                         'vec_field': vec_perp.tolist(),
                         'impedance_R':impedance_R, 'impedance_L':impedance_L, 'impedance_C':impedance_C}]
        if self._rf_port_excitation == -1 and impedance_R > 0:
            self._rf_port_excitation = len(self._ports)

    def set_port_impedance(self, port_ind, impedance_R=50, impedance_L=0, impedance_C=0):
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

    def create_RFport_CPW_groundU_Launcher_inplane(self, qObjName, thickness_side=20e-6, thickness_back=20e-6, separation_gap=0e-6, unit_conv_extra = 1):
        self._metallic_layers += [{
            'type': 'Uclip',
            'clip_type':'inplaneLauncher',
            'qObjName':qObjName,
            'thickness_side':thickness_side,
            'thickness_back':thickness_back,
            'separation_gap':separation_gap,
            'unit_conv_extra':unit_conv_extra
        }]

    def create_RFport_CPW_groundU_Route_inplane(self, route_name, pin_name, thickness_side=20e-6, thickness_back=20e-6, separation_gap=0e-6, unit_conv_extra = 1):
        self._metallic_layers += [{
            'type': 'Uclip',
            'clip_type':'inplaneRoute',
            'route_name':route_name,
            'pin_name':pin_name,
            'thickness_side':thickness_side,
            'thickness_back':thickness_back,
            'separation_gap':separation_gap,
            'unit_conv_extra':unit_conv_extra
        }]

    def _process_ports(self, ports):
        #Assumes that ports is a dictionary that contains the port names (with separate keys with suffixes a and b for multi-element ports)
        #where each value is a list of element IDs corresponding to the particular port...
        config_ports = []
        for m, cur_port in enumerate(self._ports):
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
        return config_ports

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

