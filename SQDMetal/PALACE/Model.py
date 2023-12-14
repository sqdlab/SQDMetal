from SQDMetal.Utilities.QUtilities import QUtilities
import numpy as np
import os

class PALACE_Model:
    def __init__(self, meshing):
        self.meshing = meshing
        self._metallic_layers = []
        self._ground_plane = {'omit': True}

    def create_batch_file(self):
        pass
    def create_config_file(self):
        pass
    def run(self):
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

    def _check_simulation_mode(self):
        '''method to check the type of simualtion being run and return
            the parent directory to store the simulation files'''

        parent_simulation_dir = None

        if self.mode == "HPC":
            parent_simulation_dir = self.sim_parent_directory
        elif self.mode == "simPC":
            parent_simulation_dir = self.simPC_parent_simulation_dir
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
            os.mkdir(path)
            print("Directory '% s' created" % directory)


    def _save_mesh_gmsh(self):
        '''function used to save the gmsh mesh file'''
        if self.create_files == True:
            parent_simulation_dir = self._check_simulation_mode()

            # file_name
            file_name = self.name + "/" + self.name + ".msh"
    
            # Path
            path = os.path.join(parent_simulation_dir, file_name)
            gmsh.write(path)
            
            
    def _save_mesh_comsol(self, comsol_obj):
        '''function used to save the comsol mesh file'''
        if self.create_files == True:
            parent_simulation_dir = self._check_simulation_mode()

            # file_name
            file_name = self.name + "/" + self.name + ".mphbin"
    
            # Path
            path = os.path.join(parent_simulation_dir, file_name)

            #COMSOL export commands
            comsol_obj._model.java.component("comp1").mesh("mesh1").export().set("filename", path)
            comsol_obj._model.java.component("comp1").mesh("mesh1").export(path)


class PALACE_Model_RF_Base(PALACE_Model):
    def _prepare_simulation(self, metallic_layers, ground_plane):
        '''set-up the simulation'''
        
        if self.meshing == 'GMSH':

            #Create the gmsh renderer to convert qiskit metal geomentry to gmsh geometry
            pgr = Palace_Gmsh_Renderer(self.metal_design)

            #Prepare the ports...
            assert len(self._ports) > 0, "There must be at least one port in the RF simulation - do so via the create_port_CPW_on_Launcher function."
            lePorts = []
            for cur_port in self._ports:
                port_name, qObjName, launchesA, launchesB, vec_perps = cur_port
                lePorts += [(port_name + 'a', launchesA)]
                lePorts += [(port_name + 'b', launchesB)]

            #prepare design by converting shapely geometries to Gmsh geometries
            gmsh_render_attrs = pgr._prepare_design(metallic_layers, ground_plane, lePorts, self.user_options['fillet_resolution'], 'eigenmode_simulation')

            if self.create_files == True:
                #create directory to store simulation files
                self._create_directory(self.name)

                #create config file
                self.create_config_file(gmsh_render_attrs = gmsh_render_attrs)

                # #create batch file
                # if self.mode == 'HPC':
                #     self.create_batch_file()
                
                #create mesh
                pgr._intelligent_mesh('eigenmode_simulation', 
                                min_size = self.user_options['mesh_min'], 
                                max_size = self.user_options['mesh_max'], 
                                mesh_sampling = self.user_options['mesh_sampling'])

                self._save_mesh_gmsh()

            if self.view_design_gmsh_gui == True:
                #plot design in gmsh gui
                pgr.view_design_components()
                    
        if self.meshing == 'COMSOL':

            #initialise the COMSOL engine
            COMSOL_Model.init_engine()
            cmsl = COMSOL_Model('res')

            #Create COMSOL RF sim object
            sim_sParams = COMSOL_Simulation_RFsParameters(cmsl, adaptive = 'None')
            cmsl.initialize_model(self.metal_design, [sim_sParams], bottom_grounded = True, resolution = 10)

            #Add metallic layers
            for cur_layer in metallic_layers:
                layer_id = cur_layer.get('layer_id')
                cmsl.add_metallic(layer_id, **cur_layer)
            if not ground_plane['omit']:
                cmsl.add_ground_plane(threshold=ground_plane['threshold'])
            #cmsl.fuse_all_metals()

            #assign ports
            for cur_port in self._ports:
                port_name, qObjName, launchesA, launchesB, vec_perps = cur_port
                sim_sParams.create_port_CPW_on_Launcher(qObjName)
            
            #build model
            cmsl.build_geom_mater_elec_mesh(mesh_structure = self.user_options["comsol_meshing"])

            #plot model
            cmsl.plot()

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

    def create_port_CPW_on_Launcher(self, qObjName, len_launch = 20e-6):
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
        
        launchesA, launchesB, vec_perp = QUtilities.get_RFport_CPW_coords_Launcher(self.metal_design, qObjName, len_launch, 1e3)  #Units of mm...

        #Each port is defined as: [port_name, portA-coords, portB-coords, vec_CPW2GND_1] where vec_CPW2GND_1 is a vector pointing in the direction
        #of ground from the CPW for portA (given as +X, -X, +Y or -Y for AWS Palace...).
        #See here for details: https://awslabs.github.io/palace/stable/config/boundaries/#boundaries[%22LumpedPort%22]
        self._ports += [(port_name, qObjName, launchesA + [launchesA[-1]], launchesB + [launchesB[-1]], self._check_port_orientation(vec_perp))]

    def _check_port_orientation(self, vec_perp):
        thresh = 0.9999

        xDot = np.dot(vec_perp, [1,0])
        if xDot > thresh:
            return '+X', '-X'
        elif xDot < -thresh:
            return '-X', '+X'
        
        yDot = np.dot(vec_perp, [0,1])
        if yDot > thresh:
            return '+Y', '-Y'
        elif yDot < -thresh:
            return '-Y', '+Y'
        
        assert False, f"AWS Palace requires RF Lumped Ports to be aligned with the x/y axes. Here the port is pointing: {vec_perp}."

