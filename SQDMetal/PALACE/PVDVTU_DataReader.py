import pyvista as pv
import numpy as np

class PVDVTU_DataReader:
    def __init__(self, pvd_file_name):
        self._reader = pv.get_reader(pvd_file_name)
        self.NumDatasets = len(self._reader.time_values) - 1    #It sets the last one to 99.0 and has no data; perhaps a terminating point?
        self._units = 0.001     #TODO: See if this should read/infer from config file?
    
    def get_fields(self, dataset_index = 0, field_type='E'):
        assert dataset_index >=0 and dataset_index < self.NumDatasets, f"Invalid zero-indexed index. There are only {self.NumDatasets} datasets"
        self._reader.set_active_time_value(self._reader.time_values[dataset_index]) # Load first timestep
        multi_block = self._reader.read()
        mesh = multi_block[0]   #Assuming a single mesh...
        if field_type == 'E':
            return np.array(mesh.point_data['E_real'] + 1j*mesh.point_data['E_imag'])
        elif field_type == 'B':
            return np.array(mesh.point_data['B_real'] + 1j*mesh.point_data['B_imag'])
        else:
            assert False, "Argument 'field_type' must be E or B."

    def get_mesh_coords(self, dataset_index = 0):
        assert dataset_index >=0 and dataset_index < self.NumDatasets, f"Invalid zero-indexed for 'dataset_index'. There are only {self.NumDatasets} datasets"
        self._reader.set_active_time_value(self._reader.time_values[dataset_index]) # Load first timestep
        multi_block = self._reader.read()
        mesh = multi_block[0]   #Assuming a single mesh...
        return np.array(mesh.points) * self._units

    def get_interpolated_fields(self, data_set_index, coords:np.ndarray, field_type='E'):
        assert len(coords.shape) == 2 and coords.shape[1] == 3, "The argument 'coords' must be an array of (x,y,z) triplets."
        multi_block = self._reader.read()
        mesh = multi_block[0]   #Assuming a single mesh...

        point_cloud = pv.PolyData(coords / self._units)
        sampled_mesh = point_cloud.sample(mesh)

        if field_type == 'E':
            return np.array(sampled_mesh['E_real'] + 1j*sampled_mesh['E_imag'])
        elif field_type == 'B':
            return np.array(sampled_mesh['B_real'] + 1j*sampled_mesh['B_imag'])
        else:
            assert False, "Argument 'field_type' must be E or B."
