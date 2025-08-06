import pyvista as pv
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

class PVDVTU_Slice:
    def __init__(self, pv_slice, units):
        self.slc = pv_slice
        self._units = units
    
    def get_params(self):
        return self.slc.array_names
    
    def get_data(self, param_name):
        return self.slc[param_name]
    
    def plot(self, values, cmap='viridis', equalise_histogram=False, fig=None, ax=None):
        #Values must be a N-array - basically use get_data() and mathematically manipulate it to be a scalar array

        zMin, zMax = np.min(values), np.max(values)
        if equalise_histogram:
            quantileStops = np.linspace(0,1,64)
            leQuantiles = np.quantile(values, quantileStops)
            zinterpInv = lambda x: np.interp(x, quantileStops, leQuantiles) # noqa: E731
            zinterp = lambda x: np.interp(x, leQuantiles, quantileStops) # noqa: E731
            cMapNorm = matplotlib.colors.FuncNorm((zinterp,zinterpInv), vmin=zMin, vmax=zMax)
            cBarTicks = np.quantile(values, np.linspace(0,1,5))
        else:
            cMapNorm = None
            cBarTicks = np.linspace(zMin, zMax, 3)

        pts = self.slc.points
        tri = self.slc.faces.reshape((-1,4))[:, 1:]
        if fig is None:
            fig, ax = plt.subplots(1)
        pcm = ax.tripcolor(pts[:,0]*self._units, pts[:,1]*self._units, tri, values, norm=cMapNorm, cmap=cmap)
        ax.set_xlabel('x position (m)')
        ax.set_ylabel('y position (m)')
        self._colbar_obj = fig.colorbar(pcm, shrink=0.6, ticks=cBarTicks)
        return fig
    
    def plot_mesh(self):
        pts = self.slc.points
        tri = self.slc.faces.reshape((-1,4))[:, 1:]
        fig, ax = plt.subplots(1)
        ax.triplot(pts[:,0]*self._units, pts[:,1]*self._units, tri, linewidth=0.1, color='red')
        ax.set_xlabel('x position (m)')
        ax.set_ylabel('y position (m)')
        return fig

class PVDVTU_Viewer:
    def __init__(self, pvd_file, units=0.001):
        self._reader = pv.get_reader(pvd_file)
        self.num_datasets = len(self._reader.datasets)
        self._units = units
    
    def get_data_slice(self, data_set_index, slice_plane_normal=np.array([0,0,1]), slice_plane_origin=np.array([0,0,0])):
        assert data_set_index < self.num_datasets, f"Argument \'data_set_index\' must be in: [0,{self.num_datasets-1}]."
        self._reader.set_active_time_value(self._reader.time_values[data_set_index])
        mesh = self._reader.read()
        return PVDVTU_Slice(mesh[0].slice(normal=slice_plane_normal, origin=slice_plane_origin, generate_triangles=True,  contour=False), units=self._units)
