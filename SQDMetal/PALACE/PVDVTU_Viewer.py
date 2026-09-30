# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0

import pyvista as pv
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from SQDMetal.Utilities.PlotUtilities import HistogramEqualisation

class PVDVTU_Slice:
    def __init__(self, pv_slice, units):
        self.slc = pv_slice
        self._units = units
    
    def get_params(self):
        return self.slc.array_names
    
    def get_data(self, param_name):
        return self.slc[param_name]
    
    def plot(self, values, cmap='viridis', equalise_histogram=True, equal_pos_neg=True, fig=None, ax=None):
        #Values must be a N-array - basically use get_data() and mathematically manipulate it to be a scalar array

        zMin, zMax = np.min(values), np.max(values)
        if equal_pos_neg:
            zMag = max(np.abs(zMin),np.abs(zMax))
        if equalise_histogram:
            histeq = HistogramEqualisation(values, equal_pos_neg)
            cMapNorm = histeq.cmap()
            cBarTicks = histeq.cTicks(5)
        else:
            cMapNorm = None
            if equal_pos_neg:
                cBarTicks = np.linspace(-zMag, zMag, 3)
            else:
                cBarTicks = np.linspace(zMin, zMax, 3)

        pts = self.slc.points
        tri = self.slc.faces.reshape((-1,4))[:, 1:]
        if fig is None:
            fig, ax = plt.subplots(1)
        if equal_pos_neg and not equalise_histogram:
            pcm = ax.tripcolor(pts[:,0]*self._units, pts[:,1]*self._units, tri, values, norm=cMapNorm, cmap=cmap, vmin=-zMag, vmax=zMag)
        else:
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
