# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0

import pyvista as pv
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from SQDMetal.Utilities.PlotUtilities import HistogramEqualisation

class PVDVTU_Slice:
    def __init__(self, pv_slice, units, normal, vec_x_axis=None):
        # Ensure we have triangular surface cells.
        self.slc = pv_slice.triangulate()
        self._units = units

        norm = np.linalg.norm(normal)
        assert norm > 0, "Slice-plane normal must be non-zero."

        self._normal = normal / norm

        if vec_x_axis is None:
            #Choose the Cartesian axis least aligned with the normal
            reference = np.eye(3)[np.argmin(np.abs(normal))]
        else:
            reference = np.array(vec_x_axis)*1.0
        self._v = np.cross(normal, reference)
        assert np.linalg.norm(self._v) > 0, "The supplied vec_x_axis is parallel with the normal..."
        self._v /= np.linalg.norm(self._v)
        self._u = np.cross(self._v, normal)
        self._u /= np.linalg.norm(self._u)

    def get_params(self):
        # Include both point- and cell-associated arrays.
        return sorted(set(self.slc.point_data.keys()) | set(self.slc.cell_data.keys()))

    def get_2d_coordinates(self):
        """
        Project the 3D slice onto its local 2D coordinate system.

        Returns coordinates in the same units as the original mesh.
        """
        pts = self.slc.points

        if len(pts) == 0:
            raise ValueError("The slice contains no points.")

        # Referencing the slice centre keeps the 2D coordinates compact.
        origin = pts.mean(axis=0)
        relative_pts = pts - origin

        x2d = relative_pts @ self._u
        y2d = relative_pts @ self._v

        return x2d, y2d

    def get_data(self, param_name, association="point"):
        """
        Retrieve a field with an explicit data association.

        association: 'point' or 'cell'
        """
        if association == "point":
            return self.slc.point_data[param_name]
        if association == "cell":
            return self.slc.cell_data[param_name]
        raise ValueError("association must be 'point' or 'cell'")

    def plot(self, values, cmap="viridis", equalise_histogram=True, equal_pos_neg=True, association="point", fig=None, ax=None):
        """
        Plot a scalar field on an arbitrarily oriented planar slice.

        association:
            'point' -> one scalar per mesh vertex
            'cell'  -> one scalar per triangle
        """

        values = np.asarray(values)

        # This plotting method expects a scalar field.
        if values.ndim != 1:
            raise ValueError("Expected a 1D scalar array. For vector fields, select a component or calculate a scalar quantity such as the magnitude first.")

        if association == "point":
            expected = self.slc.n_points
        elif association == "cell":
            expected = self.slc.n_cells
        else:
            raise ValueError("association must be 'point' or 'cell'")

        if len(values) != expected:
            raise ValueError(f"{association.capitalize()} data has {len(values)} values, but the slice has {expected} {'points' if association == 'point' else 'cells'}.")

        # Complex arrays must be converted to an explicitly chosen
        # real-valued quantity before plotting.
        if np.iscomplexobj(values):
            raise ValueError("Complex data must be converted to a real scalar, e.g. np.abs(values) or np.real(values).")

        values = values.astype(float, copy=False)

        if not np.all(np.isfinite(values)):
            raise ValueError("Field contains NaN or infinite values. Handle these explicitly before plotting.")

        pts = self.slc.points
        x2d, y2d = self.get_2d_coordinates()

        # triangulate() ensures triangular surface cells.
        faces = self.slc.faces

        if len(faces) % 4 != 0:
            raise ValueError("Expected triangular cells after triangulation.")

        faces = faces.reshape(-1, 4)

        if not np.all(faces[:, 0] == 3):
            raise ValueError("The slice contains non-triangular faces.")

        triangles = faces[:, 1:]

        if fig is None and ax is None:
            fig, ax = plt.subplots()
        elif ax is None:
            ax = fig.subplots()
        elif fig is None:
            fig = ax.figure

        # Determine colour limits.
        z_min = np.min(values)
        z_max = np.max(values)

        if equal_pos_neg:
            z_mag = max(abs(z_min), abs(z_max))
            z_min, z_max = -z_mag, z_mag

            # Avoid a degenerate colour range for an all-zero field.
            if z_mag == 0:
                z_min, z_max = -1.0, 1.0

        # Construct the colour normalisation.
        if equalise_histogram:
            histeq = HistogramEqualisation(values, equal_pos_neg)
            cmap_norm = histeq.cmap()
            cbar_ticks = histeq.cTicks(5)
        else:
            from matplotlib.colors import Normalize

            cmap_norm = Normalize(vmin=z_min, vmax=z_max)
            cbar_ticks = np.linspace(z_min, z_max, 5)

        # Avoid passing both norm and vmin/vmax to tripcolor.
        if association == "point":
            pcm = ax.tripcolor( x2d * self._units, y2d * self._units, triangles, values, shading="gouraud", cmap=cmap, norm=cmap_norm )
        else:
            pcm = ax.tripcolor( x2d * self._units, y2d * self._units, triangles, facecolors=values, shading="flat", cmap=cmap, norm=cmap_norm )

        ax.set_xlabel("In-plane coordinate u (m)")
        ax.set_ylabel("In-plane coordinate v (m)")
        ax.set_aspect("equal", adjustable="box")

        self._colbar_obj = fig.colorbar( pcm, ax=ax, shrink=0.6, ticks=cbar_ticks )

        return fig

    def plot_mesh(self, fig=None, ax=None):
        x2d, y2d = self.get_2d_coordinates()

        faces = self.slc.faces.reshape(-1, 4)
        triangles = faces[:, 1:]

        if fig is None and ax is None:
            fig, ax = plt.subplots()
        elif ax is None:
            ax = fig.subplots()
        elif fig is None:
            fig = ax.figure

        ax.triplot( x2d * self._units, y2d * self._units, triangles, linewidth=0.1, color="red" )

        ax.set_xlabel("In-plane coordinate u (m)")
        ax.set_ylabel("In-plane coordinate v (m)")
        ax.set_aspect("equal", adjustable="box")

        return fig

class PVDVTU_Viewer:
    def __init__(self, pvd_file, units=0.001):
        self._reader = pv.get_reader(pvd_file)
        self.num_datasets = len(self._reader.datasets)
        self._units = units
    
    def get_data_slice(self, data_set_index, slice_plane_normal=np.array([0,0,1]), slice_plane_origin=np.array([0,0,0]), vec_x_axis=None):
        assert data_set_index < self.num_datasets, f"Argument \'data_set_index\' must be in: [0,{self.num_datasets-1}]."
        self._reader.set_active_time_value(self._reader.time_values[data_set_index])
        mesh = self._reader.read()
        return PVDVTU_Slice(mesh[0].slice(normal=slice_plane_normal, origin=slice_plane_origin/self._units, generate_triangles=True,  contour=False), units=self._units, normal=slice_plane_normal, vec_x_axis=vec_x_axis)
