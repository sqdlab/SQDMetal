# -*- coding: utf-8 -*-
# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0
# Description: Components to mark wire-bond locations in a design. The footprint is drawn (for visualisation and
#              GDS export) on a dedicated layer, while the 3D geometry is generated in the PALACE simulation classes
#              via PALACE_Model_Base.add_wirebonds_from_design().

from qiskit_metal.toolbox_python.attr_dict import Dict
from qiskit_metal.qlibrary.core import QComponent
import numpy as np
import shapely
from SQDMetal.Utilities.QiskitShapelyRenderer import QiskitShapelyRenderer
from SQDMetal.Utilities.QUtilities import QUtilities


class _WireBondSimOptions:
    """
    Mixin handling the simulation options that are passed onto the PALACE wire-bond builder. Any option left as an
    empty string '' falls back onto the simulation's wire-bond defaults (see PALACE_Model_Base.set_wirebond_defaults).
    """
    sim_length_options = ['width', 'thickness', 'height', 'foot_length', 'foot_width', 'foot_thickness']
    sim_string_options = ['cross_section', 'material', 'profile', 'bond_type']
    sim_default_options = Dict(width='', thickness='', height='', foot_length='', foot_width='', foot_thickness='',
                               cross_section='', material='', profile='', bond_type='')

    def _get_sim_options(self):
        unit_conv = QUtilities.get_units(self.design)
        ret = {}
        for k in self.sim_length_options:
            if self.options[k] != '':
                ret[k] = float(self.p[k]) * unit_conv
        for k in self.sim_string_options:
            if self.options[k] != '':
                ret[k] = str(self.options[k])
        return ret

    def _footprint(self, pA, pB, width, foot_len):
        line = shapely.LineString([pA, pB]).buffer(width/2, cap_style='flat')
        u = (np.array(pB) - np.array(pA)) / np.linalg.norm(np.array(pB) - np.array(pA))
        feet = [shapely.LineString([p - u*foot_len/2, p + u*foot_len/2]).buffer(width, cap_style='flat') for p in [np.array(pA), np.array(pB)]]
        return shapely.unary_union([line] + feet)


class WireBond(_WireBondSimOptions, QComponent):
    """A single wire-bond between two points (start and end feet). The feet may land on the ground plane, a CPW
    trace (or any other conductor), or lie off-chip (e.g. a bond to the PCB/sample-box). The type is detected
    automatically when building the simulation.

    Inherits QComponent class.

    Placement options:
        * start_x, start_y - Position of the first foot (ball/first-wedge bond)
        * end_x, end_y     - Position of the second foot (wedge bond)
        * layer            - Layer on which to draw the plan-view footprint (purely for visualisation/GDS export;
                             do not pass this layer to add_metallic)
        * draw_footprint   - If True (default), draws the plan-view footprint.

    Simulation options (leave as '' to use the simulation defaults):
        * width, thickness, height, foot_length, foot_width, foot_thickness - lengths
        * cross_section - 'round', 'ribbon' or 'sheet'
        * material      - e.g. 'Al', 'Au', 'Cu'
        * profile       - 'semicircle', 'arc', 'elliptical', 'parabolic', 'trapezoid', 'ball_wedge'
        * bond_type     - '' or 'auto' (default) or an expected type ('ground', 'cpw_trace', 'ground_stitch', ...) which
                          is checked against the automatically detected type.

    This class ignores pos_x, pos_y and orientation.

    Pins:
        There are no pins in this object.

    .. meta::
        Wire bond
    """

    default_options = Dict(start_x='0um', start_y='0um', end_x='100um', end_y='0um', layer=50, draw_footprint=True,
                           **_WireBondSimOptions.sim_default_options)

    def make(self):
        p = self.p
        pA = np.array([p.start_x, p.start_y])
        pB = np.array([p.end_x, p.end_y])
        assert np.linalg.norm(pB - pA) > 0, f"WireBond {self.name}: start and end positions coincide."
        self.wirebond_feet = [(pA, pB)]
        if p.draw_footprint:
            unit_conv = QUtilities.get_units(self.design)
            w = float(p.width) if self.options.width != '' else 25e-6/unit_conv
            fl = float(p.foot_length) if self.options.foot_length != '' else 2.5*w
            self.add_qgeometry('poly', dict(wirebond=self._footprint(pA, pB, w, fl)), layer=p.layer)

    def get_wirebond_specs(self):
        """Returns the wire-bond(s) of this component as a list of dictionaries in metres."""
        unit_conv = QUtilities.get_units(self.design)
        p = self.p
        return [dict(name=f"wb_{self.name}",
                     pos_start=np.array([p.start_x, p.start_y])*unit_conv,
                     pos_end=np.array([p.end_x, p.end_y])*unit_conv,
                     **self._get_sim_options())]


class WireBondsCPW(_WireBondSimOptions, QComponent):
    """Places ground-to-ground wire-bonds (stitches) across a CPW at regular intervals - analogous to AirBridgeCPW. The
    bonds are perpendicular to the CPW with feet landing on the ground plane on either side.

    Inherits QComponent class.

    Placement options:
        * pathObj            - Name of the CPW component (e.g. RouteMeander) to stitch across.
        * dist_wirebonds     - Distance between consecutive wire-bonds along the CPW.
        * dist_wirebonds_ends - Padding at the start and end of the CPW in which no wire-bonds are placed.
        * foot_offset        - Distance from the outer edge of the CPW gap to the centre of each foot.
        * layers_obj_avoid   - Layers in which obstacles are avoided (a bond is skipped if a foot overlaps any
                               non-ground metal/gap in these layers).
        * layer              - Layer on which to draw the footprints (visualisation only).
        * draw_footprint     - If True (default), draws the plan-view footprints.

    Simulation options: as in WireBond.

    The positions are given in the attributes current_positions and current_arclengths.

    Sketch:
        ::

            ===O========O========O===  ground
            ---|--------|--------|---  gap
            ###|########|########|###  trace
            ---|--------|--------|---  gap
            ===O========O========O===  ground      O = feet, | = bond

    .. meta::
        Wire-bond stitches across a CPW
    """

    default_options = Dict(pathObj='', dist_wirebonds='300um', dist_wirebonds_ends='50um', foot_offset='40um',
                           layers_obj_avoid=[1, 2, 3], layer=50, draw_footprint=True,
                           **_WireBondSimOptions.sim_default_options)

    def make(self):
        p = self.p
        assert p.pathObj != '', "Must provide a valid component name for the option: pathObj."

        final_pts, normals, width, gap, total_dist = QUtilities.calc_points_on_path([], self._design, self.options.pathObj)
        dists = np.arange(p.dist_wirebonds_ends, total_dist - p.dist_wirebonds_ends, p.dist_wirebonds)
        if len(dists) == 0:
            dists = np.array([total_dist/2])
        final_pts, normals, width, gap, total_dist = QUtilities.calc_points_on_path(dists, self._design, self.options.pathObj)
        assert gap > 0, f"The component {self.options.pathObj} does not appear to be a valid CPW with a finite trace width and finite trace gap."

        unit_conv = QUtilities.get_units(self.design)
        w = float(p.width) if self.options.width != '' else 25e-6/unit_conv
        fl = float(p.foot_length) if self.options.foot_length != '' else 2.5*w
        fw = float(p.foot_width) if self.options.foot_width != '' else 1.5*w

        qmpl = QiskitShapelyRenderer(design=self.design, canvas=None, logger=None)
        gsdf = qmpl.get_net_coordinates()
        gsdf = gsdf[gsdf['layer'].isin(p.layers_obj_avoid)]
        gsdf = gsdf[gsdf['component'] != self.id] if 'component' in gsdf else gsdf
        obstacles = shapely.unary_union(gsdf['geometry'])

        half_span = width*0.5 + gap + p.foot_offset
        assert p.foot_offset > fl/2, f"WireBondsCPW {self.name}: foot_offset ({p.foot_offset}) must exceed half the foot length ({fl/2}) so that the feet clear the CPW gap."
        self.wirebond_feet = []
        self.current_positions = []
        self.current_arclengths = []
        footprints = []
        for m in range(final_pts.shape[0]):
            n = np.array(normals[m])
            pA = final_pts[m] - n*half_span
            pB = final_pts[m] + n*half_span
            feet = shapely.MultiPolygon([shapely.LineString([x - n*fl/2, x + n*fl/2]).buffer(fw/2, cap_style='flat') for x in [pA, pB]])
            if obstacles.is_empty or shapely.intersection(feet, obstacles).area / feet.area < 1e-6:
                self.wirebond_feet.append((pA, pB))
                self.current_positions.append(final_pts[m])
                self.current_arclengths.append(dists[m])
                footprints.append(self._footprint(pA, pB, w, fl))
        self.current_positions = np.array(self.current_positions)
        self.current_arclengths = np.array(self.current_arclengths)

        if p.draw_footprint and len(footprints) > 0:
            self.add_qgeometry('poly', dict(wirebonds=shapely.unary_union(footprints)), layer=p.layer)

    def get_wirebond_specs(self):
        unit_conv = QUtilities.get_units(self.design)
        sim_opts = self._get_sim_options()
        if 'bond_type' not in sim_opts:
            sim_opts['bond_type'] = 'ground'
        return [dict(name=f"wb_{self.name}_{m}", pos_start=pA*unit_conv, pos_end=pB*unit_conv, **sim_opts)
                for m, (pA, pB) in enumerate(self.wirebond_feet)]
