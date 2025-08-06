# -*- coding: utf-8 -*-
# Author: Prasanna Pakkiam
# Creation Date: 02/08/2023
# Description: Collection of classes to draw air bridges.

from qiskit_metal.toolbox_python.attr_dict import Dict
from qiskit_metal.qlibrary.core import QComponent
import numpy as np
import shapely
from SQDMetal.Utilities.QiskitShapelyRenderer import QiskitShapelyRenderer
from SQDMetal.Utilities.QUtilities import QUtilities

class AirBridgeCPW(QComponent):
    """Places air bridges along a given CPW path. The positions of the air-bridges are given in handy object attributes:
        * current_positions  - xy-coordinates of the air-bridges (in the Qiskit-Metal units)
        * current_arclengths - Distance (along the CPW path - i.e. arclength) of each air-bridge.

    Inherits QComponent class.

    The air bridges are are comprised of the first pad layer and the second I-bar layer that forms the actual bridge. The basic
    options are specified via:
        * pathObj - Name of component in the design to which to attach air bridges.
        * layer   - Layer in which to place the first layer pads on either side of the CPW
        * layer_bridge     - Layer in which to place the I-bar in the second layer
        * dist_air_bridges - Distance between consecutive air bridges along the CPW.
        * dist_air_bridges_ends - Padding distances on the start and end of the CPW; air bridges are not placed in this region.
        * layers_obj_avoid - Layers in which obstacles (i.e. blue or grey regions in the GUI) are avoided. That is, if any non-
                             ground-plane metals/gaps intersect with an air bridge's pads, that air bridge will not be placed.
    The first layer pad parameters are specified via:
        * pad_gap    - Gap from the CPW gap to the air bridge pad.
        * pad_width  - Width of the air bridge pad
        * pad_length - Length of the air bridge pad along the CPW
    The second layer I-bar parameters are specified via:
        * bridge_pad_gap    - Gap from the first layer pad to the pad-sections on the I-bar.
        * bridge_pad_width  - Width of the pad-sections on the I-bar.
        * bridge_pad_length - Length of the pad-sections on the I-bar along the CPW
        * bridge_width      - The width of the stem in the I-bar (the actual bridge across the CPW)

    This class ignores pos_x, pos_y and orientation...
        
    Pins:
        There are no pins in this object.

    Sketch:
        Below is a sketch of a given air bridge on a CPW along with the available options.
        ::
            --------wwwwwwwwww--------
            --------------------------
            --------OOOOOOOOOO-l------
            ------Y-OOO####OOO-l------    - = CPW
            ------Y-OOO####OOO-l------
            ------G-OOOO##OOOO-l------    O = First layer pads
            ------------##---g--------    g = pad_gap
            ------------##---g--------    w = pad_width
                        ##                l = pad_length
            ------------##------------    
                        BB                # = Second layer I-bar
            ------------##---g--------    G = bridge_pad_gap
            ------------##---g--------    Y = bridge_pad_width
            --------OOOO##OOOO-G------    L = bridge_pad_length
            --------OOO####OOO-Y------    B = bridge_width
            --------OOO####OOO-Y------
            --------OOOOOOOOOO--------
            --------------------------
            -----------LLLL-----------

    .. image::
        Cap3Interdigital.png

    .. meta::
        Air bridges along a CPW.

    Default Options:
        * pathObj=''
        * layer_bridge=5
        * dist_air_bridges='100um'
        * dist_air_bridges_ends='50um'
        * layers_obj_avoid=[1,2,3]
        * pad_gap='5um'
        * pad_width='20um'
        * pad_length='40um'
        * bridge_pad_gap='5um'
        * bridge_pad_width='5um'
        * bridge_pad_length='15um'
        * bridge_width='8um
    """

    default_options = Dict(pathObj='',
                           layer_bridge=5,
                           dist_air_bridges='100um',
                           dist_air_bridges_ends='50um',
                           layers_obj_avoid=[1,2,3],
                           pad_gap='5um',
                           pad_width='20um',
                           pad_length='40um',
                           bridge_pad_gap='5um',
                           bridge_pad_width='5um',
                           bridge_pad_length='15um',
                           bridge_width='8um')

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################

        assert p.pathObj != '', "Must provide a valid component name for the option: pathObj."

        final_pts, normals, width, gap, total_dist = QUtilities.calc_points_on_path([], self._design, self.options.pathObj)
        dists = np.arange(p.dist_air_bridges_ends, total_dist-p.dist_air_bridges_ends, p.dist_air_bridges)
        if total_dist - (dists[-1] + p.dist_air_bridges) > p.dist_air_bridges_ends:
            dists = np.concatenate([dists, [dists[-1] + p.dist_air_bridges]])
        final_pts, normals, width, gap, total_dist = QUtilities.calc_points_on_path(dists, self._design, self.options.pathObj)
        assert gap > 0, f"The component {self.options.pathObj} does not appear to be a valid CPW with a finite trace width and finite trace gap."

        qmpl = QiskitShapelyRenderer(None, self.design, None)
        gsdf = qmpl.get_net_coordinates()
        gsdf = gsdf[gsdf['layer'].isin(p.layers_obj_avoid)]
        obstacles = shapely.unary_union(gsdf['geometry'])

        assert p.layer_bridge != p.layer, "The second air bridge layer must be different to the first pad layer."

        airbridge_pads = []
        airbridge_iBars = []
        self.current_positions = []
        self.current_arclengths = []
        for m in range(final_pts.shape[0]):
            matTrans = np.array([normals[m], [-normals[m][1], normals[m][0]]])

            #Draw Pads on the first layer
            cur_pts = np.array([
                (width*0.5 + gap + p.pad_gap, p.pad_length*0.5),
                (width*0.5 + gap + p.pad_gap, -p.pad_length*0.5),
                (width*0.5 + gap + p.pad_gap + p.pad_width, -p.pad_length*0.5),
                (width*0.5 + gap + p.pad_gap + p.pad_width, p.pad_length*0.5)])
            cur_pts = cur_pts@matTrans + final_pts[m]
            cur_pts2 = np.array([
                (-width*0.5 - gap - p.pad_gap, -p.pad_length*0.5),
                (-width*0.5 - gap - p.pad_gap, p.pad_length*0.5),
                (-width*0.5 - gap - p.pad_gap - p.pad_width, p.pad_length*0.5),
                (-width*0.5 - gap - p.pad_gap - p.pad_width, -p.pad_length*0.5)])
            cur_pts2 = cur_pts2@matTrans + final_pts[m]

            #Check that the air bridge pads do not intersect with any obstacles...
            thresh = 1e-6
            chk_geom = shapely.MultiPolygon([(cur_pts,[]), (cur_pts2,[])])
            if shapely.intersection(chk_geom, obstacles).area / chk_geom.area < thresh:
                airbridge_pads += [cur_pts, cur_pts2]
                self.current_positions += [final_pts[m]]
                self.current_arclengths += [dists[m]]
            else:
                continue
            
            #Draw the I-bar for the bridge on the second layer
            iBar = [
                (width*0.5 + gap + p.pad_gap + p.bridge_pad_gap, -p.bridge_width*0.5),
                (width*0.5 + gap + p.pad_gap + p.bridge_pad_gap, -p.bridge_pad_length*0.5),
                (width*0.5 + gap + p.pad_gap + p.bridge_pad_gap + p.bridge_pad_width, -p.bridge_pad_length*0.5),
                (width*0.5 + gap + p.pad_gap + p.bridge_pad_gap + p.bridge_pad_width, p.bridge_pad_length*0.5),
                (width*0.5 + gap + p.pad_gap + p.bridge_pad_gap, p.bridge_pad_length*0.5),
                (width*0.5 + gap + p.pad_gap + p.bridge_pad_gap, p.bridge_width*0.5),
                #
                (-width*0.5 - gap - p.pad_gap - p.bridge_pad_gap, p.bridge_width*0.5),
                (-width*0.5 - gap - p.pad_gap - p.bridge_pad_gap, p.bridge_pad_length*0.5),
                (-width*0.5 - gap - p.pad_gap - p.bridge_pad_gap - p.bridge_pad_width, p.bridge_pad_length*0.5),
                (-width*0.5 - gap - p.pad_gap - p.bridge_pad_gap - p.bridge_pad_width, -p.bridge_pad_length*0.5),
                (-width*0.5 - gap - p.pad_gap - p.bridge_pad_gap, -p.bridge_pad_length*0.5),
                (-width*0.5 - gap - p.pad_gap - p.bridge_pad_gap, -p.bridge_width*0.5)
            ]
            airbridge_iBars += [np.array(iBar)@matTrans + final_pts[m]]


        self.current_positions = np.array(self.current_positions)
        self.current_arclengths = np.array(self.current_arclengths)

        airbridge_pads = [(x,[]) for x in airbridge_pads]
        airbridge_iBars = [(x,[]) for x in airbridge_iBars]

        # Adds the objects to the qgeometry table
        self.add_qgeometry('poly',
                           dict(pads=shapely.MultiPolygon(airbridge_pads)),
                           layer=p.layer)
        self.add_qgeometry('poly',
                           dict(iBars=shapely.MultiPolygon(airbridge_iBars)),
                           layer=p.layer_bridge)
