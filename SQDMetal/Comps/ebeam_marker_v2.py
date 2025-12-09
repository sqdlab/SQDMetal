# -*- coding: utf-8 -*-

# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2021.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

# This class was created by Figen YILMAZ

# Imports required for drawing

# import numpy as np # (currently not used, may be needed later for component customization)
from qiskit_metal import draw
from qiskit_metal.draw.basic import subtract
from qiskit_metal.toolbox_python.attr_dict import Dict
from qiskit_metal.qlibrary.core import QComponent

# Define class and options for the markers geometry

class Markers(QComponent):
    """Markers for the EBPG marker search; positive markers.

    Inherits 'QComponent' class.

    .. image:
        Markers.png

    .. meta::
        Markers

    Creates 4 sets of markers for EBPG marker search. 
    The (0,0) point is the centre of these 4 markers.
    For how their named one can check the drawing below. 
    2 Dimensional Cartesian Coordinate System is the way how they are named.

    Sketch:
        Below is a sketch of one the marker set.
        ::
         __________________________________________________     y
        |                                                  |    ^                            |
        |      __       __       __       __       __      |    |                        2   |   1
        |     |__|     |__|     |__|     |__|     |__|     |    |                      ______|______ 
        |      1        2        3        4        5       |    |-------> x                  |
        |__________________________________________________|                             3   |   4
                                                                                             |   
        

    .. image::
        Markers.png

    Default Options:
        * pos_x = '0.0mm' -- the x position of the marker
        * pos_y = '0.0mm' -- the y position of the marker
        * marker_sep = '100um' -- the distance between each marker. Better to have at least 100um (Ref: Ivan)
        * marker_w = '20um' -- marker width. Better way to have succesful marker search (Ref: Ivan)
        * marker_h = '20um' -- marker height
        * markers_gap = '200um' -- Size of the pocket (cut out in ground)
    """

    default_options = Dict(
        pos_x = '0.0mm',
        pos_y = '0.0mm',
        marker_sep = '100um',
        marker_w = '20um',
        marker_h = '20um',
        markers_gap = '200um',
         )
    """Default options"""

    TOOLTIP = """Markers for the chip"""

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""

        p = self.p
        pos_x = p.pos_x
        pos_y = p.pos_y
        marker_sep = p.marker_sep
        marker_w = p.marker_w
        marker_h = p.marker_h
        markers_gap = p.markers_gap
        
        # Square shape markers for E-Beam marker search, see the drawing above for their positions
        # on the cartesian coordinate system
        marker_1 = draw.rectangle(marker_w, marker_h, pos_x/2 - marker_sep*2 - marker_w*2, pos_y/2)
        marker_2 = draw.rectangle(marker_w, marker_h, pos_x/2 - marker_sep - marker_w, pos_y/2)
        marker_3 = draw.rectangle(marker_w, marker_h, pos_x/2, pos_y/2)
        marker_4 = draw.rectangle(marker_w, marker_h, pos_x/2 + marker_sep + marker_w, pos_y/2)
        marker_5 = draw.rectangle(marker_w, marker_h, pos_x/2 + marker_sep*2 + marker_w*2, pos_y/2)

        markers = draw.union(marker_1, marker_2, marker_3, marker_4, marker_5)

        # Create the pocket for the markers to distinguish as positive markers on the substrate
        markers_pk = draw.rectangle(markers_gap*4, markers_gap*1.5, pos_x/2, pos_y/2)

        # Create polygon object list
        polys = [markers, markers_pk]
        # Rotates and translates all the objects as requested. Uses package functions
        # in 'draw_utility' for easy rotation/translation
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, xoff=pos_x/2, yoff=pos_y/2)
        [markers, markers_pk] = polys

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly', dict(markers=markers),layer=p.layer)
        # Subtracts out ground plane on the layer its on
        self.add_qgeometry('poly',
                           dict(markers_pk=markers_pk),
                           subtract=True,
                           layer=p.layer)