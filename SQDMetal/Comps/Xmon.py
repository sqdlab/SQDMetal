# -*- coding: utf-8 -*-

# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0
# Author: Prasanna Pakkiam
# Creation Date: 29/05/2023
# Description: Collection of classes to draw Xmon qubits.

from qiskit_metal import draw
from qiskit_metal.toolbox_python.attr_dict import Dict
from qiskit_metal.qlibrary.core import QComponent
import shapely

class Xmon(QComponent):
    """Create an Xmon cross

    Inherits QComponent class.

    The Xmon consists of the cross and allocation to remove the surrounding ground plane.

    Dolan Bridge Josephson Junction Metal Geometry (no Ground Cutout):
        * cross_width  - Overall horizontal width of the cross
        * cross_height - Overall vertical height of the cross
        * vBar_width - Width of the horizontal bar of the cross
        * hBar_width - Width of the vertical bar of the cross
        * vBar_gap  - Gap to leave to the ground plane on the left and right of the vertical bar of the cross
        * hBar_gap  - Gap to leave to the ground plane on the top and bottom of the horizontal bar of the cross
        * gap_up  - Gap to extend to the ground plane on the top of the cross
        * gap_down  - Gap to extend to the ground plane on the bottom of the cross
        * gap_left  - Gap to extend to the ground plane on the left hand side of the cross
        * gap_right  - Gap to extend to the ground plane on the right hand side of the cross

    The positioning can be done via position (pos_x,pos_y) and angle given by 'orientation'.

    Pins:
        There are pins up, down, left and right to link to the different edges of the Xmon cross.

    Sketch:
        Below is a sketch of the Xmon along with its ground cut-out
        ::

            <.......CW.......>
          @@@@@@@@@@@@@@@@@@@@@@
          @@@@@@@  u__   @@@@@@@         @  = Ground Plane
          @@@@@@@  |  |  @@@@@@@   /|\   CW = cross_width
          @        |  |        @    |    CH = cross_height
          @l ______|  |______ r@    |    V  = vBar_width
          @ |                |H@   CH    H  = hBar_width
          @ |______    ______|H@    |    v  = vBar_gap
          @        |  |  h     @    |    h  = hBar_gap
          @        |  |  h     @    |    u  = gap_up
          @@@@@@@  |__|vv@@@@@@@   \|/   d  = gap_down
          @@@@@@@ d VV   @@@@@@@         l  = gap_left
          @@@@@@@@@@@@@@@@@@@@@@         r  = gap_right

    .. image::
        Cap3Interdigital.png

    .. meta::
        Xmon cross

    Default Options:
        * pos_x='0um',pos_y='0um'
        * orientation=0
        * cross_width='120um'
        * cross_height='100um'
        * vBar_width='20um'
        * hBar_width='20um'
        * vBar_gap='20um'
        * hBar_gap='20um'
        * gap_up='2um'
        * gap_down='2um'
        * gap_left='2um'
        * gap_right='2um'
    """

    default_options = Dict(pos_x='0um',pos_y='0um',
                           orientation=0,
                           cross_width='120um',
                           cross_height='100um',
                           vBar_width='20um',
                           hBar_width='20um',
                           vBar_gap='20um',
                           hBar_gap='20um',
                           gap_up='2um',
                           gap_down='2um',
                           gap_left='2um',
                           gap_right='2um',
                           )

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################

        cross = [
                (p.cross_width*0.5, p.hBar_width*0.5),
                (p.vBar_width*0.5, p.hBar_width*0.5),
                (p.vBar_width*0.5, p.cross_height*0.5),
                (-p.vBar_width*0.5, p.cross_height*0.5),
                (-p.vBar_width*0.5, p.hBar_width*0.5),
                (-p.cross_width*0.5, p.hBar_width*0.5),
                (-p.cross_width*0.5, -p.hBar_width*0.5),
                (-p.vBar_width*0.5, -p.hBar_width*0.5),
                (-p.vBar_width*0.5, -p.cross_height*0.5),
                (p.vBar_width*0.5, -p.cross_height*0.5),
                (p.vBar_width*0.5, -p.hBar_width*0.5),
                (p.cross_width*0.5, -p.hBar_width*0.5)
        ]

        gap = [
                (p.cross_width*0.5+p.gap_right, p.hBar_width*0.5+p.hBar_gap),
                (p.vBar_width*0.5+p.vBar_gap, p.hBar_width*0.5+p.hBar_gap),
                (p.vBar_width*0.5+p.vBar_gap, p.cross_height*0.5+p.gap_up),
                (-p.vBar_width*0.5-p.vBar_gap, p.cross_height*0.5+p.gap_up),
                (-p.vBar_width*0.5-p.vBar_gap, p.hBar_width*0.5+p.hBar_gap),
                (-p.cross_width*0.5-p.gap_left, p.hBar_width*0.5+p.hBar_gap),
                (-p.cross_width*0.5-p.gap_left, -p.hBar_width*0.5-p.hBar_gap),
                (-p.vBar_width*0.5-p.vBar_gap, -p.hBar_width*0.5-p.hBar_gap),
                (-p.vBar_width*0.5-p.vBar_gap, -p.cross_height*0.5-p.gap_down),
                (p.vBar_width*0.5+p.vBar_gap, -p.cross_height*0.5-p.gap_down),
                (p.vBar_width*0.5+p.vBar_gap, -p.hBar_width*0.5-p.hBar_gap),
                (p.cross_width*0.5+p.gap_right, -p.hBar_width*0.5-p.hBar_gap)
        ]

        cross = shapely.Polygon(cross)
        gap = shapely.Polygon(gap)
        
        pin_up = shapely.LineString([[p.vBar_width*0.5, p.cross_height*0.5], [-p.vBar_width*0.5, p.cross_height*0.5]])
        pin_down = shapely.LineString([[-p.vBar_width*0.5, -p.cross_height*0.5], [p.vBar_width*0.5, -p.cross_height*0.5]])
        pin_left = shapely.LineString([[-p.cross_width*0.5, p.hBar_width*0.5], [-p.cross_width*0.5, -p.hBar_width*0.5]])
        pin_right = shapely.LineString([[p.cross_width*0.5, -p.hBar_width*0.5], [p.cross_width*0.5, p.hBar_width*0.5]])

        polys = [cross, gap, pin_up, pin_down, pin_left, pin_right]
        polys = draw.rotate(polys, p.orientation, origin=(0, 0), use_radians=False)
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [cross, gap, pin_up, pin_down, pin_left, pin_right] = polys

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(cross=cross),
                           layer=p.layer)

        #subtracts out ground plane on the layer it's on
        self.add_qgeometry('poly',
                           dict(gap=gap),
                           subtract=True,
                           layer=p.layer)

        # Generates its own pins
        self.add_pin('up', pin_up.coords[::-1], width=p.vBar_width)
        self.add_pin('down', pin_down.coords[::-1], width=p.vBar_width)
        self.add_pin('left', pin_left.coords[::-1], width=p.hBar_width)
        self.add_pin('right', pin_right.coords[::-1], width=p.hBar_width)
