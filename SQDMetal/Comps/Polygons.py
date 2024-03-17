# -*- coding: utf-8 -*-
# Author: Prasanna Pakkiam
# Creation Date: 13/03/2024
# Description: Collection of classes to draw simple polygons. These do NOT replace markers etc. These are general
#              purpose classes designed for niche applications.

from qiskit_metal import draw
from qiskit_metal.toolbox_python.attr_dict import Dict
from qiskit_metal.qlibrary.core import QComponent
import shapely

class PolyBorderRectangle(QComponent):
    """Create a rectangle with a border of some thickness. Can either be solid metal or just a ground-plane cut-out.

    Inherits QComponent class.

    Square marker either has a Metal or ground cutout Geometry as specified by is_ground_cutout:
        * width  - Inner width of square along x-axis
        * height - Inner height of square along y-axis
        * border_thickness_X - Thickness of the border outside the specified rectangle along the x-axis
        * border_thickness_Y - Thickness of the border outside the specified rectangle along the y-axis
        * is_ground_cutout - If True, the border is a ground cutout, otherwise it is a metallic square.

    As usual, the centre-positioning can be done dynamically via (pos_x,pos_y) with orientation specifying an optional rotation.
        
    Pins:
        There are pins N, E, W and S on the rectangle on the centre across every edges of the outer border
        There are pins n, e, w and s on the rectangle on the centre across every edges of the inner border

    Sketch:
        Below is a sketch of the solitary square marker
        ::

                <-------L------->
            @@@@@@@@@@@@N@@@@@@@@@@@@
            @@@@@@@@@@@@n@@@@@@@@@@@@              @ = ground-plane cutout or metal depending on is_ground_cutout
            @@@@/|\              UUUU
            @@@@ |               @@@@              x = (pos_x, pos_y)
            W@@w h      x        e@@E              U,V = border_thickness_X, border_thickness_Y
            @@@@ |               @@@@              h = height
            @@@@\|/              @@@@              L = width
            @@@@@@@@@@@@s@@@@@V@@@@@@              NEWS are pins on the outer rectangle
            @@@@@@@@@@@@S@@@@@V@@@@@@              news are pins on the inner rectangle

    .. image::
        Cap3Interdigital.png

    .. meta::
        Solitary square marker

    Default Options:
        * pos_x='0um',pos_y='0um'
        * width='50um'
        * height='50um'
        * border_thickness_X='20um'
        * border_thickness_Y='20um'
        * is_ground_cutout=False
    """

    default_options = Dict(pos_x='0um',pos_y='0um', orientation=0,
                           width='50um',
                           height='50um',
                           border_thickness_X='20um',
                           border_thickness_Y='20um',
                           is_ground_cutout=False,
                           cpw_width='10um')

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################

        rectangle_inner = [
                (+p.width*0.5, +p.height*0.5),
                (-p.width*0.5, +p.height*0.5),
                (-p.width*0.5, -p.height*0.5),
                (+p.width*0.5, -p.height*0.5),
        ]

        rectangle_outer = [
                (+p.width*0.5+p.border_thickness_X, +p.height*0.5+p.border_thickness_Y),
                (-p.width*0.5-p.border_thickness_X, +p.height*0.5+p.border_thickness_Y),
                (-p.width*0.5-p.border_thickness_X, -p.height*0.5-p.border_thickness_Y),
                (+p.width*0.5+p.border_thickness_X, -p.height*0.5-p.border_thickness_Y),
        ]

        pinN = shapely.LineString(rectangle_outer[0:2])
        pinS = shapely.LineString(rectangle_outer[2:])
        pinE = shapely.LineString([rectangle_outer[-1], rectangle_outer[0]])
        pinW = shapely.LineString(rectangle_outer[1:3])

        pinn = shapely.LineString(rectangle_inner[0:2])
        pins = shapely.LineString(rectangle_inner[2:])
        pine = shapely.LineString([rectangle_inner[-1], rectangle_inner[0]])
        pinw = shapely.LineString(rectangle_inner[1:3])

        rectangle_inner = shapely.Polygon(rectangle_inner)
        rectangle_outer = shapely.Polygon(rectangle_outer)
        rectangle = rectangle_outer.difference(rectangle_inner)

        polys = [rectangle, pinN, pinS, pinE, pinW, pinn, pins, pine, pinw]
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        rectangle, pinN, pinS, pinE, pinW, pinn, pins, pine, pinw = polys

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(rectangle=rectangle),
                           layer=p.layer,
                           subtract=p.is_ground_cutout)
        self.add_pin('N', pinN.coords[::-1], width=p.cpw_width)
        self.add_pin('S', pinS.coords[::-1], width=p.cpw_width)
        self.add_pin('E', pinE.coords[::-1], width=p.cpw_width)
        self.add_pin('W', pinW.coords[::-1], width=p.cpw_width)
        #
        self.add_pin('n', pinn.coords[::-1], width=p.cpw_width)
        self.add_pin('s', pins.coords[::-1], width=p.cpw_width)
        self.add_pin('e', pine.coords[::-1], width=p.cpw_width)
        self.add_pin('w', pinw.coords[::-1], width=p.cpw_width)
