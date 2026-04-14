# -*- coding: utf-8 -*-

# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0
# Author: Prasanna Pakkiam
# Creation Date: 03/06/2023
# Description: Collection of classes to draw markers.

from qiskit_metal import draw
from qiskit_metal.toolbox_python.attr_dict import Dict
from qiskit_metal.qlibrary.core import QComponent
import numpy as np
import shapely


class MarkerDicingCross(QComponent):
    """Create a cross that can be used as a dicing marker.

    Inherits QComponent class.

    Cross Metal Geometry (with Ground Cutout):
        * cross_size - Total width and height of the metallic cross
        * bar_width  - Width of the horizontal and vertical bars forming the cross
        * cross_gap  - Size of the uniform gap from cross to the ground plane.

    As usual, the centre-positioning can be done dynamically via (pos_x,pos_y) with orientation specifying an optional rotation.
        
    Pins:
        There are no pins given that overlap and precise positioning is usually the concern...

    Sketch:
        Below is a sketch of the Cross Marker along with its ground cut-out
        ::

            <.......C.......>
          @@@@@@@@@@@@@@@@@@@@@
          @@@@@@@  G_   @@@@@@@         @ = Ground Plane
          @@@@@@@  | |  @@@@@@@   /|\   x = (pox_x, pos_y)
          @@@@@@@  |B|  @@@@@@@    |    C = cross_size
          @  ______| |______  @    |    B = bar_width
          @ |______ x ______|G@    C    G = cross_gap
          @        | |   G    @    |    
          @@@@@@@  | |  h@@@@@@    |    
          @@@@@@@  |_|G @@@@@@@   \|/   
          @@@@@@@   G   @@@@@@@         
          @@@@@@@@@@@@@@@@@@@@@         

    .. image::
        Cap3Interdigital.png

    .. meta::
        Dicing Marker Cross

    Default Options:
        * pos_x='0um',pos_y='0um',orientation=0,
        * cross_size='400um'
        * bar_width='40um'
        * cross_gap='50um'
    """

    default_options = Dict(pos_x='0um',pos_y='0um', orientation=0,
                           cross_size='400um',
                           bar_width='40um',
                           cross_gap='50um')

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################

        cross = [
                (p.cross_size*0.5, p.bar_width*0.5),
                (p.bar_width*0.5, p.bar_width*0.5),
                (p.bar_width*0.5, p.cross_size*0.5),
                (-p.bar_width*0.5, p.cross_size*0.5),
                (-p.bar_width*0.5, p.bar_width*0.5),
                (-p.cross_size*0.5, p.bar_width*0.5),
                (-p.cross_size*0.5, -p.bar_width*0.5),
                (-p.bar_width*0.5, -p.bar_width*0.5),
                (-p.bar_width*0.5, -p.cross_size*0.5),
                (p.bar_width*0.5, -p.cross_size*0.5),
                (p.bar_width*0.5, -p.bar_width*0.5),
                (p.cross_size*0.5, -p.bar_width*0.5)
        ]

        gap = [
                (p.cross_size*0.5+p.cross_gap, p.bar_width*0.5+p.cross_gap),
                (p.bar_width*0.5+p.cross_gap, p.bar_width*0.5+p.cross_gap),
                (p.bar_width*0.5+p.cross_gap, p.cross_size*0.5+p.cross_gap),
                (-p.bar_width*0.5-p.cross_gap, p.cross_size*0.5+p.cross_gap),
                (-p.bar_width*0.5-p.cross_gap, p.bar_width*0.5+p.cross_gap),
                (-p.cross_size*0.5-p.cross_gap, p.bar_width*0.5+p.cross_gap),
                (-p.cross_size*0.5-p.cross_gap, -p.bar_width*0.5-p.cross_gap),
                (-p.bar_width*0.5-p.cross_gap, -p.bar_width*0.5-p.cross_gap),
                (-p.bar_width*0.5-p.cross_gap, -p.cross_size*0.5-p.cross_gap),
                (p.bar_width*0.5+p.cross_gap, -p.cross_size*0.5-p.cross_gap),
                (p.bar_width*0.5+p.cross_gap, -p.bar_width*0.5-p.cross_gap),
                (p.cross_size*0.5+p.cross_gap, -p.bar_width*0.5-p.cross_gap)
        ]

        cross=shapely.Polygon(cross)
        gap=shapely.Polygon(gap)

        polys = [cross, gap]
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [cross, gap] = polys

        #subtracts out ground plane on the layer it's on
        self.add_qgeometry('poly',
                           dict(gap=gap),
                           subtract=True,
                           layer=p.layer)
        
        # Adds the object to the qgeometry table
        # TODO: this block can cause deletion of other qgeometry components?? 
        self.add_qgeometry('poly',
                           dict(cross=cross),
                           # subtract=True,
                           layer=p.layer)



class MarkerSquare(QComponent):
    """Create a solitary square marker.

    Inherits QComponent class.

    Square marker either has a Metal or ground cutout Geometry as specified by is_ground_cutout:
        * square_width  - Width of square along x-axis
        * square_height - Height of square along y-axis
        * is_ground_cutout - If True, the square is a ground cutout, otherwise it is a metallic square.

    As usual, the centre-positioning can be done dynamically via (pos_x,pos_y) with orientation specifying an optional rotation.
        
    Pins:
        There are no pins given that overlap and precise positioning is usually the concern...

    Sketch:
        Below is a sketch of the solitary square marker
        ::

             <--W-->    
             _______
            |       |  /|\    W = square_width
            |   X   |   |H    H = square_height
            |_______|  \|/    X = (pos_x, pos_y)

    .. image::
        Cap3Interdigital.png

    .. meta::
        Solitary square marker

    Default Options:
        * pos_x='0um',pos_y='0um'
        * square_width='20um'
        * square_height='20um'
        * is_ground_cutout=False
    """

    default_options = Dict(pos_x='0um',pos_y='0um', orientation=0,
                           square_width='20um',
                           square_height='20um',
                           is_ground_cutout=False)

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################

        square = [
                (+p.square_width*0.5, +p.square_height*0.5),
                (-p.square_width*0.5, +p.square_height*0.5),
                (-p.square_width*0.5, -p.square_height*0.5),
                (+p.square_width*0.5, -p.square_height*0.5),
        ]
        square = shapely.Polygon(square)

        polys = [square]
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [square] = polys

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(square=square),
                           layer=p.layer,
                           subtract=p.is_ground_cutout)


class MarkerSquare4(QComponent):
    """Create a solitary square marker.

    Inherits QComponent class.

    Square marker either has a Metal or ground cutout Geometry as specified by is_ground_cutout:
        * square_spacing_x - Spacing of the 4 squares along the x-axis (from centre to centre of the squares)
        * square_spacing_y - Spacing of the 4 squares along the y-axis (from centre to centre of the squares)
        * square_width  - Width of each individual square along x-axis
        * square_height - Height of each individual square along y-axis
        * is_ground_cutout - If True, the squares are ground cutouts, otherwise they are metallic squares.

    As usual, the centre-positioning can be done dynamically via (pos_x,pos_y) with orientation specifying an optional rotation.
        
    Pins:
        There are no pins given that overlap and precise positioning is usually the concern...

    Sketch:
        Below is a sketch of the 4-square marker
        ::

             <--W-->    
             _______             _______ 
            |       |  /|\      |       |              W = square_width
            |       |   |H      |       |  /|\         H = square_height
            |_______|  \|/      |_______|   |          X = (pos_x, pos_y)
                                            |          w = square_spacing_x
                          X                 h          h = square_spacing_y
             _______             _______    |
            |       |           |       |   |
            |       |           |       |  \|/
            |_______|           |_______|
                <---------w--------->

    .. image::
        Cap3Interdigital.png

    .. meta::
        4-square marker

    Default Options:
        * pos_x='0um',pos_y='0um'
        * square_spacing_x='320um'
        * square_spacing_y='320um'
        * square_width='20um'
        * square_height='20um'
        * is_ground_cutout=False
    """

    default_options = Dict(pos_x='0um',pos_y='0um', orientation=0,
                           square_spacing_x='320um',
                           square_spacing_y='320um',
                           square_width='20um',
                           square_height='20um',
                           is_ground_cutout=False)

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################

        square = [
                (+p.square_width*0.5, +p.square_height*0.5),
                (-p.square_width*0.5, +p.square_height*0.5),
                (-p.square_width*0.5, -p.square_height*0.5),
                (+p.square_width*0.5, -p.square_height*0.5),
        ]
        square = np.array(square)

        square1 = shapely.Polygon(square + np.array([[-p.square_spacing_x*0.5, p.square_spacing_y*0.5]]*4))
        square2 = shapely.Polygon(square + np.array([[p.square_spacing_x*0.5, p.square_spacing_y*0.5]]*4))
        square3 = shapely.Polygon(square + np.array([[-p.square_spacing_x*0.5, -p.square_spacing_y*0.5]]*4))
        square4 = shapely.Polygon(square + np.array([[p.square_spacing_x*0.5, -p.square_spacing_y*0.5]]*4))

        polys = [square1, square2, square3, square4]
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [square1, square2, square3, square4] = polys

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(square1=square1, square2=square2, square3=square3, square4=square4),
                           layer=p.layer,
                           subtract=p.is_ground_cutout)


class MarkerSquarePocket(QComponent):
    """Create a solitary square marker in a pocket. Deters cheese holes.
        * square_width  - Width of square along x-axis
        * square_height - Height of square along y-axis
        * pocket_distance - Distance edge of marker is from edge of pocket
    """

    default_options = Dict(pos_x='0um',pos_y='0um', orientation=0,
                           square_width='20um',
                           square_height='20um',
                           pocket_distance='5um')

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################


        square = draw.rectangle(p.square_width, p.square_height, 0, 0)
        pocket = draw.rectangle(p.square_width+p.pocket_distance, p.square_width+p.pocket_distance, 0, 0)

        polys = [square, pocket]
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [square, pocket] = polys

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(square=square),
                           layer=p.layer)
        
        self.add_qgeometry('poly',
                           dict(pocket=pocket),
                           layer=p.layer,
                           subtract=True)
        
class Pam(QComponent):
    """Pam markers, the default values should be the same default values in cJob

    Inherits QComponent class.

    Square marker either has a Metal (with a pocket) or ground cutout Geometry as specified by is_ground_cutout:
        * array_x - number of squares along x axis
        * array_y - number of squares along y axis
        * square_spacing_x - Initial spacing of the squares along the x-axis (from centre to centre of the squares)
        * square_spacing_y - Initial spacing of the 4 squares along the y-axis (from centre to centre of the squares)
        * square_width  - Width of each individual square along x-axis
        * square_height - Height of each individual square along y-axis
        * is_ground_cutout - If True, the squares are ground cutouts, otherwise they are metallic squares.
        * pocket_distance - Distance edge of marker is from edge of pocket. if is_ground_cutout=True should be 0
        * inc_x - spacing change after each square in x direction. will making spacing smaller if positive
        * inc_y - spacing change after each square in y direction. will making spacing smaller if positive

    As usual, the centre-positioning can be done dynamically via (pos_x,pos_y) with orientation specifying an optional rotation.
        
    Pins:
        There are no pins given that overlap and precise positioning is usually the concern...

    Sketch:
        Below is a sketch of the 4-square marker
        ::

             <--W-->    
             _______             _______       _______
            |       |  /|\      |       |     |       | W = square_width
            |       |   |H      |       |  /|\|       | H = square_height
            |_______|  \|/      |_______|   | |_______| X = (pos_x, pos_y)
                                            |           w = square_spacing_x
                          X                 h           h = square_spacing_y
             _______             _______    |  _______  l = inc_x
            |       |           |       |   | |       |
            |       |           |       |  \|/|       |
            |_______|           |_______|     |_______|
                <---------w---------><----(w-l)--->


    Default Options:
        * pos_x='0um',pos_y='0um'
        * array_x=5, array_y=5
        * square_spacing_x='75um'
        * square_spacing_y='75um'
        * square_width='20um'
        * square_height='20um'
        * is_ground_cutout=False
        * pocket_distance='5um'
        * inc_x='1um',inc_y='1um'
    """

    default_options = Dict(pos_x='0um',pos_y='0um', orientation=0,
                           array_x=5, array_y=5,
                           square_spacing_x='75um',
                           square_spacing_y='75um',
                           square_width='20um',
                           square_height='20um',
                           is_ground_cutout=False,
                           pocket_distance='5um',
                           inc_x='1um', inc_y='1um')

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################
        if p.is_ground_cutout and p.pocket_distance != 0:
            self.logger.info(
                'Warning: You are making pocket inside a pocket. Are you sure you want to do this?'
            )
        polys=[]#fill with squares of the Pam
        pockets=[]
        polys_dict=dict()
        pockets_dict=dict()
        y=0#start with top left square at origin. translate to center Pam at origin later
        names=[]
        names_pocket=[]
        for row in np.arange(p.array_y):#do first row then move to next row
            x=0#start with left square at x=0. translate to center Pam at origin later
            for column in np.arange(p.array_x):
               name='square'+str(row)+'_'+str(column)#name of the square
               name_pocket='pocket_'+name
               names.append(name)
               names_pocket.append(name_pocket)
               locals().update({name: draw.rectangle(p.square_width, p.square_height, x, y)})
               exec('polys.append('+name+')')
               pocket=draw.rectangle(p.square_width+p.pocket_distance, p.square_height+p.pocket_distance, x, y)
               pockets.append(pocket)
               x=x+p.square_spacing_x-(p.inc_x*column) 
            y=y-p.square_spacing_y+(p.inc_y*row)


        polys = draw.translate(polys, -(x/2), -y/2)
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, p.pos_x, p.pos_y)

        pockets = draw.translate(pockets, -(x/2), -y/2)
        pockets = draw.rotate(pockets, p.orientation, origin=(0, 0))
        pockets = draw.translate(pockets, p.pos_x, p.pos_y) 

        for index in np.arange(len(names)):
            polys_dict.update({names[index]: polys[index]})
            pockets_dict.update({names_pocket[index]: pockets[index]})
        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           polys_dict,
                           layer=p.layer,
                           subtract=p.is_ground_cutout)
        if p.pocket_distance != 0:
            self.add_qgeometry('poly',
                            pockets_dict,
                            layer=p.layer,
                            subtract=True)