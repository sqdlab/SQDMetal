# -*- coding: utf-8 -*-
# Author: Prasanna Pakkiam
# Creation Date: 02/06/2023
# Description: Collection of classes to draw bandages onto pins.

from qiskit_metal import draw
from qiskit_metal.toolbox_python.attr_dict import Dict
from qiskit_metal.qlibrary.core import QComponent
import numpy as np
import shapely

class BandageRectPin(QComponent):
    """Create a rectangular bandage on a pin

    Inherits QComponent class.

    Square marker either has a Metal or ground cutout Geometry as specified by is_ground_cutout:
        * square_width  - Width of square along x-axis
        * square_height - Height of square along y-axis

    The positioning can be done dynamically via:
        * pin_inputs=Dict(start_pin=Dict(component=f'...',pin='...')) - Specifying centre position via a component pin
    The resulting bandage is right on this pin. This class ignores pos_x, pos_y and orientation...
        
    Pins:
        There are no pins given that overlap and precise positioning is usually the concern...

    Sketch:
        Below is a sketch of the solitary square marker
        ::

             <--W-->    
             _______
            |       |  /|\    W = width
            |   X   |   |H    H = height
            |_______|  \|/    X = location of target pin

    .. image::
        Cap3Interdigital.png

    .. meta::
        Rectangular bandage on a pin

    Default Options:
        * width='10um'
        * height='10um'
    """

    default_options = Dict(width='10um',
                           height='10um'
                           )

    def __init__(self, design,
                    name: str = None,
                    options: Dict = None,
                    type: str = "CPW",
                    **kwargs):
        #QRoute forces an end-pin to exist... So make it artificial...
        assert 'pin_inputs' in options, "Must provide a starting pin input via \'pin_inputs\'."
        assert 'start_pin' in options.pin_inputs, "Must provide \'start_pin\' in \'pin_inputs\'."
        assert 'component' in options.pin_inputs.start_pin, "Must provide \'component\' in \'start_pin\'."
        assert 'pin' in options.pin_inputs.start_pin, "Must provide \'pin\' in \'start_pin\'."
        super().__init__(design, name, options, **kwargs)
        #TODO: Perhaps a pull request to add poppable options?

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################

        start_point = self.design.components[self.options.pin_inputs.start_pin.component].pins[self.options.pin_inputs.start_pin.pin]
        startPt = start_point['middle']
        normal = -start_point['normal']
        rot_angle = np.arctan2(normal[1], normal[0])

        bandaid = [
                  (-p.width*0.5, -p.height*0.5),
                  (-p.width*0.5, p.height*0.5),
                  (p.width*0.5, p.height*0.5),
                  (p.width*0.5, -p.height*0.5)]

        bandaid = shapely.Polygon(bandaid)

        polys = [bandaid]
        polys = draw.rotate(polys, rot_angle, origin=(0, 0), use_radians=True)
        polys = draw.translate(polys, *startPt)
        [bandaid] = polys

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(bandaid=bandaid),
                           layer=p.layer)


