 # -*- coding: utf-8 -*-
# Author: Prasanna Pakkiam and Divita Gautam
# Creation Date: 02/06/2023
# Modified Date: 17/05/2024 (tapered bandages)
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
        * width  - Width of square along x-axis
        * height - Height of square along y-axis
        * offset_distance  - distance to attach the bandage away from the target component
        * offset_direction - direction in degrees to move the bandage away from the target component
        * fillet_radius - The radius of the curved edges of the bandage 
    The positioning can be done dynamically via:
        * target_comp - Name of the Qiskit metal component which has the pin to be bandaged
        * target_pin  - Name of the pin of the target Qiskit metal component where this bandage sits
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

        fillet_radius and fillet_resolution round off the corners...

    .. image::
        Cap3Interdigital.png

    .. meta::
        Rectangular bandage on a pin

    Default Options:
        * width='10um'
        * height='10um'
        * fillet_radius='0um',
        * fillet_resolution=4,
        * offset_distance='0um',
        * offset_direction=0
    """

    default_options = Dict(width='10um',
                           height='10um',
                           fillet_radius='0um',
                           fillet_resolution=4,
                           offset_distance='0um',
                           offset_direction=0
                           )

    def __init__(self, design,
                    name: str = None,
                    options: Dict = None,
                    type: str = "CPW",
                    **kwargs):
        #QRoute forces an end-pin to exist... So make it artificial...
        assert 'target_comp' in options, "Must provide a target component via \'target_comp\'."
        assert 'target_pin' in options, "Must provide a target component pin via \'target_pin\'."
        super().__init__(design, name, options, **kwargs)
        #TODO: Perhaps a pull request to add poppable options?

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################

        start_point = self.design.components[self.options.target_comp].pins[self.options.target_pin]
        startPt = start_point['middle']+(p.offset_distance*np.cos(p.offset_direction/180*np.pi),p.offset_distance*np.sin(p.offset_direction/180*np.pi))
        normal = -start_point['normal']
        rot_angle = np.arctan2(normal[1], normal[0])

        bandaid = [
                  (-p.width*0.5+p.fillet_radius, -p.height*0.5+p.fillet_radius),
                  (-p.width*0.5+p.fillet_radius, p.height*0.5-p.fillet_radius),
                  (p.width*0.5-p.fillet_radius, p.height*0.5-p.fillet_radius),
                  (p.width*0.5-p.fillet_radius, -p.height*0.5+p.fillet_radius)]

        bandaid = shapely.Polygon(bandaid)
        bandaid = bandaid.buffer(p.fillet_radius, cap_style=1, join_style=1, mitre_limit=2.0, quad_segs=p.fillet_resolution)
        
        polys = [bandaid]
        polys = draw.rotate(polys, rot_angle, origin=(0, 0), use_radians=True)
        polys = draw.translate(polys, *startPt)
        [bandaid] = polys

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(bandaid=bandaid),
                           layer=p.layer)

class BandageTaperedPin(QComponent):
    """Creates a tapered bandage on a pin

        Inherits QComponent class.

        Square marker either has a Metal or ground cutout Geometry as specified by is_ground_cutout:
            * taper_width_top  - top width of tapered bandage
            * taper_width_base - bottom width of tapered bandage
            * taper_height     - hieght of tapered bandage
            * fillet_radius    - fillet radius to curve edges
            * fillet_resolution- number of points used to calculate the fillet/curve
            * offset_distance      - distance to attach the bandage away from the target component
            * offset_direction - direction in degrees to move the bandage away from the target component
            * rotate_angle     - angle of rotation 
            * base_width       - width of the rectangular base
            * base_height      - height of the rectangular base
            * ignore_slope_check - (default False) If True, it performs the Martinis/Divita check of 0.4 for the slope...


        The positioning can be done dynamically via:
            * target_comp - Name of the Qiskit metal component which has the pin to be bandaged
            * target_pin  - Name of the pin of the target Qiskit metal component where this bandage sits
        The resulting bandage is right on this pin. This class ignores pos_x, pos_y and orientation...
            
        Pins:
            There are no pins given that overlap and precise positioning is usually the concern...

        Sketch:
            Below is a sketch of the solitary square marker
            ::

            |x
        ____|   Note the axis...
        y   
        
                   <--t-->    
                    _____
                   |     |        /|\    t = taper_width_top
                 |         |       |h    h = taper_height
             ___|_____x_____|___  \|/    x = location of target pin
         /|\|    <----b---->    |        b = taper_width_base
         H| |                   |        B = base_width
         \|/|___________________|        H = base_height
             <--------B-------->        
            fillet_radius and fillet_resolution round off the corners... Orientation rotates the patch about its pin coordinate...

        .. image::
            Cap3Interdigital.png

        .. meta::
            Rectangular bandage on a pin

        Default Options:
            * taper_width_top='10um',
            * taper_width_base='20um',
            * taper_height='30um',
            * base_width='30um',
            * base_height='20um',
            * fillet_radius='1um',
            * fillet_resolution=4,
            * offset_distance='0um',
            * offset_direction=0
            * ignore_slope_check=False
        """

    default_options = Dict(taper_width_top='10um',
                           taper_width_base='20um',
                           taper_height='30um',
                           base_width='30um',
                           base_height='20um',
                           fillet_radius='1um',
                           fillet_resolution=4,
                           offset_distance='0um',
                           offset_direction=0,
                           ignore_slope_check=False
                           )

    def __init__(self, design,
                    name: str = None,
                    options: Dict = None,
                    type: str = "CPW",
                    **kwargs):
        assert 'target_comp' in options, "Must provide a target component via \'target_comp\'."
        assert 'target_pin' in options, "Must provide a target component pin via \'target_pin\'."


        if not options.get('ignore_slope_check'):
            # Compute and validate slope
            taper_width_top = float(options.get('taper_width_top', '26.8um')[:-2])
            taper_width_base = float(options.get('taper_width_base', '10um')[:-2])
            taper_height = float(options.get('taper_height', '30um')[:-2])

            computed_slope = (2 * taper_height) / (taper_width_base - taper_width_top)
            assert computed_slope < 0.4, f"Slope {computed_slope:.2f} is too steep! Adjust taper dimensions for the slope to be less than 0.4."
            
        super().__init__(design, name, options, **kwargs)

    def make(self):
        p = self.p
        start_point = self.design.components[self.options.target_comp].pins[self.options.target_pin]
        startPt = start_point['middle']+(p.offset_distance*np.cos(p.offset_direction/180*np.pi),p.offset_distance*np.sin(p.offset_direction/180*np.pi))

        bandaid = [(p.taper_width_base/2, 0),
                   (p.taper_width_top/2, p.taper_height),
                   (-p.taper_width_top/2, p.taper_height),
                   (-p.taper_width_base/2, 0),
                   (-p.base_width/2, 0),
                   (-p.base_width/2, -p.base_height),
                   (p.base_width/2, -p.base_height),
                   (p.base_width/2, 0)]
        bandaid = shapely.Polygon(bandaid)
        
        bandaid = bandaid.buffer(p.fillet_radius, join_style=2, cap_style=3).buffer(-p.fillet_radius, cap_style=1, join_style=1, mitre_limit=2.0, quad_segs=p.fillet_resolution)
        bandaid = bandaid.buffer(-p.fillet_radius, join_style=2, cap_style=3).buffer(p.fillet_radius, cap_style=1, join_style=1, mitre_limit=2.0, quad_segs=p.fillet_resolution)
        
        polys = [bandaid]
        polys = draw.rotate(polys,p.orientation-90, origin=(0, 0), use_radians=False)
        polys = draw.translate(polys, *startPt)
        [bandaid] = polys

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(bandaid=bandaid),
                           layer=p.layer)


