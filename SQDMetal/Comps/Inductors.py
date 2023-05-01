# -*- coding: utf-8 -*-
# Author: Prasanna Pakkiam
# Creation Date: 28/04/2023
# Description: Collection of classes to dynamically route inductors.

from qiskit_metal import draw
from qiskit_metal.toolbox_python.attr_dict import Dict
from qiskit_metal.qlibrary.core import QComponent
import numpy as np
import shapely
from SQDMetal.Utilities.CpwParams import CpwParams
from SQDMetal.Utilities.QUtilities import QUtilities

class InductorMeander(QComponent):
    """Draws a meandering inductor with sharp corners as defined by the meander width, spacing and number
    of long sections (see 10.2298/sjee0403057s).

    Inherits QComponent class.

    Inductor Metal Geometry and Ground Cutout Pocket:
        * track_width - Width of the track forming the inductor
        * meander_spacing - Spacing of the adjacent meanders (measured from the track centres)
        * meander_width - Width of the meandering inductor (measured from the track centres)
        * num_long_sections - Number of long sections (i.e. on each meander)
    
    The spacing (i.e. cuts into the ground plane) can be controlled via:
        * mean_gap - The ground plane is already cut out from the outer start to outer end of the meandering section. Its width away
                     from the inductor cut in extra is given by side_gap.
        * init_gap - This adds spacing to the inlet and outlet sections to the ground plane. If it is set to zero, it is assumed to be
                     in line with the cut-out used in the meandering section. Note that if bounding box of the inlets is less than the
                     inductor width, there is an init_gap padding on the start and end sections of the meander so that the inductor does
                     not touch the ground plane.
    
    As usual, the positioning can be done dynamically as a vector given by the supplied parameters: (pos_x,pos_y) to (end_x,end_y). The
    resulting inductor is right in the centre...
        
    Pins:
        There are two pins on the capacitor at either end
        The pins attach directly to the built in lead length and only needs a width defined
        * cpw_width - center trace width of the trace lead line and cap fingers

    Sketch:
        Below is a sketch of the capacitor
        ::
                _____       _____
               |     |     |  /\ |
        _______|<-d->|     |   : |      _______
                     |     |  h: |     |
                     |_____|  \/ |_____|

        d = meander_spacing
        h = meander_width
        In the diagram above, num_long_sections is 3 (i.e. 3 long sections)

    .. image::
        Cap3Interdigital.png

    .. meta::
        Cap 3 Interdigital

    Default Options:
        * pos_x='0um',pos_y='0um'
        * end_x='400um',end_y='0um'
        * track_width='2um',
        * meander_spacing='4um'
        * meander_width='30um'
        * num_long_sections=67
        * mean_gap='2um'
        * init_gap='0um'
    """

    #  Define structure functions

    default_options = Dict(pos_x='0um',pos_y='0um',
                           end_x='400um',end_y='0um',
                           track_width='2um',
                           meander_spacing='4um',
                           meander_width='30um',
                           num_long_sections=67,
                           mean_gap='2um',
                           init_gap='0um')
    """Default drawing options"""

    TOOLTIP = """Create a three finger planar capacitor with a ground pocket cuttout."""

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################

        # Make the shapely polygons for the main cap structure
        poly, padGap, pin1, pin2 = InductorMeander._draw_inductor(p)

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(inductor=poly),
                           layer=p.layer)

        # subtracts out ground plane on the layer its on
        self.add_qgeometry('poly',
                           dict(padGap=padGap),
                           subtract=True,
                           layer=p.layer)

        # Generates its own pins
        self.add_pin('a', pin1.coords[::], width=p.track_width)
        self.add_pin('b', pin2.coords[::], width=p.track_width)

    @staticmethod
    def _draw_inductor(p):
        len_trace = np.sqrt((p.end_x-p.pos_x)**2+(p.end_y-p.pos_y)**2)

        meander = [(0,0), (0, p.meander_width*0.5), (p.meander_spacing, p.meander_width*0.5)]
        for m in range(p.num_long_sections):
            meander += [(meander[-1][0], -meander[-1][1])]
            meander += [(meander[-1][0]+p.meander_spacing, meander[-1][1])]
        meander += [(meander[-1][0], 0)]

        len_meander = meander[-1][0]
        meander = [((len_meander-len_trace)*0.5, 0)] + meander + [((len_meander+len_trace)*0.5, 0)]
        meander = np.array(meander)
        meander[:,0] -= meander[0,0]

        poly = shapely.LineString(meander).buffer(p.track_width*0.5, join_style=2, cap_style=2)

        pin1 = [(0,-0.5*p.track_width), (0,0.5*p.track_width)]
        pin2 = [(len_trace,0.5*p.track_width), (len_trace,-0.5*p.track_width)]
        pin1 = shapely.LineString(pin1)
        pin2 = shapely.LineString(pin2)

        mean_start_x = (len_trace-len_meander)*0.5-p.track_width*0.5
        mean_stop_x = (len_trace+len_meander)*0.5+p.track_width*0.5
        mean_stop_y = (p.meander_width+p.track_width)*0.5 + p.mean_gap
        init_stop_y = p.track_width*0.5+p.init_gap
        if p.init_gap > 0 and init_stop_y < mean_stop_y:
            mean_start_x -= p.init_gap
            mean_stop_x += p.init_gap
        padGap = shapely.Polygon([(mean_start_x,mean_stop_y), (mean_start_x,-mean_stop_y), (mean_stop_x,-mean_stop_y), (mean_stop_x,mean_stop_y)])
        if p.init_gap == 0:
            init_stop_y = mean_stop_y
        padGap1 = shapely.Polygon([(0,init_stop_y), (0,-init_stop_y), (mean_start_x,-init_stop_y), (mean_start_x,init_stop_y)])
        padGap2 = shapely.Polygon([(mean_stop_x,init_stop_y), (mean_stop_x,-init_stop_y), (len_trace,-init_stop_y), (len_trace,init_stop_y)])
        padGap = shapely.unary_union([padGap, padGap1, padGap2])

        polys = [poly, padGap, pin1, pin2]
        polys = draw.rotate(polys, np.arctan2(p.end_y-p.pos_y,p.end_x-p.pos_x), origin=(0, 0), use_radians=True)
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [poly, padGap, pin1, pin2] = polys

        return poly, padGap, pin1, pin2

class InductorMeanderPinStretch(QComponent):
    """Draws a meandering inductor with sharp corners as defined by the meander width, spacing and number
    of long sections (see 10.2298/sjee0403057s).

    Inherits QComponent class.

    Inductor Metal Geometry and Ground Cutout Pocket:
        * track_width - Width of the track forming the inductor
        * meander_spacing - Spacing of the adjacent meanders (measured from the track centres)
        * meander_width - Width of the meandering inductor (measured from the track centres)
        * num_long_sections - Number of long sections (i.e. on each meander)
    
    The spacing (i.e. cuts into the ground plane) can be controlled via:
        * mean_gap - The ground plane is already cut out from the outer start to outer end of the meandering section. Its width away
                     from the inductor cut in extra is given by side_gap.
        * init_gap - This adds spacing to the inlet and outlet sections to the ground plane. If it is set to zero, it is assumed to be
                     in line with the cut-out used in the meandering section. Note that if bounding box of the inlets is less than the
                     inductor width, there is an init_gap padding on the start and end sections of the meander so that the inductor does
                     not touch the ground plane.
    
    The positioning can be done dynamically via:
        * pin_inputs=Dict(start_pin=Dict(component=f'...',pin='...')) - Specifying start position via a component pin
        * dist_extend - Distance upon to stretch away from the start pin.
    The resulting inductor is right in the centre. This class ignores pos_x, pos_y and orientation...
        
    Pins:
        There are two pins on the capacitor at either end
        The pins attach directly to the built in lead length and only needs a width defined
        * cpw_width - center trace width of the trace lead line and cap fingers

    Sketch:
        Below is a sketch of the capacitor
        ::
                _____       _____
               |     |     |  /\ |
        _______|<-d->|     |   : |      _______
                     |     |  h: |     |
                     |_____|  \/ |_____|

        d = meander_spacing
        h = meander_width
        In the diagram above, num_long_sections is 3 (i.e. 3 long sections)

    .. image::
        Cap3Interdigital.png

    .. meta::
        Cap 3 Interdigital

    Default Options:
        * dist_extend='400um'
        * track_width='2um',
        * meander_spacing='4um'
        * meander_width='30um'
        * num_long_sections=67
        * mean_gap='2um'
        * init_gap='0um'
    """

    #  Define structure functions

    default_options = Dict(dist_extend='400um',
                           track_width='2um',
                           meander_spacing='4um',
                           meander_width='30um',
                           num_long_sections=67,
                           mean_gap='2um',
                           init_gap='0um')
    """Default drawing options"""

    TOOLTIP = """Create a three finger planar capacitor with a ground pocket cuttout."""

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
        norm = start_point['normal']
        p.pos_x = startPt[0]
        p.pos_y = startPt[1]
        endPt = startPt + norm*p.dist_extend
        p.end_x = endPt[0]
        p.end_y = endPt[1]

        self.options.pos_x, self.options.pos_y = startPt

        # Make the shapely polygons for the main cap structure
        poly, padGap, pin1, pin2 = InductorMeander._draw_inductor(p)

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(inductor=poly),
                           layer=p.layer)

        # subtracts out ground plane on the layer its on
        self.add_qgeometry('poly',
                           dict(padGap=padGap),
                           subtract=True,
                           layer=p.layer)

        # Generates its own pins
        self.add_pin('a', pin1.coords[::], width=p.track_width)
        self.add_pin('b', pin2.coords[::], width=p.track_width)
