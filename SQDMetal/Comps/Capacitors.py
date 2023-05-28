# -*- coding: utf-8 -*-
# Author: Prasanna Pakkiam
# Creation Date: 24/04/2023
# Description: Collection of classes to dynamically route capacitors.

from qiskit_metal import draw
from qiskit_metal.toolbox_python.attr_dict import Dict
from qiskit_metal.qlibrary.core import QComponent
import numpy as np
import shapely
from SQDMetal.Utilities.CpwParams import CpwParams
from SQDMetal.Utilities.QUtilities import QUtilities

class CapacitorInterdigital(QComponent):
    """Create an interdigital capacitor on a CPW.
    The width of the fingers is determined by cpw_width.

    Inherits QComponent class.

    Capacitor Metal Geometry and Ground Cutout Pocket:
        * len_diag - Length of the staggered region that bridges from CPW to the capacitor
        * len_flat - Length of the flat region before starting onto the fingers
        * fing_len - Actual finger length
        * fing_len_gap - Gap between finger and the opposite capacitor conductor
        * fing_wid - Finger width
        * fing_wid_gap - Gap between adjacent fingers
        * N_total - Total number of fingers
        * larger_first - If True, for odd N_total, the larger number of fingers (i.e. (N+1)/2) will be on the first pad's conductor
    
    The spacing (i.e. cuts into the ground plane) can be controlled via:
        * side_gap - If this is zero, then the gap on the sides of the capacitor is calculated via a 50ohm impedance CPW line. Otherwise,
                     it is set via the given side gap.
        * init_pad - This adds spacing to the ground plane on the feed lines. This is useful when the diagonal section is steep (e.g.
                     when len_diag is zero) to ensure that the ground plane does not intersect with the main capacitor conductors. The
                     ground plane spacing typically starts to change on meeting with the LD or LF sections (see below). If init_pad > 0,
                     it starts from an earlier point.
        * side_gap_uniform - If True, the spacing on the sides are made uniform by taking the maximum extent of either the gap around
                             the inlet feed or the main capacitor body with the fingers.
    
    As usual, the positioning can be done dynamically as a vector given by the supplied parameters: (pos_x,pos_y) to (end_x,end_y)
        
    Pins:
        There are two pins on the capacitor at either end
        The pins attach directly to the built in lead length and only needs a width defined
        * cpw_width - center trace width of the trace lead line and cap fingers

    Sketch:
        Below is a sketch of the capacitor
        ::

        @@@@@   |   @@@@@                               When setting init_pad > 0:
        @@@@@   |   @@@@@                                   @@@@@   |   @@@@@
        @@@@@   |   @@@@@                                   @@@@@   |   @@@@@
        @@@@@   |   @@@@@   @   = Ground Plane              @@      |   IP @@
        @@@@@   -   LD  @   -,| = Conductors                @      ---  LD  @
        @@@@   ---  LD  @                                   @    -------LF  @
        @@@   ----- LD  @   LD = len_diag                   @    -------LF  @ 
        @    -------LF  @                                   @    | | | |    @
        @    -------LF  @   LF = len_flat                   @    |||||||    @   IP = init_pad
        @ FL | | | |FLG @   FLG = fing_len_gap              @    |||||||    @
        @ FL |||||||    @                                   @     | | |     @ 
        @ FL |||||||SSSS@   S = side_gap                    @    -------LF  @
        @     | | | FLG @   FL = fing_len                   @    -------LF  @ 
        @    -------LF  @                                   @      ---  LD  @
        @    -------LF  @   FP = Front Pad                  @@      |   IP @@
        @@@   ----- LD  @                                   @@@@@   |   @@@@@
        @@@@   ---  LD  @                                   @@@@@   |   @@@@@
        @@@@@   -   LD  @
        @@@@@   |       @   
        @@@@@   |       @
        @@@@@   |       @
        @@@@@   |   @@@@@
        @@@@@   |   @@@@@

    .. image::
        Cap3Interdigital.png

    .. meta::
        Cap 3 Interdigital

    Default Options:
        * trace_width: '10um
        * pos_x='0um',pos_y='0um'
        * end_x='50um',end_y='0um'
        * cpw_width='10um'
        * len_diag='5um'
        * len_flat='5um'
        * fing_len='10um'
        * fing_len_gap='1um'
        * fing_wid='2um'
        * fing_wid_gap='1um'
        * N_total=5
        * larger_first=True
        * side_gap_uniform=False
        * side_gap='0um'
        * init_pad='0um'
    """

    #  Define structure functions

    default_options = Dict(pos_x='0um',pos_y='0um',
                           end_x='50um',end_y='0um',
                           cpw_width='10um',
                           len_diag='5um',
                           len_flat='5um',
                           fing_len='10um',
                           fing_len_gap='1um',
                           fing_wid='2um',
                           fing_wid_gap='1um',
                           N_total=5,
                           larger_first=True,
                           side_gap_uniform=False,
                           side_gap='0um',
                           init_pad='0um')
    """Default drawing options"""

    TOOLTIP = """Create a three finger planar capacitor with a ground pocket cuttout."""

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################

        pad1, pad2, padGap, pin1, pin2 = CapacitorInterdigital._draw_capacitor(p, self._design)

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(pad1=pad1, pad2=pad2),
                           layer=p.layer)

        #subtracts out ground plane on the layer it's on
        self.add_qgeometry('poly',
                           dict(padGap=padGap),
                           subtract=True,
                           layer=p.layer)

        # Generates its own pins
        self.add_pin('a', pin1.coords[::-1], width=p.cpw_width)
        self.add_pin('b', pin2.coords[::-1], width=p.cpw_width)
    
    @staticmethod
    def _draw_capacitor(p, design, discard_leads=False):
        Num = int(p.N_total)

        # Make the shapely polygons for the main cap structure
        len_fings_plus_gap = p.fing_len + p.fing_len_gap
        cap_width = Num*p.fing_wid + (Num-1)*p.fing_wid_gap
        len_trace = np.sqrt((p.end_x-p.pos_x)**2+(p.end_y-p.pos_y)**2)
        #Tracing the left wall and the bottom-side of the pad
        startX = 0
        if discard_leads:
            startX = len_trace*0.5-len_fings_plus_gap*0.5-p.len_flat-p.len_diag
        pad1a = [
                (startX, p.cpw_width*0.5),
                (startX, -p.cpw_width*0.5),
                (len_trace*0.5-len_fings_plus_gap*0.5-p.len_flat-p.len_diag, -p.cpw_width*0.5),
                (len_trace*0.5-len_fings_plus_gap*0.5-p.len_flat, -cap_width*0.5),
                (len_trace*0.5-len_fings_plus_gap*0.5, -cap_width*0.5)]
        #Tracing the top-side of the pad
        pad1b = [
                (len_trace*0.5-len_fings_plus_gap*0.5, cap_width*0.5),
                (len_trace*0.5-len_fings_plus_gap*0.5-p.len_flat, cap_width*0.5),
                (len_trace*0.5-len_fings_plus_gap*0.5-p.len_flat-p.len_diag, p.cpw_width*0.5),
                (startX, p.cpw_width*0.5)]
        pad2a = np.array(pad1a)
        pad2b = np.array(pad1b)
        #Do the fingers
        if Num % 2 == 0:
            num_fings_1 = Num // 2
            num_fings_2 = num_fings_1
        else:
            if p.larger_first:
                num_fings_1 = (Num+1) // 2
                num_fings_2 = (Num-1) // 2
            else:
                num_fings_1 = (Num-1) // 2
                num_fings_2 = (Num+1) // 2
        fing_coords1 = []
        cur_x = pad1a[-1][0]
        cur_y = pad1a[-1][1]
        if num_fings_1 >= num_fings_2:
            pad1a.pop(-1)
        else:
            cur_y += p.fing_wid + p.fing_wid_gap
        for m in range(num_fings_1):
            #Basically drawing:
            #  e
            #  d.......c
            #  a.......b
            fing_coords1 += [(cur_x, cur_y), (cur_x+p.fing_len, cur_y), (cur_x+p.fing_len, cur_y+p.fing_wid), (cur_x, cur_y+p.fing_wid), (cur_x, cur_y+p.fing_wid*2+p.fing_wid_gap*2)]
            cur_y = fing_coords1[-1][1]
        pad1 = pad1a + fing_coords1[:-1] + pad1b

        if Num % 2 == 0:
            pad2 = np.array(pad1)
            pad2[:,0] = len_trace - pad2[:,0]
            pad2[:,1] = -pad2[:,1]
            pad2 = pad2[::-1]
        else:
            pad2a[:,0] = len_trace - pad2a[:,0]
            pad2b[:,0] = len_trace - pad2b[:,0]
            pad2a = pad2a.tolist()
            
            fing_coords2 = []
            cur_x = pad2a[-1][0]
            cur_y = pad2a[-1][1]
            if num_fings_2 > num_fings_1:
                pad2a.pop(-1)
            else:
                cur_y += p.fing_wid + p.fing_wid_gap
            for m in range(num_fings_2):
                #Basically drawing:
                #          e
                #  c.......d
                #  b.......a
                fing_coords2 += [(cur_x, cur_y), (cur_x-p.fing_len, cur_y), (cur_x-p.fing_len, cur_y+p.fing_wid), (cur_x, cur_y+p.fing_wid), (cur_x, cur_y+p.fing_wid*2+p.fing_wid_gap*2)]
                cur_y = fing_coords2[-1][1]
            pad2 = np.concatenate([pad2a, fing_coords2[:-1], pad2b])
            pad2 = pad2[::-1]

        units = QUtilities.get_units(design)
        cpwP = CpwParams.fromQDesign(design)
        gap_cpw_line = cpwP.get_gap_from_width(p.cpw_width*units)/units
        if p.side_gap == 0:
            gap_cpw_cap = cpwP.get_gap_from_width(cap_width*units)/units
        else:
            gap_cpw_cap = p.side_gap
        if p.side_gap_uniform:
            max_extent = max(p.cpw_width*0.5+gap_cpw_line, cap_width*0.5+gap_cpw_cap)
            gap_cpw_line = max_extent - p.cpw_width*0.5
            gap_cpw_cap = max_extent - cap_width*0.5
        #
        padGap = np.array([
                 (len_trace*0.5-len_fings_plus_gap*0.5-p.len_flat, cap_width*0.5 + gap_cpw_cap),
                 (len_trace*0.5-len_fings_plus_gap*0.5-p.len_flat-p.len_diag, p.cpw_width*0.5 + gap_cpw_line),
                 (startX, p.cpw_width*0.5 + gap_cpw_line),

                 (startX, -p.cpw_width*0.5 - gap_cpw_line),
                 (len_trace*0.5-len_fings_plus_gap*0.5-p.len_flat-p.len_diag, -p.cpw_width*0.5 - gap_cpw_line),
                 (len_trace*0.5-len_fings_plus_gap*0.5-p.len_flat, -cap_width*0.5 - gap_cpw_cap),

                 (len_trace*0.5+len_fings_plus_gap*0.5+p.len_flat, -cap_width*0.5 - gap_cpw_cap),
                 (len_trace*0.5+len_fings_plus_gap*0.5+p.len_flat+p.len_diag, -p.cpw_width*0.5 - gap_cpw_line),
                 (len_trace-startX, -p.cpw_width*0.5 - gap_cpw_line),

                 (len_trace-startX, p.cpw_width*0.5 + gap_cpw_line),
                 (len_trace*0.5+len_fings_plus_gap*0.5+p.len_flat+p.len_diag, p.cpw_width*0.5 + gap_cpw_line),
                 (len_trace*0.5+len_fings_plus_gap*0.5+p.len_flat, cap_width*0.5 + gap_cpw_cap),
                 ])
        if p.len_diag == 0:
            padGap[[0,1,4,5],0] -= p.init_pad
            padGap[[6,7,10,11],0] += p.init_pad
        else:
            padGap[[1,4],0] -= p.init_pad
            padGap[[7,10],0] += p.init_pad
        if discard_leads:
            padGap[[2,3],0] -= p.init_pad
            padGap[[8,9],0] += p.init_pad

        pin1 = pad1[:2]
        pin2 = pad2[[-2,-1]]
        #
        pad1 = shapely.Polygon(pad1)
        pad2 = shapely.Polygon(pad2)
        padGap = shapely.Polygon(padGap)
        pin1 = shapely.LineString(pin1)
        pin2 = shapely.LineString(pin2)

        polys = [pad1, pad2, padGap, pin1, pin2]
        polys = draw.rotate(polys, np.arctan2(p.end_y-p.pos_y,p.end_x-p.pos_x), origin=(0, 0), use_radians=True)
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [pad1, pad2, padGap, pin1, pin2] = polys

        return pad1, pad2, padGap, pin1, pin2

class CapacitorInterdigitalPinStretch(QComponent):
    """Create an interdigital capacitor on a CPW.
    The width of the fingers is determined by cpw_width.

    Inherits QComponent class.

    Capacitor Metal Geometry and Ground Cutout Pocket:
        * len_diag - Length of the staggered region that bridges from CPW to the capacitor
        * len_flat - Length of the flat region before starting onto the fingers
        * fing_len - Actual finger length
        * fing_len_gap - Gap between finger and the opposite capacitor conductor
        * fing_wid - Finger width
        * fing_wid_gap - Gap between adjacent fingers
        * N_total - Total number of fingers
        * larger_first - If True, for odd N_total, the larger number of fingers (i.e. (N+1)/2) will be on the first pad's conductor
    
    The spacing (i.e. cuts into the ground plane) can be controlled via:
        * side_gap - If this is zero, then the gap on the sides of the capacitor is calculated via a 50ohm impedance CPW line. Otherwise,
                     it is set via the given side gap.
        * init_pad - This adds spacing to the ground plane on the feed lines. This is useful when the diagonal section is steep (e.g.
                     when len_diag is zero) to ensure that the ground plane does not intersect with the main capacitor conductors. The
                     ground plane spacing typically starts to change on meeting with the LD or LF sections (see below). If init_pad > 0,
                     it starts from an earlier point.
        * side_gap_uniform - If True, the spacing on the sides are made uniform by taking the maximum extent of either the gap around
                             the inlet feed or the main capacitor body with the fingers.
    
    The positioning can be done dynamically via:
        * pin_inputs=Dict(start_pin=Dict(component=f'...',pin='...')) - Specifying start position via a component pin
        * dist_extend - Distance upon to stretch away from the start pin.
    The resulting capacitor is right in the centre. This class ignores pos_x, pos_y and orientation...
        
    Pins:
        There are two pins on the capacitor at either end
        The pins attach directly to the built in lead length and only needs a width defined
        * cpw_width - center trace width of the trace lead line and cap fingers

    Sketch:
        Below is a sketch of the capacitor
        ::

        @@@@@   |   @@@@@                               When setting init_pad > 0:
        @@@@@   |   @@@@@                                   @@@@@   |   @@@@@
        @@@@@   |   @@@@@                                   @@@@@   |   @@@@@
        @@@@@   |   @@@@@   @   = Ground Plane              @@      |   IP @@
        @@@@@   -   LD  @   -,| = Conductors                @      ---  LD  @
        @@@@   ---  LD  @                                   @    -------LF  @
        @@@   ----- LD  @   LD = len_diag                   @    -------LF  @ 
        @    -------LF  @                                   @    | | | |    @
        @    -------LF  @   LF = len_flat                   @    |||||||    @   IP = init_pad
        @ FL | | | |FLG @   FLG = fing_len_gap              @    |||||||    @
        @ FL |||||||    @                                   @     | | |     @ 
        @ FL |||||||SSSS@   S = side_gap                    @    -------LF  @
        @     | | | FLG @   FL = fing_len                   @    -------LF  @ 
        @    -------LF  @                                   @      ---  LD  @
        @    -------LF  @   FP = Front Pad                  @@      |   IP @@
        @@@   ----- LD  @                                   @@@@@   |   @@@@@
        @@@@   ---  LD  @                                   @@@@@   |   @@@@@
        @@@@@   -   LD  @
        @@@@@   |       @   
        @@@@@   |       @
        @@@@@   |       @
        @@@@@   |   @@@@@
        @@@@@   |   @@@@@

    .. image::
        Cap3Interdigital.png

    .. meta::
        Cap 3 Interdigital

    Default Options:
        * trace_width: '10um
        * dist_extend='50um'
        * cpw_width='10um'
        * len_diag='5um'
        * len_flat='5um'
        * fing_len='10um'
        * fing_len_gap='1um'
        * fing_wid='2um'
        * fing_wid_gap='1um'
        * N_total=5
        * larger_first=True
        * side_gap_uniform=False
        * side_gap='0um'
        * init_pad='0um'
    """

    #  Define structure functions

    default_options = Dict(dist_extend='50um',
                           cpw_width='10um',
                           len_diag='5um',
                           len_flat='5um',
                           fing_len='10um',
                           fing_len_gap='1um',
                           fing_wid='2um',
                           fing_wid_gap='1um',
                           N_total=5,
                           larger_first=True,
                           side_gap_uniform=False,
                           side_gap='0um',
                           init_pad='0um')
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

        pad1, pad2, padGap, pin1, pin2 = CapacitorInterdigital._draw_capacitor(p, self._design)

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(pad1=pad1, pad2=pad2),
                           layer=p.layer)

        #subtracts out ground plane on the layer it's on
        self.add_qgeometry('poly',
                           dict(padGap=padGap),
                           subtract=True,
                           layer=p.layer)

        # Generates its own pins
        self.add_pin('a', pin1.coords[::-1], width=p.cpw_width)
        self.add_pin('b', pin2.coords[::-1], width=p.cpw_width)

class CapacitorInterdigitalPinPin(QComponent):
    """Create an interdigital capacitor on a CPW.
    The width of the fingers is determined by cpw_width.

    Inherits QComponent class.

    Capacitor Metal Geometry and Ground Cutout Pocket:
        * len_diag - Length of the staggered region that bridges from CPW to the capacitor
        * len_flat - Length of the flat region before starting onto the fingers
        * fing_len - Actual finger length
        * fing_len_gap - Gap between finger and the opposite capacitor conductor
        * fing_wid - Finger width
        * fing_wid_gap - Gap between adjacent fingers
        * N_total - Total number of fingers
        * larger_first - If True, for odd N_total, the larger number of fingers (i.e. (N+1)/2) will be on the first pad's conductor
    
    The spacing (i.e. cuts into the ground plane) can be controlled via:
        * side_gap - If this is zero, then the gap on the sides of the capacitor is calculated via a 50ohm impedance CPW line. Otherwise,
                     it is set via the given side gap.
        * init_pad - This adds spacing to the ground plane on the feed lines. This is useful when the diagonal section is steep (e.g.
                     when len_diag is zero) to ensure that the ground plane does not intersect with the main capacitor conductors. The
                     ground plane spacing typically starts to change on meeting with the LD or LF sections (see below). If init_pad > 0,
                     it starts from an earlier point.
        * side_gap_uniform - If True, the spacing on the sides are made uniform by taking the maximum extent of either the gap around
                             the inlet feed or the main capacitor body with the fingers.
    
    The positioning can be done dynamically via:
        * pin_inputs=Dict(start_pin=Dict(component=f'...',pin='...'), end_pin=Dict(component=f'...',pin='...')) - Specifying start and end
          positions via a component pins
    The resulting capacitor is right in the centre. This class ignores pos_x, pos_y and orientation...
        
    Pins:
        There are two pins on the capacitor at either end
        The pins attach directly to the built in lead length and only needs a width defined
        * cpw_width - center trace width of the trace lead line and cap fingers

    Note that the pins are not drawn to ensure good compatibility with routing/wiring constructs.

    Sketch:
        Below is a sketch of the capacitor
        ::

        @@@@@   |   @@@@@                               When setting init_pad > 0:
        @@@@@   |   @@@@@                                   @@@@@   |   @@@@@
        @@@@@   |   @@@@@                                   @@@@@   |   @@@@@
        @@@@@   |   @@@@@   @   = Ground Plane              @@      |   IP @@
        @@@@@   -   LD  @   -,| = Conductors                @      ---  LD  @
        @@@@   ---  LD  @                                   @    -------LF  @
        @@@   ----- LD  @   LD = len_diag                   @    -------LF  @ 
        @    -------LF  @                                   @    | | | |    @
        @    -------LF  @   LF = len_flat                   @    |||||||    @   IP = init_pad
        @ FL | | | |FLG @   FLG = fing_len_gap              @    |||||||    @
        @ FL |||||||    @                                   @     | | |     @ 
        @ FL |||||||SSSS@   S = side_gap                    @    -------LF  @
        @     | | | FLG @   FL = fing_len                   @    -------LF  @ 
        @    -------LF  @                                   @      ---  LD  @
        @    -------LF  @   FP = Front Pad                  @@      |   IP @@
        @@@   ----- LD  @                                   @@@@@   |   @@@@@
        @@@@   ---  LD  @                                   @@@@@   |   @@@@@
        @@@@@   -   LD  @
        @@@@@   |       @   
        @@@@@   |       @
        @@@@@   |       @
        @@@@@   |   @@@@@
        @@@@@   |   @@@@@

    .. image::
        Cap3Interdigital.png

    .. meta::
        Cap 3 Interdigital

    Default Options:
        * trace_width: '10um
        * dist_extend='50um'
        * cpw_width='10um'
        * len_diag='5um'
        * len_flat='5um'
        * fing_len='10um'
        * fing_len_gap='1um'
        * fing_wid='2um'
        * fing_wid_gap='1um'
        * N_total=5
        * larger_first=True
        * side_gap_uniform=False
        * side_gap='0um'
        * init_pad='0um'
    """

    #  Define structure functions

    default_options = Dict(dist_extend='50um',
                           cpw_width='10um',
                           len_diag='5um',
                           len_flat='5um',
                           fing_len='10um',
                           fing_len_gap='1um',
                           fing_wid='2um',
                           fing_wid_gap='1um',
                           N_total=5,
                           larger_first=True,
                           side_gap_uniform=False,
                           side_gap='0um',
                           init_pad='0um')
    """Default drawing options"""

    TOOLTIP = """Create a three finger planar capacitor with a ground pocket cuttout."""

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################

        start_point = self.design.components[self.options.pin_inputs.start_pin.component].pins[self.options.pin_inputs.start_pin.pin]
        startPt = start_point['middle']
        p.pos_x = startPt[0]
        p.pos_y = startPt[1]

        end_point = self.design.components[self.options.pin_inputs.end_pin.component].pins[self.options.pin_inputs.end_pin.pin]
        endPt = end_point['middle']
        p.end_x = endPt[0]
        p.end_y = endPt[1]

        pad1, pad2, padGap, pin1, pin2 = CapacitorInterdigital._draw_capacitor(p, self._design, True)

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(pad1=pad1, pad2=pad2),
                           layer=p.layer)

        #subtracts out ground plane on the layer it's on
        self.add_qgeometry('poly',
                           dict(padGap=padGap),
                           subtract=True,
                           layer=p.layer)

        # Generates its own pins
        self.add_pin('a', pin1.coords[::-1], width=p.cpw_width)
        self.add_pin('b', pin2.coords[::-1], width=p.cpw_width)

class CapacitorGap(QComponent):
    """Creates a gap capacitor on a CPW with an optional bisectional ground plane.

    Inherits QComponent class.

    Capacitor Metal Geometry and Ground Cutout Pocket:
        * len_diag - Length of the staggered region that bridges from CPW to the capacitor
        * len_flat - Length of the flat region before starting onto the fingers
        * cap_width - Width of the main capacitor
        * cap_gap  - Distance between the two conductors of the capacitor
        * gnd_width - Width of ground plane that bisects the two conductors of the capacitor (can be zero)
        * offset_lead1 - Offsets the first lead (positive being to the right when facing into the capacitor) along the capacitor.
        * offset_lead2 - Offsets the second lead (positive being to the right when facing into the capacitor) along the capacitor.

    The spacing (i.e. cuts into the ground plane) can be controlled via:
        * side_gap - If this is zero, then the gap on the sides of the capacitor is calculated via a 50ohm impedance CPW line. Otherwise,
                     it is set via the given side gap.
        * init_pad - This adds spacing to the ground plane on the feed lines. This is useful when the diagonal section is steep (e.g.
                     when len_diag is zero) to ensure that the ground plane does not intersect with the main capacitor conductors. The
                     ground plane spacing typically starts to change on meeting with the LD or LF sections (see below). If init_pad > 0,
                     it starts from an earlier point.

    As usual, the positioning can be done dynamically as a vector given by the supplied parameters: (pos_x,pos_y) to (end_x,end_y)

    Pins:
        There are two pins on the capacitor at either end
        The pins attach directly to the built in lead length and only needs a width defined
        * trace_width - center trace width of the trace lead line and cap fingers

    Sketch:
        Below is a sketch of the capacitor
        ::

            @@@@   |   @@@@                                         When setting len_diag=0:
            @@@@   |   @@@@            @   = Ground Plane               @@@@@   |   @@@@@
            @@@@   |   @@@@            -,| = Conductors                 @@@@@   |   @@@@@
            @@@@   -    @@@ LD                                          @@      |   IP @@
            @@@   ---    @@ LD         LD = len_diag                    @      ---  LD  @
            @@   -----    @ LD                                          @    -------LF  @
            @   -------   @ LF         LF = len_flat                    @    -------LF  @     IP = init_pad
            @   -------   @ LF                                          @           CG  @
            @   WWWWWWW   @ CG         W = cap_width                    @           CG  @
            @             @ CG                                       GW @@@@@@@@@@@@@@@@@
         GW @@@@@@@@@@@@@@@ CG         CG = cap_gap                     @           CG  @
         GW @@@@@@@@@@@@@@@ CG         GW = gnd_width                   @           CG  @
         GW @@@@@@@@@@@@@@@ CG                                          @    -------LF  @
            @             @ CG                                          @    -------LF  @
            @             @ CG                                          @      ---  LD  @
            @   -------   @ LF                                          @@      |   IP @@
            @   -------   @ LF                                          @@@@@   |   @@@@@
            @@   -----    @ LD                                          @@@@@   |   @@@@@
            @@@   ---    @@ LD
            @@@@   -    @@@ LD
            @@@@   |   @@@@
            @@@@   |   @@@@
            @@@@   |   @@@@

    .. image::
        Cap3Interdigital.png

    .. meta::
        Cap 3 Interdigital

    Default Options:
        * trace_width: '10um
        * pos_x='0um',pos_y='0um'
        * end_x='50um',end_y='0um'
        * len_diag='5um'
        * len_flat='5um'
        * cap_gap='3um'
        * gnd_width='1um'
        * side_gap='0um'
        * init_pad='0um'
        * offset_lead1='0um'
        * offset_lead2='0um'
    """

    #  Define structure functions

    default_options = Dict(pos_x='0um',pos_y='0um',
                           end_x='50um',end_y='0um',
                           cpw_width='10um',
                           cap_width='20um',
                           len_diag='5um',
                           len_flat='5um',
                           cap_gap='3um',
                           gnd_width='1um',
                           side_gap='0um',
                           init_pad='0um',
                           offset_lead1='0um',
                           offset_lead2='0um')
    """Default drawing options"""

    TOOLTIP = """Create a three finger planar capacitor with a ground pocket cuttout."""

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################

        pad1, pad2, padGap1, padGap2, pin1, pin2, polyBndMid = CapacitorGap._draw_capacitor(p, self._design)

        self.bound_poly_mid = polyBndMid

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(pad1=pad1,
                                pad2=pad2),
                           layer=p.layer)

        #subtracts out ground plane on the layer it's on
        self.add_qgeometry('poly',
                           dict(padGap1=padGap1,
                                padGap2=padGap2),
                           subtract=True,
                           layer=p.layer)

        # Generates its own pins
        self.add_pin('a', pin1.coords[::-1], width=p.cpw_width)
        self.add_pin('b', pin2.coords[::-1], width=p.cpw_width) 
    
    @staticmethod
    def _draw_capacitor(p, design):
        # Make the shapely polygons for the main cap structure
        len_trace = np.sqrt((p.end_x-p.pos_x)**2+(p.end_y-p.pos_y)**2)
        pad1 = [(0, p.cpw_width*0.5),
                (0, -p.cpw_width*0.5),
                (len_trace*0.5-p.cap_gap*0.5-p.len_flat-p.len_diag, -p.cpw_width*0.5),
                (len_trace*0.5-p.cap_gap*0.5-p.len_flat, -p.cap_width*0.5),
                (len_trace*0.5-p.cap_gap*0.5, -p.cap_width*0.5),
                (len_trace*0.5-p.cap_gap*0.5, p.cap_width*0.5),
                (len_trace*0.5-p.cap_gap*0.5-p.len_flat, p.cap_width*0.5),
                (len_trace*0.5-p.cap_gap*0.5-p.len_flat-p.len_diag, p.cpw_width*0.5),
                (0, p.cpw_width*0.5)]
        pad2 = np.array(pad1)
        pad2[:,0] = len_trace - pad2[:,0]
        #
        #Handle lead offsets...
        pad1 = np.array(pad1)
        pad1[[0,1,2,7,8],1] -= p.offset_lead1
        # if np.abs(p.len_diag) < 1e-12 and np.abs(-p.cap_width*0.5+p.cpw_width*0.5 - p.offset_lead1) < 1e-12:     
        #     pad1 = np.delete(pad1, 7, 0)
        # elif np.abs(p.cap_width*0.5-p.cpw_width*0.5 - p.offset_lead1) < 1e-12:
        #     pad1 = np.delete(pad1, 2, 0)
        pad2[[0,1,2,7,8],1] -= p.offset_lead2
        # if np.abs(p.len_diag) < 1e-12 and np.abs(-p.cap_width*0.5+p.cpw_width*0.5 - p.offset_lead2) < 1e-12:     
        #     pad2 = np.delete(pad2, 7, 0)
        # elif np.abs(p.cap_width*0.5-p.cpw_width*0.5 - p.offset_lead2) < 1e-12:
        #     pad2 = np.delete(pad2, 2, 0)
        #
        units = QUtilities.get_units(design)
        cpwP = CpwParams.fromQDesign(design)
        gap_cpw_line = cpwP.get_gap_from_width(p.cpw_width*units)/units
        if p.side_gap == 0:
            gap_cpw_cap = cpwP.get_gap_from_width(p.cap_width*units)/units
        else:
            gap_cpw_cap = p.side_gap
        #
        padGap1 = pad1*1.0
        padGap2 = pad2*1.0
        padGap1[[-2,-1,0],1] += gap_cpw_line
        padGap1[[1,2],1] -= gap_cpw_line
        padGap1[[5,6],1] += gap_cpw_cap
        padGap1[[3,4],1] -= gap_cpw_cap
        padGap2[[-2,-1,0],1] += gap_cpw_line
        padGap2[[1,2],1] -= gap_cpw_line
        padGap2[[5,6],1] += gap_cpw_cap
        padGap2[[3,4],1] -= gap_cpw_cap
        if p.len_diag == 0:
            padGap1[[2,3,6,7],0] -= p.init_pad
            padGap2[[2,3,6,7],0] += p.init_pad
        else:
            padGap1[[2,7],0] -= p.init_pad
            padGap1[[2,7],0] += p.init_pad
        padGap1[4:6,0] += 0.5*(p.cap_gap-p.gnd_width)
        padGap2[4:6,0] -= 0.5*(p.cap_gap-p.gnd_width)
        #
        pad2 = pad2[::-1]
        padGap2 = padGap2[::-1]
        #
        pin1 = pad1[:2]
        pin2 = pad2[[-2,-1]]

        polyBndMid = [(len_trace*0.5-p.cap_gap*0.5, -p.cap_width*0.5-gap_cpw_cap),
                      (len_trace*0.5+p.cap_gap*0.5, -p.cap_width*0.5-gap_cpw_cap),
                      (len_trace*0.5+p.cap_gap*0.5, p.cap_width*0.5+gap_cpw_cap),
                      (len_trace*0.5-p.cap_gap*0.5, p.cap_width*0.5+gap_cpw_cap)]
        polyBndMid = shapely.Polygon(polyBndMid)

        pad1 = shapely.Polygon(pad1)
        pad2 = shapely.Polygon(pad2)
        padGap1 = shapely.Polygon(padGap1)
        padGap2 = shapely.Polygon(padGap2)
        pin1 = shapely.LineString(pin1)
        pin2 = shapely.LineString(pin2)

        polys = [pad1, pad2, padGap1, padGap2, pin1, pin2, polyBndMid]
        polys = draw.rotate(polys, np.arctan2(p.end_y-p.pos_y,p.end_x-p.pos_x), origin=(0, 0), use_radians=True)
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [pad1, pad2, padGap1, padGap2, pin1, pin2, polyBndMid] = polys

        return pad1, pad2, padGap1, padGap2, pin1, pin2, polyBndMid


class CapacitorProngPin(QComponent):
    """Creates a forked gap capacitorthat wraps around the target lead. Usually used for Xmons or any line couplings.

    Inherits QComponent class.

    Capacitor Metal Geometry and Ground Cutout Pocket:
        * prong_width  - Width of the two prongs
        * prong_length - Length of the two prongs
        * pin_gap      - Capacitor gap from the fork to the target pin
        * pin_gap_side - Capacitor gap of side prongs to the sides of the target pin
        * pad_thickness - Length of the Capacitor fork pad

    The structure automatically cuts a box around the entire fork; further spacing (i.e. cuts into the ground plane) can be controlled via:
        * gap_side  - Spacing to ground plane on the outer sides of the fork
        * gap_front - Spacing to ground plane in the direction of the target pin
        * gap_back - Spacing to ground plane away from the direction of the target pin

    The positioning can be done dynamically via:
        * pin_inputs=Dict(start_pin=Dict(component=f'...',pin='...')) - Specifying the target component pin

    Pins:
        There is one pin 'a' to link to the fork (there is no trace drawn in this structure intrinsically). Only needs a width defined
        * trace_width - center trace width of the trace that attaches to this capacitor.

    Sketch:
        Below is a sketch of the capacitor
        ::

            @@@@@@@@@@@@@@@@@@@     @  = Ground Plane
            @              GB @     #  = Target pin
            @   _____x_____GB @     x = pin location with width trace_width
            @ T|           |  @     W = prong_width
            @ T|  _______  |  @     L = prong_length
            @ L| |   P   | |GS@     P = pin_gap
            @ L|W|  ###  | |  @     S = pin_gap_side
            @ L|_|SS###  |_|  @     T = pad_thickness
            @       ###   GF  @     BG = gap_back
            @       ###   GF  @     GF = gap_front
            @@@@@@  ###  @@@@@@     GS = gap_side

    .. image::
        Cap3Interdigital.png

    .. meta::
        Cap 3 Interdigital

    Default Options:
        * prong_width='4um'
        * prong_length='10um'
        * pin_gap='6um'
        * pin_gap_side='4um'
        * pad_thickness='10um'
        * gap_side='5um'
        * gap_front='5um'
        * gap_back='5um'
        * trace_width='10um'
    """

    #  Define structure functions

    default_options = Dict(prong_width='4um',
                           prong_length='10um',
                           pin_gap='6um',
                           pin_gap_side='4um',
                           pad_thickness='10um',
                           gap_side='5um',
                           gap_front='5um',
                           gap_back='5um',
                           trace_width='10um')
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
        startPtNorm = start_point['normal']

        wid = start_point['width']

        pad = [
              (p.pin_gap-p.prong_length, wid*0.5+p.pin_gap_side),
              (p.pin_gap-p.prong_length, wid*0.5+p.pin_gap_side+p.prong_width),
              (p.pin_gap+p.pad_thickness, wid*0.5+p.pin_gap_side+p.prong_width),
              (p.pin_gap+p.pad_thickness, -wid*0.5-p.pin_gap_side-p.prong_width),
              (p.pin_gap-p.prong_length, -wid*0.5-p.pin_gap_side-p.prong_width),
              (p.pin_gap-p.prong_length, -wid*0.5-p.pin_gap_side),
              (p.pin_gap, -wid*0.5-p.pin_gap_side),
              (p.pin_gap, wid*0.5+p.pin_gap_side)]

        gap = [
              (p.pin_gap+p.pad_thickness+p.gap_back, -wid*0.5-p.pin_gap_side-p.prong_width-p.gap_side),
              (p.pin_gap+p.pad_thickness+p.gap_back, wid*0.5+p.pin_gap_side+p.prong_width+p.gap_side),
              (p.pin_gap-p.prong_length-p.gap_front, wid*0.5+p.pin_gap_side+p.prong_width+p.gap_side),
              (p.pin_gap-p.prong_length-p.gap_front, -wid*0.5-p.pin_gap_side-p.prong_width-p.gap_side)]

        pad = shapely.Polygon(pad[::-1])
        gap = shapely.Polygon(gap)
        pin = shapely.LineString([(p.pin_gap+p.pad_thickness, wid*0.5+p.pin_gap_side+p.prong_width),
                                  (p.pin_gap+p.pad_thickness, -wid*0.5-p.pin_gap_side-p.prong_width)])

        polys = [pad, gap, pin]
        polys = draw.rotate(polys, np.arctan2(startPtNorm[1], startPtNorm[0]), origin=(0, 0), use_radians=True)
        polys = draw.translate(polys, *startPt)
        [pad, gap, pin] = polys

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(pad=pad),
                           layer=p.layer)

        #subtracts out ground plane on the layer it's on
        self.add_qgeometry('poly',
                           dict(padGap=gap),
                           subtract=True,
                           layer=p.layer)

        # Generates its own pin
        self.add_pin('a', pin.coords, width=p.trace_width)
