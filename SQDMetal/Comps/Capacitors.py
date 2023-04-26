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
                           side_gap='0um',
                           init_pad='0um')
    """Default drawing options"""

    TOOLTIP = """Create a three finger planar capacitor with a ground pocket cuttout."""

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################

        # Make the shapely polygons for the main cap structure
        len_fings_plus_gap = p.fing_len + p.fing_len_gap
        cap_width = p.N_total*p.fing_wid + (p.N_total-1)*p.fing_wid_gap
        len_trace = np.sqrt((p.end_x-p.pos_x)**2+(p.end_y-p.pos_y)**2)
        pad1 = [(0, p.cpw_width*0.5),
                (0, -p.cpw_width*0.5),
                (len_trace*0.5-len_fings_plus_gap*0.5-p.len_flat-p.len_diag, -p.cpw_width*0.5),
                (len_trace*0.5-len_fings_plus_gap*0.5-p.len_flat, -cap_width*0.5),
                (len_trace*0.5-len_fings_plus_gap*0.5, -cap_width*0.5),
                (len_trace*0.5-len_fings_plus_gap*0.5, cap_width*0.5),
                (len_trace*0.5-len_fings_plus_gap*0.5-p.len_flat, cap_width*0.5),
                (len_trace*0.5-len_fings_plus_gap*0.5-p.len_flat-p.len_diag, p.cpw_width*0.5),
                (0, p.cpw_width*0.5)]
        pad2 = np.array(pad1)
        pad2[:,0] = len_trace - pad2[:,0]
        pad2 = pad2[::-1]
        #
        units = QUtilities.get_units(self._design)
        cpwP = CpwParams.fromQDesign(self._design)
        gap_cpw_line = cpwP.get_gap_from_width(p.cpw_width*units)/units
        if p.side_gap == 0:
            gap_cpw_cap = cpwP.get_gap_from_width(cap_width*units)/units
        else:
            gap_cpw_cap = p.side_gap
        #
        padGap1 = np.array(pad1)
        padGap1[[-2,-1,0],1] += gap_cpw_line
        padGap1[[1,2],1] -= gap_cpw_line
        padGap1[[5,6],1] += gap_cpw_cap
        padGap1[[3,4],1] -= gap_cpw_cap
        if p.len_diag == 0:
            padGap1[[2,3,6,7],0] -= p.init_pad
        else:
            padGap1[[2,7],0] -= p.init_pad
        padGap1[4:6,0] += 0.5*(len_fings_plus_gap)
        padGap2 = padGap1*1.0
        padGap2[:,0] = len_trace - padGap2[:,0]
        padGap2 = padGap2[::-1]
        #
        pin1 = pad1[:2]
        pin2 = pad2[[-2,-1]]
        #
        pad1 = shapely.Polygon(pad1)
        pad2 = shapely.Polygon(pad2)
        padGap = shapely.unary_union([shapely.Polygon(padGap1), shapely.Polygon(padGap2)])
        pin1 = shapely.LineString(pin1)
        pin2 = shapely.LineString(pin2)
        #Do the fingers
        cur_x = len_trace*0.5-len_fings_plus_gap*0.5
        cur_y = cap_width*0.5
        poly_fings_1 = []
        poly_fings_2 = []
        flip = False
        for m in range(p.N_total):
            if flip:
                poly_fings_2.append(shapely.Polygon([[cur_x+p.fing_len_gap, cur_y], [cur_x+p.fing_len_gap, cur_y-p.fing_wid], [cur_x+p.fing_len_gap+p.fing_len, cur_y-p.fing_wid], [cur_x+p.fing_len_gap+p.fing_len, cur_y]]))
            else:
                poly_fings_1.append(shapely.Polygon([[cur_x, cur_y], [cur_x, cur_y-p.fing_wid], [cur_x+p.fing_len, cur_y-p.fing_wid], [cur_x+p.fing_len, cur_y]]))
            flip = not flip
            cur_y -= p.fing_wid + p.fing_wid_gap
        pad1 = shapely.unary_union([pad1] + poly_fings_1)
        pad2 = shapely.unary_union([pad2] + poly_fings_2)

        polys = [pad1, pad2, padGap, pin1, pin2]
        polys = draw.rotate(polys, np.arctan2(p.end_y-p.pos_y,p.end_x-p.pos_x), origin=(0, 0), use_radians=True)
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [pad1, pad2, padGap, pin1, pin2] = polys        

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(pad1=pad1, pad2=pad2),
                           layer=p.layer)

        #subtracts out ground plane on the layer its on
        self.add_qgeometry('poly',
                           dict(padGap=padGap),
                           subtract=True,
                           layer=p.layer)

        # Generates its own pins
        self.add_pin('a', pin1.coords[::-1], width=p.cpw_width)
        self.add_pin('b', pin2.coords[::-1], width=p.cpw_width) 

class CapacitorGap(QComponent):
    """Creates a gap capacitor on a CPW with an optional bisectional ground plane.
    The width of the fingers is determined by cpw_width.

    Inherits QComponent class.

    Capacitor Metal Geometry and Ground Cutout Pocket:
        * len_diag - Length of the staggered region that bridges from CPW to the capacitor
        * len_flat - Length of the flat region before starting onto the fingers
        * cap_width - Width of the main capacitor
        * cap_gap  - Distance between the two conductors of the capacitor
        * gnd_width - Width of ground plane that bisects the two conductors of the capacitor (can be zero)

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
                           init_pad='0um')
    """Default drawing options"""

    TOOLTIP = """Create a three finger planar capacitor with a ground pocket cuttout."""

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################

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
        pad2 = pad2[::-1]
        #
        units = QUtilities.get_units(self._design)
        cpwP = CpwParams.fromQDesign(self._design)
        gap_cpw_line = cpwP.get_gap_from_width(p.cpw_width*units)/units
        if p.side_gap == 0:
            gap_cpw_cap = cpwP.get_gap_from_width(p.cap_width*units)/units
        else:
            gap_cpw_cap = p.side_gap
        #
        padGap1 = np.array(pad1)
        padGap1[[-2,-1,0],1] += gap_cpw_line
        padGap1[[1,2],1] -= gap_cpw_line
        padGap1[[5,6],1] += gap_cpw_cap
        padGap1[[3,4],1] -= gap_cpw_cap
        if p.len_diag == 0:
            padGap1[[2,3,6,7],0] -= p.init_pad
        else:
            padGap1[[2,7],0] -= p.init_pad
        padGap1[4:6,0] += 0.5*(p.cap_gap-p.gnd_width)
        padGap2 = padGap1*1.0
        padGap2[:,0] = len_trace - padGap2[:,0]
        padGap2 = padGap2[::-1]
        #
        pin1 = pad1[:2]
        pin2 = pad2[[-2,-1]]
        #
        pad1 = shapely.Polygon(pad1)
        pad2 = shapely.Polygon(pad2)
        padGap1 = shapely.Polygon(padGap1)
        padGap2 = shapely.Polygon(padGap2)
        pin1 = shapely.LineString(pin1)
        pin2 = shapely.LineString(pin2)

        polys = [pad1, pad2, padGap1, padGap2, pin1, pin2]
        polys = draw.rotate(polys, np.arctan2(p.end_y-p.pos_y,p.end_x-p.pos_x), origin=(0, 0), use_radians=True)
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [pad1, pad2, padGap1, padGap2, pin1, pin2] = polys

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(pad1=pad1,
                                pad2=pad2),
                           layer=p.layer)

        #subtracts out ground plane on the layer its on
        self.add_qgeometry('poly',
                           dict(padGap1=padGap1,
                                padGap2=padGap2),
                           subtract=True,
                           layer=p.layer)

        # Generates its own pins
        self.add_pin('a', pin1.coords[::-1], width=p.cpw_width)
        self.add_pin('b', pin2.coords[::-1], width=p.cpw_width) 