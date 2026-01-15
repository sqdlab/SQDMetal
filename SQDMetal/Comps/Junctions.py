# -*- coding: utf-8 -*-

# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0
# Author: Prasanna Pakkiam
# Creation Date: 24/05/2023
# Description: Collection of classes to draw Josephson Junctions.

from qiskit_metal import draw
from qiskit_metal.toolbox_python.attr_dict import Dict
from qiskit_metal.qlibrary.core import QComponent
import numpy as np
import shapely
# for Manhattan JJ
from SQDMetal.Utilities.ShapelyEx import ShapelyEx

class JunctionDolan(QComponent):
    """Create a Dolan Bridge Josephson Junction

    Inherits QComponent class.

    The Dolan Bridge consists of a T-Section followed by a Fork Section.

    Dolan Bridge Josephson Junction Metal Geometry (no Ground Cutout):
        * bridge_gap - Gap between the T-Section and the Fork Section
        * finger_length - Length of the thin section on the end of the Fork Section
        * finger_width  - Width of the thin section on the end of the Fork Section
        * squid_width  - Overall width of the SQUID loop
        * stem_width  - Width of the leads connecting to the Josephson Junction
        * prong_width  - Length of the thicker section on the end of the Fork Section
        * prong_length  - Width of the thicker section on the end of the Fork Section
        * t_pad_width  - Thickness of the T-Section pad connecting to the leads
        * t_finger_width - Thicknes of the finger of the T-section pad connecting to the leads
        * t_pad_length - Length of the central pad of the T-section 
        * fork_pad_size  - Thickness of the Fork Section pad connecting to the leads
        * t_pad_extra - Extra length added to the sides of the T-Section to ensure proper overlap with fingers
        * reverse    - Default False. If True, the direction of the prongs is flipped by 180 degrees (that is, the prongs will
                       be near the first pin)

    As usual, the positioning can be done dynamically as a vector given by the supplied parameters: (pos_x,pos_y) to (end_x,end_y)
        
    Pins:
        There are pins given on either side to help position the bandage. They are called 't' and 'f' for the T-pad and
        fork respectively. Pin width is stem_width.

    Sketch:
        Below is a sketch of the Josephson Junction Shadow Evaporation masking template (there is no ground cut-out)
        ::

            SW...................SW
                      | |
                      |S|                      SW  = squid_width
             _________| |_________             S   = stem_width 
            |      _________      |     FPS    FPS = fork_pad_size
            |     |         |PW.PW|     PL     PL  = prong_length
            |_   _|         |_   _|     PL     PW  = prong_width
              | |             | |       FW  
              | |             |F|       FW     F  = finger_width
              |_|             |_|       FW     FW = finger_length
                                        BG
                  |<--TPL-->|           BG     BG = bridge_gap
                   _________            BG
         _________|    ^    |_________  BG                 
        |_________   TPW     _________| TFW    TFW = t_finger_width  
                  |___ | ___|                  TPL = t_pad_length
        <-->          | |         <-->          TPW = t_pad_width
         TPE          |S|          TPE           TPE = t_pad_extra
                      | |
    .. image::
        Cap3Interdigital.png

    .. meta::
        Dolan Bridge Josephson Junction

    Default Options:
        * trace_width: '10um
        * pos_x='0um',pos_y='0um'
        * end_x='50um',end_y='0um'
        * bridge_gap='0.2um'
        * finger_length='1.75um'
        * finger_width='0.235um'
        * squid_width='5.4um'
        * stem_width='2um'
        * prong_width='0.5um'
        * prong_length='1.75um'
        * t_pad_width='0.3um'
        * t_finger_width='0.3um'
        * t_pad_length= '3um'
        * fork_pad_size='0.5um'
        * t_pad_extra='0.0um'
    """

    default_options = Dict(pos_x='0um',pos_y='0um',
                           end_x='50um',end_y='0um',
                           bridge_gap='0.2um',
                           finger_length='1.75um',
                           finger_width='0.235um',
                           squid_width='5.4um',
                           stem_width='2um',
                           prong_width='0.5um',
                           prong_length='1.75um',
                           t_pad_width='0.3um',
                           t_finger_width='0.3um',
                           t_pad_length='3um',
                           fork_pad_size='0.5um',
                           t_pad_extra='0.0um',
                           reverse=False)

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################

        pad_T, pad_Fork, pin1, pin2, sim_JJ = JunctionDolan.draw_junction(p)

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(pad1=pad_T, pad_Fork=pad_Fork),
                           layer=p.layer)

        #subtracts out ground plane on the layer it's on
        # self.add_qgeometry('poly',
        #                    dict(padGap=padGap),
        #                    subtract=True,
        #                    layer=p.layer)

        self.add_qgeometry('junction',
                    {'design': sim_JJ},
                    layer=p.layer,
                    subtract=False,
                    width=p.squid_width+2*p.t_pad_extra)

        # Generates its own pins
        self.add_pin('t', pin1.coords[::-1], width=p.stem_width)
        self.add_pin('f', pin2.coords[::-1], width=p.stem_width)

    @staticmethod
    def draw_junction(p):
        struct_width = p.t_pad_width/2+ p.t_finger_width/2 + p.bridge_gap + p.finger_length + p.prong_length + p.fork_pad_size
        len_comp = np.sqrt((p.end_x-p.pos_x)**2+(p.end_y-p.pos_y)**2)
        len_stem = (len_comp - struct_width)/2
        t_finger_length=(p.squid_width-p.t_pad_length)/2

        #The T-Section and Stem
        pad_T = [
                (-p.t_pad_width*0.5, p.stem_width*0.5), #1
                (-len_stem-p.t_pad_width*0.5, p.stem_width*0.5), #2
                (-len_stem-p.t_pad_width*0.5, -p.stem_width*0.5), #3
                (-p.t_pad_width*0.5, -p.stem_width*0.5), #4
                (-p.t_pad_width*0.5,-p.t_pad_length*0.5), #4'
                (-p.t_finger_width*0.5,-p.t_pad_length*0.5), #5
                (-p.t_finger_width*0.5,-p.t_pad_length*0.5-t_finger_length-p.t_pad_extra), #6
                (p.t_finger_width*0.5,-p.t_pad_length*0.5-t_finger_length-p.t_pad_extra), #7
                (p.t_finger_width*0.5,-p.t_pad_length*0.5), #8
                (p.t_pad_width*0.5,-p.t_pad_length*0.5), #9
                #
                (p.t_pad_width*0.5,p.t_pad_length*0.5), #10
                (p.t_finger_width*0.5,p.t_pad_length*0.5), #11
                (p.t_finger_width*0.5,p.t_pad_length*0.5+t_finger_length+p.t_pad_extra), #12
                (-p.t_finger_width*0.5,p.t_pad_length*0.5+t_finger_length+p.t_pad_extra), #13
                (-p.t_finger_width*0.5,p.t_pad_length*0.5), #14
                (-p.t_pad_width*0.5,p.t_pad_length*0.5), #15
                ] 
    
        #The Fork Section and Stem
        pad_Fork = [
                   (p.finger_length+p.prong_length+p.fork_pad_size+len_stem, p.stem_width*0.5),
                   (p.finger_length+p.prong_length+p.fork_pad_size, p.stem_width*0.5),
                   (p.finger_length+p.prong_length+p.fork_pad_size, p.squid_width*0.5),
                   (p.finger_length+p.prong_length, p.squid_width*0.5),
                   (p.finger_length, p.squid_width*0.5),
                   (p.finger_length, (p.squid_width-p.prong_width+p.finger_width)*0.5),
                   (0, (p.squid_width-p.prong_width+p.finger_width)*0.5),
                   (0, (p.squid_width-p.prong_width-p.finger_width)*0.5),
                   (p.finger_length, (p.squid_width-p.prong_width-p.finger_width)*0.5),
                   (p.finger_length, p.squid_width*0.5-p.prong_width),
                   (p.finger_length+p.prong_length, p.squid_width*0.5-p.prong_width),
                   #
                   (p.finger_length+p.prong_length, -p.squid_width*0.5+p.prong_width),
                   (p.finger_length, -p.squid_width*0.5+p.prong_width),
                   (p.finger_length, (-p.squid_width+p.prong_width+p.finger_width)*0.5),
                   (0, (-p.squid_width+p.prong_width+p.finger_width)*0.5),
                   (0, (-p.squid_width+p.prong_width-p.finger_width)*0.5),
                   (p.finger_length, (-p.squid_width+p.prong_width-p.finger_width)*0.5),
                   (p.finger_length, -p.squid_width*0.5),
                   (p.finger_length+p.prong_length, -p.squid_width*0.5),
                   (p.finger_length+p.prong_length+p.fork_pad_size, -p.squid_width*0.5),
                   (p.finger_length+p.prong_length+p.fork_pad_size, -p.stem_width*0.5),
                   (p.finger_length+p.prong_length+p.fork_pad_size+len_stem, -p.stem_width*0.5)]

        pad_T = np.array(pad_T)
        pad_T[:,0] += len_stem + p.t_pad_width/2
        pad_Fork = np.array(pad_Fork)
        pad_Fork[:,0] += p.t_finger_width*0.5+ p.bridge_gap + len_stem + p.t_pad_width/2
        
        sim_JJ = shapely.LineString([
            np.mean(pad_T[1:3],axis=0),
            np.mean([pad_Fork[-1], pad_Fork[0]],axis=0)
        ])

        pin1 = shapely.LineString(pad_T[1:3])
        pin2 = shapely.LineString([pad_Fork[-1], pad_Fork[0]])
        pad_T = shapely.Polygon(pad_T)
        pad_Fork = shapely.Polygon(pad_Fork)

        polys = [pad_T, pad_Fork, pin1, pin2, sim_JJ]
        if p.reverse:
            polys = draw.translate(polys, -len_comp, 0)
            polys = draw.rotate(polys, np.pi, origin=(0, 0), use_radians=True)
        polys = draw.rotate(polys, np.arctan2(p.end_y-p.pos_y,p.end_x-p.pos_x), origin=(0, 0), use_radians=True)
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [pad_T, pad_Fork, pin1, pin2, sim_JJ] = polys
        return pad_T, pad_Fork, pin1, pin2, sim_JJ

class JunctionSingleDolan(QComponent):
    """Create a single Dolan Bridge Josephson Junction

    Inherits QComponent class.

    Consists of a vertical prong with a finger, and a horizontal half-T.

    Dolan Bridge Josephson Junction Metal Geometry (no Ground Cutout):
        * bridge_gap - Gap between the half-T and the prong Section
        * finger_length - Length of the thin section on the end of the prong Section
        * finger_width  - Width of the thin section on the end of the prong Section
        * prong_width  - Length of the thicker section on the end of the prong Section
        * prong_length  - Width of the thicker section on the end of the prong Section
        * t_pad_width  - Thickness of the half-T pad connecting to the leads
        * t_pad_length  - Length of the half-T pad connecting to the leads
        * t_finger_width - Thickness of the finger of the half-T pad connecting to the leads
        * t_pad_extra - Extra length added to the sides of the T-Section to ensure proper overlap with fingers
        * reverse    - Default False. If True, the direction of the prongs is flipped by 180 degrees (that is, the prongs will
                       be near the first pin)

    As usual, the positioning can be done dynamically as a vector given by the supplied parameters: (pos_x,pos_y) to (end_x,end_y)
        
    Pins:
        There are pins given on either side to help position the bandage. They are called 't' and 'f' for the T-pad and
        fork respectively. Pin width is stem_width.

    The junction area is given by F x TFW (finger_width x t_finger_width)

    Sketch:
        Below is a sketch of the Josephson Junction Shadow Evaporation masking template (there is no ground cut-out)
        ::
                    'f' pin
                    |     |
             _______|<-SW>|                 
            |             |     SW   
            |      _______|     SW            
            |<PW->|             PL     PL  = prong_length
            |_   _|             PL     PW  = prong_width
              | |               FW  
              |F|               FW     F  = finger_width
              |_|               FW     FW = finger_length
                                BG     BG = bridge_gap
             TFL  <--TPL-->     BG
            <----> ________     BG
         _________|    ^  |     BG     TFL = t_finger_length                
        |_________    TPW |     TFW    TFW = t_finger_width  
                  |__  |  |            TPL = t_pad_length
        <-->        |     |            TPW = t_pad_width
         TPE        |<SW->|            TPE = t_pad_extra
                    |     |            SW  = stem_width
                    't' pin
                    
    .. image::

    .. meta::
        Fixed Frequency Dolan Bridge Josephson Junction

    Default Options:
        * pos_x='0um', pos_y='0um'
        * end_x='20um', end_y='0um'
        * dist_extend='40um'
        * bridge_gap='0.2um'
        * finger_length='1.75um'
        * finger_width='0.23um'
        * prong_width='0.5um'
        * prong_length='1.75um'
        * t_pad_width='0.6um'
        * t_finger_width='0.23um'
        * t_finger_length='1um'
        * t_pad_length= '1um'
        * t_pad_extra='0.1um'
        * stem_width=`1um`
    """
    default_options = Dict(pos_x='0um', pos_y='0um',
                           end_x='20um', end_y='0um', 
                           dist_extend='40um',
                           bridge_gap='0.2um',
                           finger_length='1.75um',
                           finger_width='0.230um',
                           prong_width='0.5um',
                           prong_length='1.75um',
                           t_pad_width='0.6um',
                           t_finger_width='0.230um',
                           t_finger_length='1um',
                           t_pad_length='3um',
                           t_pad_extra='0.1um',
                           stem_width='1um',
                           reverse=False)

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################

        pad_T, pad_Fork, pin1, pin2, sim_JJ = JunctionSingleDolan.draw_junction(p)

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(pad1=pad_T, pad_Fork=pad_Fork),
                           layer=p.layer)
        
        self.add_qgeometry('junction',
                           {'design': sim_JJ},
                           layer=p.layer,
                           subtract=False,
                           width=p.t_pad_length+p.t_pad_extra+p.prong_width)

        # Generates its own pins
        self.add_pin('t', pin1.coords[::-1], width=p.t_pad_length)
        self.add_pin('f', pin2.coords[::-1], width=p.prong_width)
    
    @staticmethod
    def draw_junction(p):
        struct_width = p.t_pad_width/2 + p.t_finger_width/2 + \
            p.bridge_gap + p.finger_length + p.prong_length
        len_comp = np.sqrt((p.end_x-p.pos_x)**2+(p.end_y-p.pos_y)**2)
        len_stem = (len_comp - struct_width)/2 # extension from start/end to start of junction structure

        # The T-Section and Stem
        pad_T = [
            (-p.t_pad_width*0.5, p.stem_width*0.5),  # 1
            (-len_stem-p.t_pad_width*0.5, p.stem_width*0.5),  # 2
            (-len_stem-p.t_pad_width*0.5, -p.stem_width*0.5),  # 3
            (-p.t_pad_width*0.5, -p.stem_width*0.5),  # 4
            (-p.t_pad_width*0.5, -p.t_pad_length*0.5),  # 4'
            (-p.t_finger_width*0.5, -p.t_pad_length*0.5),  # 5
            (-p.t_finger_width*0.5, -p.t_pad_length *
             0.5-p.t_finger_length-p.t_pad_extra),  # 6
            (p.t_finger_width*0.5, -p.t_pad_length *
             0.5-p.t_finger_length-p.t_pad_extra),  # 7
            (p.t_finger_width*0.5, -p.t_pad_length*0.5),  # 8
            (p.t_pad_width*0.5, -p.t_pad_length*0.5),  # 9
            (p.t_pad_width*0.5, p.stem_width*0.5),  # 10
        ]

        # The Fork Section and Stem
        p.squid_width = (p.t_finger_length + (p.t_pad_length/2))*2
        p.fork_pad_size = p.stem_width
        pad_Fork = [
            (p.finger_length+p.prong_length +
             p.fork_pad_size+len_stem, p.stem_width*0.5),
            (p.finger_length+p.prong_length, p.stem_width*0.5),
            (p.finger_length+p.prong_length, -
             p.squid_width*0.5+p.prong_width),
            (p.finger_length, -p.squid_width*0.5+p.prong_width),
            (p.finger_length, (-p.squid_width+p.prong_width+p.finger_width)*0.5),
            (0, (-p.squid_width+p.prong_width+p.finger_width)*0.5),
            (0, (-p.squid_width+p.prong_width-p.finger_width)*0.5),
            (p.finger_length, (-p.squid_width+p.prong_width-p.finger_width)*0.5),
            (p.finger_length, -p.squid_width*0.5),
            (p.finger_length+p.prong_length, -p.squid_width*0.5),
            (p.finger_length+p.prong_length +
             p.fork_pad_size, -p.squid_width*0.5),
            (p.finger_length+p.prong_length +
             p.fork_pad_size, -p.stem_width*0.5),
            (p.finger_length+p.prong_length+p.fork_pad_size+len_stem, -p.stem_width*0.5)]

        pad_T = np.array(pad_T)
        pad_T[:, 0] += len_stem + p.t_pad_width/2 # adjust in x
        pad_Fork = np.array(pad_Fork)
        pad_Fork[:, 0] += p.t_finger_width*0.5 + \
            p.bridge_gap + len_stem + p.t_pad_width/2

        sim_JJ = shapely.LineString([
            np.mean(pad_T[1:3], axis=0),
            np.mean([pad_Fork[-1], pad_Fork[0]], axis=0)
        ])

        # draw pins and shapes
        pin1 = shapely.LineString(pad_T[1:3])
        pin2 = shapely.LineString([pad_Fork[-1], pad_Fork[0]])
        pad_T = shapely.Polygon(pad_T)
        pad_Fork = shapely.Polygon(pad_Fork)

        # move into position relative to given coords pos_x, pos_y
        polys = [pad_T, pad_Fork, pin1, pin2, sim_JJ]
        if p.reverse:
            polys = draw.translate(polys, -len_comp, 0)
            polys = draw.rotate(polys, np.pi, origin=(0, 0), use_radians=True)
        polys = draw.rotate(polys, np.arctan2(
            p.end_y-p.pos_y, p.end_x-p.pos_x), origin=(0, 0), use_radians=True)
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [pad_T, pad_Fork, pin1, pin2, sim_JJ] = polys
        return pad_T, pad_Fork, pin1, pin2, sim_JJ

class JunctionSingleDolanPinStretch(QComponent):
    """
    Create a single Dolan Bridge Josephson Junction

    Inherits QComponent class.

    Consists of a vertical prong with a finger, and a horizontal half-T.

    Dolan Bridge Josephson Junction Metal Geometry(no Ground Cutout):
        * bridge_gap - Gap between the half-T and the prong Section
        * dist_extend = distance to extend from start pin
        * finger_length - Length of the thin section on the end of the prong Section
        * finger_width - Width of the thin section on the end of the prong Section
        * prong_width - Length of the thicker section on the end of the prong Section
        * prong_length - Width of the thicker section on the end of the prong Section
        * t_pad_width - Thickness of the half-T pad connecting to the leads
        * t_pad_length - Length of the half-T pad connecting to the leads
        * t_finger_width - Thickness of the finger of the half-T pad connecting to the leads
        * t_pad_extra - Extra length added to the sides of the T-Section to ensure proper overlap with fingers
        * reverse    - Default False. If True, the direction of the prongs is flipped by 180 degrees (that is , the prongs will
                       be near the first pin)

    The positioning can be done dynamically via:
        * pin_inputs=Dict(start_pin=Dict(component=f'...',pin='...')) - Specifying start position via a component pin
        * dist_extend - Distance upon to stretch away from the start pin.
    The resulting Josephson junction is right in the centre. This class ignores pos_x, pos_y and orientation...
        
    Pins:
        There are pins given on either side to help position the bandage. They are called 't' and 'f' for the T-pad and
        fork respectively. Pin width is stem_width.

    The junction area is given by F x TFW(finger_width x t_finger_width)

    Sketch:
        Below is a sketch of the Josephson Junction Shadow Evaporation masking template (there is no ground cut-out)
        : :
                    'f' pin
                    |     |                                     ^
             _______|<-SW>|                                     |
            |             |     SW                         dist_extend
            |      _______|     SW                              |
            |<PW->|             PL     PL  = prong_length       |
            |_   _|             PL     PW  = prong_width        |
              | |               FW                              |
              |F|               FW     F  = finger_width        |
              |_|               FW     FW = finger_length       |
                                BG     BG = bridge_gap          |
             TFL  <--TPL-->     BG                              |
            <----> ________     BG                              |
         _________|    ^  |     BG     TFL = t_finger_length    |            
        |_________    TPW |     TFW    TFW = t_finger_width     |
                  |__  |  |            TPL = t_pad_length       |
        <-->        |     |            TPW = t_pad_width        |
         TPE        |<SW->|            TPE = t_pad_extra        |
                    |     |            SW  = stem_width         |
                    't' pin
                    
    .. image::

    .. meta::
        Fixed Frequency Dolan Bridge Josephson Junction

    Default Options:
        * dist_extend='40um'
        * bridge_gap='0.2um'
        * finger_length='1.75um'
        * finger_width='0.23um'
        * prong_width='0.5um'
        * prong_length='1.75um'
        * t_pad_width='0.6um'
        * t_finger_width='0.23um'
        * t_finger_length='1um'
        * t_pad_length= '1um'
        * t_pad_extra='0.1um'
        * stem_width=`1um`
    """

    default_options = Dict(dist_extend='40um',
                           bridge_gap='0.2um',
                           finger_length='1.75um',
                           finger_width='0.230um',
                           prong_width='0.5um',
                           prong_length='1.75um',
                           t_pad_width='0.6um',
                           t_finger_width='0.230um',
                           t_finger_length='1um',
                           t_pad_length='3um',
                           t_pad_extra='0.1um',
                           stem_width='1um',
                           reverse=False)
    
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
        print(startPt)
        endPt = startPt + norm*p.dist_extend
        print(endPt)
        p.end_x = endPt[0]
        p.end_y = endPt[1]

        pad_T, pad_Fork, pin1, pin2, sim_JJ = JunctionSingleDolan.draw_junction(p)

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(pad1=pad_T, pad_Fork=pad_Fork),
                           layer=p.layer)

        self.add_qgeometry('junction',
                           {'design': sim_JJ},
                           layer=p.layer,
                           subtract=False,
                           width=p.squid_width+2*p.t_pad_extra)

        # Generates its own pins
        self.add_pin('t', pin1.coords[::-1], width=p.stem_width)
        self.add_pin('f', pin2.coords[::-1], width=p.stem_width)

class JunctionDolanPinStretch(QComponent):
    """Create a Dolan Bridge Josephson Junction

    Inherits QComponent class.

    The Dolan Bridge consists of a T-Section followed by a Fork Section.

    Dolan Bridge Josephson Junction Metal Geometry (no Ground Cutout):
        * bridge_gap - Gap between the T-Section and the Fork Section
        * finger_length - Length of the thin section on the end of the Fork Section
        * finger_width  - Width of the thin section on the end of the Fork Section
        * squid_width  - Overall width of the SQUID loop
        * stem_width  - Width of the leads connecting to the Josephson Junction
        * prong_width  - Length of the thicker section on the end of the Fork Section
        * prong_length  - Width of the thicker section on the end of the Fork Section
        * t_pad_width  - Thickness of the T-Section pad connecting to the leads
        * fork_pad_size  - Thickness of the Fork Section pad connecting to the leads
        * t_pad_width  - Thickness of the T-Section pad connecting to the leads
        * t_finger_width - Thicknes of the finger of the T-section pad connecting to the leads
        * t_pad_length - Length of the central pad of the T-section 
        * fork_pad_size  - Thickness of the Fork Section pad connecting to the leads
        * t_pad_extra - Extra length added to the sides of the T-Section to ensure proper overlap with fingers
        * reverse    - Default False. If True, the direction of the prongs is flipped by 180 degrees (that is, the prongs will
                       be near the first pin)


    The positioning can be done dynamically via:
        * pin_inputs=Dict(start_pin=Dict(component=f'...',pin='...')) - Specifying start position via a component pin
        * dist_extend - Distance upon to stretch away from the start pin.
    The resulting Josephson junction is right in the centre. This class ignores pos_x, pos_y and orientation...
        
    Pins:
        There are pins given on either side to help position the bandage. They are called 't' and 'f' for the T-pad and
        fork respectively. Pin width is stem_width.

    Sketch:
        Below is a sketch of the Josephson Junction Shadow Evaporation masking template (there is no ground cut-out)
        ::

                SW.............SW
                      | |
                      |S|                      SW  = squid_width
             _________| |_________             S   = stem_width 
            |      _________      |     FPS    FPS = fork_pad_size
            |     |         |PW.PW|     PL     PL  = prong_length
            |_   _|         |_   _|     PL     PW  = prong_width
              | |             | |       FW  
              | |             |F|       FW     F  = finger_width
              |_|             |_|       FW     FW = finger_length
                    
                  |<--TPL-->|           BG     BG = bridge_gap
                   _________            BG
         _________|    ^    |_________  BG                 
        |_________   TPW   __________|  TFW    TFW = t_finger_width  
                  |___ | ___|                  TPL = t_pad_length
        <-->          | |         <-->         TPW = t_pad_width
         TPE          |S|          TPE         TPE = t_pad_extra
                      | |

    .. image::
        Cap3Interdigital.png

    .. meta::
        Dolan Bridge Josephson Junction

    Default Options:
        * trace_width: '10um
        * dist_extend='40um'
        * bridge_gap='0.2um'
        * finger_length='1.75um'
        * finger_width='0.235um'
        * squid_width='5.4um'
        * stem_width='2um'
        * prong_width='0.5um'
        * prong_length='1.75um'
        * t_pad_width='0.3um'
        * t_finger_width='0.3um'
        * t_pad_length= '3um'
        * fork_pad_size='0.5um'
        * t_pad_extra='0.0um'
    """

    default_options = Dict(dist_extend='40um',
                           bridge_gap='0.2um',
                           finger_length='1.75um',
                           finger_width='0.235um',
                           squid_width='5.4um',
                           stem_width='2um',
                           prong_width='0.5um',
                           prong_length='1.75um',
                           t_pad_width='0.3um',
                           t_finger_width='0.3um',
                           t_pad_length='3um',
                           fork_pad_size='0.5um',
                           t_pad_extra='0.0um',
                           reverse=False)

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

        pad_T, pad_Fork, pin1, pin2, sim_JJ = JunctionDolan.draw_junction(p)

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(pad1=pad_T, pad_Fork=pad_Fork),
                           layer=p.layer)

        #subtracts out ground plane on the layer it's on
        # self.add_qgeometry('poly',
        #                    dict(padGap=padGap),
        #                    subtract=True,
        #                    layer=p.layer)

        self.add_qgeometry('junction',
                    {'design': sim_JJ},
                    layer=p.layer,
                    subtract=False,
                    width=p.squid_width+2*p.t_pad_extra)

        # Generates its own pins
        self.add_pin('t', pin1.coords[::-1], width=p.stem_width)
        self.add_pin('f', pin2.coords[::-1], width=p.stem_width)

# -*- coding: utf-8 -*-
# Author: Pradeep
# Extended and modified by: Stefan
# Creation Date: 2023
# Description: Class to draw Manhattan junctions.
class JunctionManhattan(QComponent):
    """Create a Manhattan junction

    Inherits QComponent class.

    The Manhattan junction consists of two electrodes and a junction area plus extension electrodes.
    The junction area is a rectangle with width jj_width and length jj_width2.
    The electrodes are rectangles with width width and length length.
    The extension electrodes are rectangles with width width and length extension.

    The positioning can be done via position (pos_x,pos_y) and angle given by 'orientation'.

    # more details on the class will be added shortly

    """

    default_options = Dict(pos_x='0mm',
                           pos_y='0mm',
                           jj_width = '0um',
                           width='200nm',
                           width2='-1',
                           length='10um',
                           length2='-1',
                           extension='2um',
                           extension2='-1',
                           pin1_offset='0um',
                           pin2_offset ='0um',
                           taper_width1 = '2um',
                           taper_width2 = '2um',
                           extension_taper = False,
                           start_point_shift = ('0um', '0um'),
                           four_legs = False)

    component_metadata = Dict(short_name='component')
    """Component metadata"""

    def make(self):
        """Convert self.options into QGeometry."""

        extension_taper = self.options.extension_taper
        p = self.parse_options()  # Parse the string options into numbers

        # if, for some reason, non-square junctions are needed, specify jj_width2
        if p.jj_width2:
            jj_width2 = p.jj_width2
        else:
            jj_width2 = p.jj_width

        if p.width2 == -1:
            width2 = p.width
        else:
            width2 = p.width2

        if p.length2 == -1:
            length2 = p.length
        else:
            length2 = p.length2

        if p.extension2 == -1:
            extension2 = p.extension
        else:
            extension2 = p.extension2

        # tapering for the extension electrodes: this will probably not be needed often
        if extension_taper:
            extension_taper_width1 = p.taper_width1
            extension_taper_width2 = p.taper_width2
        else:
            extension_taper_width1 = p.width  # noqa: F841 # abhishekchak52: unused variable extension_taper_width1
            extension_taper_width2 = p.width  # noqa: F841 # abhishekchak52: unused variable extension_taper_width2

        if p.target_comp:
            start_point = self.design.components[self.options.target_comp].pins[self.options.target_pin]
            start_point = start_point['middle'] + p.start_point_shift
            # print(start_point)
        else:
            start_point = (p.pos_x, p.pos_y)

        # draw the first electrode as a rectangle and add the tapering rectangles
        electrode1 = ShapelyEx.fuse_polygons_threshold([draw.rectangle(p.width,
                                    p.length + p.extension,
                                    0,
                                    -(p.length - p.extension)/2), # main electrode
                                    draw.rectangle(p.taper_width1, # taper 1
                                    p.length*0.9,
                                    0,
                                    (-p.length-p.length*0.1)/2),
                                    draw.rectangle(2*p.taper_width1, # taper 2
                                    p.length*0.6,
                                    0,
                                    (-p.length-p.length*0.4)/2), # extension taper 1
                                    # draw.rectangle(extension_taper_width1,
                                    # p.extension*0.9,
                                    # 0,
                                    # (p.extension+0.1*p.extension)/2),
                                    # draw.rectangle(2*extension_taper_width1, # extension taper 2
                                    # p.extension*0.6,
                                    # 0,
                                    # (p.extension+0.4*p.extension)/2)
                                    ])


        pin1_coords = draw.LineString([(0,0), (0,-p.length+p.pin1_offset)])

        # draw the first electrode as a rectangle and add the tapering rectangles
        electrode2 = ShapelyEx.fuse_polygons_threshold([draw.rectangle(width2,
                                    length2 + extension2,
                                    0,
                                    -(length2 - extension2)/2),
                                    # draw.rectangle(p.taper_width2,
                                    # length2*0.9,
                                    # 0,
                                    # (-p.length2-p.length2*0.1)/2),
                                    # draw.rectangle(2*p.taper_width2,
                                    # length2*0.6,
                                    # 0,
                                    # (-p.length2-p.length2*0.4)/2),
                                    # draw.rectangle(extension_taper_width2,
                                    # p.extension2*0.9,
                                    # 0,
                                    # (p.extension2+0.1*p.extension2)/2),
                                    # draw.rectangle(2*extension_taper_width2,
                                    # p.extension2*0.6,
                                    # 0,
                                    # (p.extension2+0.4*p.extension2)/2)
                                    ])


        # draw the second electrode as a rectangle
        electrode2 = draw.rotate_position(electrode2, 90, (0,0), (0,0))
        pin21_coords = draw.LineString([(0,0), (0,-p.length2+p.pin2_offset)])
        pin2_coords = draw.rotate_position(pin21_coords, 90, (0,0), (0,0))


        if p.four_legs:
            # draw junction area
            junction = draw.rectangle(p.jj_width,
                                    jj_width2,
                                    0,
                                    0)
        else: 
            junction = draw.rectangle(p.jj_width,
                        jj_width2,)
                        # (p.jj_width/2-p.width/2),
                        # -(jj_width2/2-p.width2/2))
        # extension pins
        extension_pin11 = draw.LineString([(0,0), (0,-p.extension)])
        extension_pin1 = draw.rotate_position(extension_pin11, 180, (0,0), (0,0))
        extension_pin21 = draw.LineString([(0,0), (0,-p.extension2)])
        extension_pin2 = draw.rotate_position(extension_pin21, -90, (0,0), (0,0))

        electrode1 = draw.rotate_position(electrode1, p.orientation, start_point)
        electrode2 = draw.rotate_position(electrode2, p.orientation, start_point)

        junction = draw.rotate_position(junction, p.orientation, start_point)
        # final_design = draw.rotate(final_design, p.orientation)
        pin1_coords = draw.rotate_position(pin1_coords, p.orientation, start_point)
        pin2_coords = draw.rotate_position(pin2_coords, p.orientation, start_point)

        extension_pin1 = draw.rotate_position(extension_pin1, p.orientation, start_point)
        extension_pin2 = draw.rotate_position(extension_pin2, p.orientation, start_point)



        geom = {'electrode1': electrode1,
                'electrode2': electrode2,
                'junction': junction,}
        self.add_qgeometry('poly', geom, layer=p.layer, subtract=False)
        self.add_pin('pin1', pin1_coords.coords, width=p.width, input_as_norm=True)
        self.add_pin('pin2', pin2_coords.coords, width=p.width, input_as_norm=True)
        self.add_pin('extension_pin1', extension_pin1.coords, width=p.width, input_as_norm=True)
        self.add_pin('extension_pin2', extension_pin2.coords, width=p.width, input_as_norm=True)

class JunctionDolanAsymmetric(QComponent):
    """Create a Dolan Bridge Josephson Junction

    Inherits QComponent class.

    The Dolan Bridge consists of a T-Section followed by a Fork Section.

    Dolan Bridge Josephson Junction Metal Geometry (no Ground Cutout):
        * bridge_gap - Gap between the T-Section and the Fork Section
        * finger_length - Length of the thin section on the end of the Fork Section
        * finger_width_top - Width of the thin section on the end of the top Fork Section
        * finger_width_bottom - Width of the thin section on the end of the bottom Fork Section
        * squid_width  - Overall width of the SQUID loop
        * stem_width  - Width of the leads connecting to the Josephson Junction
        * prong_width_top - Length of the thicker section on the end of the top Fork Section
        * prong_width_bottom - Length of the thicker section on the end of the botom Fork Section
        * prong_length  - Width of the thicker section on the end of the Fork Section
        * t_pad_width  - Thickness of the T-Section pad connecting to the leads
        * t_finger_width - Thicknes of the finger of the T-section pad connecting to the leads
        * t_pad_length - Length of the central pad of the T-section 
        * fork_pad_size  - Thickness of the Fork Section pad connecting to the leads
        * t_pad_extra - Extra length added to the sides of the T-Section to ensure proper overlap with fingers
        * reverse    - Default False. If True, the direction of the prongs is flipped by 180 degrees (that is, the prongs will
                       be near the first pin)

    As usual, the positioning can be done dynamically as a vector given by the supplied parameters: (pos_x,pos_y) to (end_x,end_y)
        
    Pins:
        There are pins given on either side to help position the bandage. They are called 't' and 'f' for the T-pad and
        fork respectively. Pin width is stem_width.

    Sketch:
        Below is a sketch of the Josephson Junction Shadow Evaporation masking template (there is no ground cut-out)
        ::

                SW.............SW
                       | |
                       |S|                     SW  = squid_width
             __________| |_________            S   = stem_width 
            |     __________      |     FPS    FPS = fork_pad_size
            |    |          |PW.PW|     PL     PL  = prong_length
            |_  _|          |_   _|     PL     PW  = prong_width
              ||              | |       FW  
              ||              |F|       FW     F  = finger_width
              ||              |_|       FW     FW = finger_length (rop/bottom)
                                        BG     BG = bridge_gap
                  |<--TPL-->|           BG     BG = bridge_gap
                   _________            BG
           _______|    ^    |_________  BG                 
          |________   TPW   __________| TFW    TFW = t_finger_width  
                  |___ | ___|                  TPL = t_pad_length
          <-->        | |        <-->          TPW = t_pad_width
           TPE        |S|        TPE           TP  = t_pad_extra
                      | |

    .. image::
        Cap3Interdigital.png

    .. meta::
        Dolan Bridge Josephson Junction

    Default Options:
        * trace_width: '10um
        * pos_x='0um',pos_y='0um'
        * end_x='50um',end_y='0um'
        * bridge_gap='0.2um'
        * finger_length='1.75um'
        * finger_width_top='0.235um'
        * finger_width_botom='0.135um'
        * squid_width='5.4um'
        * stem_width='2um'
        * prong_width_top='0.5um'
        * prong_width_bottom='0.45um'
        * prong_length='1.75um'
        * t_pad_width='0.3um'
        * t_finger_width='0.3um'
        * t_pad_length= '0.5um'
        * fork_pad_size='0.5um'
        * t_pad_extra='0.0um'
    """

    default_options = Dict(pos_x='0um',pos_y='0um',
                           end_x='50um',end_y='0um',
                           bridge_gap='0.2um',
                           finger_length='1.75um',
                           finger_width_top='0.235um',
                           finger_width_bottom='0.135um',
                           squid_width='5.4um',
                           stem_width='2um',
                           prong_width_top='0.5um',
                           prong_width_bottom='0.45um',
                           prong_length='1.75um',
                           t_pad_width='0.3um',
                           t_finger_width='0.3um',
                           t_pad_length='0.5um',
                           fork_pad_size='0.5um',
                           t_pad_extra='0.0um',
                           reverse=False)

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################
        #call the draw_junction from the JunctionDolanAsymmetric class in order to generate the T-Section, Fork Section, and pins 
        #for the Josephson Junction based on the provided parameters p.
        pad_T, pad_Fork, pin1, pin2, sim_JJ = JunctionDolanAsymmetric.draw_junction(p)

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(pad1=pad_T, pad_Fork=pad_Fork),
                           layer=p.layer)

        #subtracts out ground plane on the layer it's on
        # self.add_qgeometry('poly',
        #                    dict(padGap=padGap),
        #                    subtract=True,
        #                    layer=p.layer)
        self.add_qgeometry('junction',
                    {'design': sim_JJ},
                    layer=p.layer,
                    subtract=False,
                    width=p.squid_width+2*p.t_pad_extra)

        # Generates its own pins
        self.add_pin('t', pin1.coords[::-1], width=p.stem_width)
        self.add_pin('f', pin2.coords[::-1], width=p.stem_width)

    @staticmethod
    def draw_junction(p):
        #To calculate the total width of the structure such that Stem length can be calculated and adjusted from jj width and total component width
        struct_width = p.t_pad_width/2+ p.t_finger_width/2 + p.bridge_gap + p.finger_length + p.prong_length + p.fork_pad_size
        len_comp = np.sqrt((p.end_x-p.pos_x)**2+(p.end_y-p.pos_y)**2)
        len_stem = (len_comp - struct_width)/2
        t_finger_length=(p.squid_width-p.t_pad_length)/2

        #The T-Section and Stem
        pad_T = [
                (-p.t_pad_width*0.5, p.stem_width*0.5), #1
                (-len_stem-p.t_pad_width*0.5, p.stem_width*0.5), #2
                (-len_stem-p.t_pad_width*0.5, -p.stem_width*0.5), #3
                (-p.t_pad_width*0.5, -p.stem_width*0.5), #4
                (-p.t_pad_width*0.5,-p.t_pad_length*0.5), #4'
                (-p.t_finger_width*0.5,-p.t_pad_length*0.5), #5
                (-p.t_finger_width*0.5,-p.t_pad_length*0.5-t_finger_length-p.t_pad_extra), #6
                (p.t_finger_width*0.5,-p.t_pad_length*0.5-t_finger_length-p.t_pad_extra), #7
                (p.t_finger_width*0.5,-p.t_pad_length*0.5), #8
                (p.t_pad_width*0.5,-p.t_pad_length*0.5), #9
                #
                (p.t_pad_width*0.5,p.t_pad_length*0.5), #10
                (p.t_finger_width*0.5,p.t_pad_length*0.5), #11
                (p.t_finger_width*0.5,p.t_pad_length*0.5+t_finger_length+p.t_pad_extra), #12
                (-p.t_finger_width*0.5,p.t_pad_length*0.5+t_finger_length+p.t_pad_extra), #13
                (-p.t_finger_width*0.5,p.t_pad_length*0.5), #14
                (-p.t_pad_width*0.5,p.t_pad_length*0.5), #15
                ] 
    
        #The Fork Section and Stem
        pad_Fork = [
                   (p.finger_length+p.prong_length+p.fork_pad_size+len_stem, p.stem_width*0.5),
                   (p.finger_length+p.prong_length+p.fork_pad_size, p.stem_width*0.5),
                   (p.finger_length+p.prong_length+p.fork_pad_size, p.squid_width*0.5),
                   (p.finger_length+p.prong_length, p.squid_width*0.5),
                   (p.finger_length, p.squid_width*0.5),
                   (p.finger_length, (p.squid_width-p.prong_width_top+p.finger_width_top)*0.5),
                   (0, (p.squid_width-p.prong_width_top+p.finger_width_top)*0.5),
                   (0, (p.squid_width-p.prong_width_top-p.finger_width_top)*0.5),
                   (p.finger_length, (p.squid_width-p.prong_width_top-p.finger_width_top)*0.5),#play with this finger width and prong width
                   (p.finger_length, p.squid_width*0.5-p.prong_width_top), 
                   (p.finger_length+p.prong_length, p.squid_width*0.5-p.prong_width_top),
                   #
                   (p.finger_length+p.prong_length, -p.squid_width*0.5+p.prong_width_bottom),
                   (p.finger_length, -p.squid_width*0.5+p.prong_width_bottom),
                   (p.finger_length, (-p.squid_width+p.prong_width_bottom+p.finger_width_bottom)*0.5),
                   (0, (-p.squid_width+p.prong_width_bottom+p.finger_width_bottom)*0.5),
                   (0, (-p.squid_width+p.prong_width_bottom-p.finger_width_bottom)*0.5),
                   (p.finger_length, (-p.squid_width+p.prong_width_bottom-p.finger_width_bottom)*0.5),
                   (p.finger_length, -p.squid_width*0.5),
                   (p.finger_length+p.prong_length, -p.squid_width*0.5),
                   (p.finger_length+p.prong_length+p.fork_pad_size, -p.squid_width*0.5),
                   (p.finger_length+p.prong_length+p.fork_pad_size, -p.stem_width*0.5),
                   (p.finger_length+p.prong_length+p.fork_pad_size+len_stem, -p.stem_width*0.5)]
    
        pad_T = np.array(pad_T)
        pad_T[:,0] += len_stem + p.t_pad_width/2
        pad_Fork = np.array(pad_Fork)
        pad_Fork[:,0] += p.t_finger_width*0.5+ p.bridge_gap + len_stem + p.t_pad_width/2

        sim_JJ = shapely.LineString([
            np.mean(pad_T[1:3],axis=0),
            np.mean([pad_Fork[-1], pad_Fork[0]],axis=0)
        ])


        pin1 = shapely.LineString(pad_T[1:3])
        pin2 = shapely.LineString([pad_Fork[-1], pad_Fork[0]])
        pad_T = shapely.Polygon(pad_T)
        pad_Fork = shapely.Polygon(pad_Fork)

        polys = [pad_T, pad_Fork, pin1, pin2, sim_JJ]
        if p.reverse:
            polys = draw.translate(polys, -len_comp, 0)
            polys = draw.rotate(polys, np.pi, origin=(0, 0), use_radians=True)
        polys = draw.rotate(polys, np.arctan2(p.end_y-p.pos_y,p.end_x-p.pos_x), origin=(0, 0), use_radians=True)
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [pad_T, pad_Fork, pin1, pin2, sim_JJ] = polys
        return pad_T, pad_Fork, pin1, pin2, sim_JJ

class JunctionDolanAsymmetricPinStretch(QComponent):
    """Create a Dolan Bridge Josephson Junction

    Inherits QComponent class.

    The Dolan Bridge consists of a T-Section followed by a Fork Section.

    Dolan Bridge Josephson Junction Metal Geometry (no Ground Cutout):
        * bridge_gap - Gap between the T-Section and the Fork Section
        * finger_length - Length of the thin section on the end of the Fork Section
        * finger_width_top - Width of the thin section on the end of the top Fork Section
        * finger_width_bottom - Width of the thin section on the end of the bottom Fork Section
        * squid_width  - Overall width of the SQUID loop
        * stem_width  - Width of the leads connecting to the Josephson Junction
        * prong_width_top - Length of the thicker section on the end of the top Fork Section
        * prong_width_bottom - Length of the thicker section on the end of the botom Fork Section
        * prong_length  - Width of the thicker section on the end of the Fork Section
        * t_pad_width  - Thickness of the T-Section pad connecting to the leads
        * t_finger_width - Thicknes of the finger of the T-section pad connecting to the leads
        * t_pad_length - Length of the central pad of the T-section 
        * fork_pad_size  - Thickness of the Fork Section pad connecting to the leads
        * t_pad_extra - Extra length added to the sides of the T-Section to ensure proper overlap with fingers
        * reverse    - Default False. If True, the direction of the prongs is flipped by 180 degrees (that is, the prongs will
                       be near the first pin)


    The positioning can be done dynamically via:
        * pin_inputs=Dict(start_pin=Dict(component=f'...',pin='...')) - Specifying start position via a component pin
        * dist_extend - Distance upon to stretch away from the start pin.
    The resulting Josephson junction is right in the centre. This class ignores pos_x, pos_y and orientation...
        
    Pins:
        There are pins given on either side to help position the bandage. They are called 't' and 'f' for the T-pad and
        fork respectively. Pin width is stem_width.

    Sketch:
        Below is a sketch of the Josephson Junction Shadow Evaporation masking template (there is no ground cut-out)
        ::

                SW.............SW
                       | |
                       |S|                     SW  = squid_width
             __________| |_________            S   = stem_width 
            |      _________      |     FPS    FPS = fork_pad_size
            |    |          |PW.PW|     PL     PL  = prong_length
            |_  _|          |_   _|     PL     PW  = prong_width
              ||              | |       FW  
              ||              |F|       FW     F  = finger_width
              ||              |_|       FW     FW = finger_length
                                        BG     BG = bridge_gap
                  |<--TPL-->|           BG
                   _________            BG
           _______|    ^    |_________  BG                 
          |________   TPW   __________| TFW    TFW = t_finger_width  
                  |___ | ___|                  TPL = t_pad_length
          <-->        | |        <-->          TPW = t_pad_width
           TPE        |S|        TPE           TP  = t_pad_extra
                      | |
    .. image::
        Cap3Interdigital.png

    .. meta::
        Dolan Bridge Asymmetric Josephson Junction

    Default Options:
        * trace_width: '10um
        * pos_x='0um',pos_y='0um'
        * end_x='50um',end_y='0um'
        * bridge_gap='0.2um'
        * finger_length='1.75um'
        * finger_width_top='0.235um'
        * finger_width_botom='0.135um'
        * squid_width='5.4um'
        * stem_width='2um'
        * prong_width_top='0.5um'
        * prong_width_bottom='0.45um'
        * prong_length='1.75um'
        * t_pad_width='0.3um'
        * t_finger_width='0.3um'
        * t_pad_length= '0.5um'
        * fork_pad_size='0.5um'
        * t_pad_extra='0.0um'
    """

    default_options = Dict(pos_x='0um',pos_y='0um',
                           end_x='50um',end_y='0um',
                           bridge_gap='0.2um',
                           finger_length='1.75um',
                           finger_width_top='0.235um',
                           finger_width_bottom='0.135um',
                           squid_width='5.4um',
                           stem_width='2um',
                           prong_width_top='0.5um',
                           prong_width_bottom='0.45um',
                           prong_length='1.75um',
                           t_pad_width='0.3um',
                           t_finger_width='0.3um',
                           t_pad_length='0.5um',
                           fork_pad_size='0.5um',
                           t_pad_extra='0.0um',
                           reverse=False)


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

        pad_T, pad_Fork, pin1, pin2, sim_JJ = JunctionDolanAsymmetric.draw_junction(p)

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(pad1=pad_T, pad_Fork=pad_Fork),
                           layer=p.layer)

        #subtracts out ground plane on the layer it's on
        # self.add_qgeometry('poly',
        #                    dict(padGap=padGap),
        #                    subtract=True,
        #                    layer=p.layer)
        self.add_qgeometry('junction',
                    {'design': sim_JJ},
                    layer=p.layer,
                    subtract=False,
                    width=p.squid_width+2*p.t_pad_extra)

        # Generates its own pins
        self.add_pin('t', pin1.coords[::-1], width=p.stem_width)
        self.add_pin('f', pin2.coords[::-1], width=p.stem_width)

def get_square_JJ_width(J_C_uA_um2, target_EJ_GHz=None, target_LJ_nH=None, rounding=True):
    """
    Function to calculate the dimensions for a Josephson junction fabricated with a
    certain critical current density (J_C), where I_C = J_C * A. The critical current 
    density must be supplied in units of micro-Amperes per square micro-metre (uA/um^2).

    A target Josephson inductance (L_J, unbits of nH) or Josephson energy (E_J, units of GHz)
    must be given. Accepts a list of values, or a single value.

    Returns width = height in um of the required JJ area. 
    """
    #assert (target_EJ_GHz is not None) or (target_LJ_nH is not None), "Must supply target EJ or LJ."
    assert not ((target_LJ_nH is not None) and (target_EJ_GHz is not None)), "Only supply either EJ or LJ, not both."
    assert isinstance(J_C_uA_um2, (float, int))
    phi_0 = 2.067833848 * 1e-15 # Wb
    h = 6.62607015 * 1e-34 # J s
    J_C_A_m2 = J_C_uA_um2 * 1e6
    # calculate areas for supplied LJ values
    if target_LJ_nH is not None:
        if not isinstance(target_LJ_nH, np.ndarray):
            target_LJ_nH = np.atleast_1d(np.array(target_LJ_nH, dtype=float))
        target_LJ_H = target_LJ_nH * 1e-9 # convert to H
        A_m2 = np.array([phi_0/(2 * np.pi * J_C_A_m2 * i) for i in target_LJ_H])
    # calculate areas for supplied EJ values
    elif target_EJ_GHz is not None:
        if not isinstance(target_EJ_GHz, np.ndarray):
            target_EJ_GHz = np.atleast_1d(np.array(target_EJ_GHz, dtype=float))
        target_EJ_Hz = target_EJ_GHz * 1e9
        A_m2 = np.array([(i * h * 2 * np.pi)/(phi_0 * J_C_A_m2) for i in target_EJ_Hz])
    else:
        raise ValueError("No areas were calculated since no target values were supplied.")
    width_JJ_nm = np.sqrt(A_m2) * 1e9
    if rounding:
        width_JJ_nm = np.array([(round(i / 1.0) * 1) for i in width_JJ_nm])
    width_JJ_um = width_JJ_nm * 1e-3
    return width_JJ_um.item() if width_JJ_um.size == 1 else width_JJ_um

