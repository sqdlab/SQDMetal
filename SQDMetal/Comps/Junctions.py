# -*- coding: utf-8 -*-
# Author: Prasanna Pakkiam
# Creation Date: 24/05/2023
# Description: Collection of classes to draw Josephson Junctions.

from qiskit_metal import draw
from qiskit_metal.toolbox_python.attr_dict import Dict
from qiskit_metal.qlibrary.core import QComponent
import numpy as np
import shapely

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
        * t_pad_size  - Thickness of the T-Section pad connecting to the leads
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
                       |S|                   SW  = squid_width
                 ______| |______             S   = stem_width 
                |      ___      |     FPS    FPS = fork_pad_size
                |     |   |PW.PW|     PL     PL  = prong_length
                |_   _|   |_   _|     PL     PW  = prong_width
                  | |       | |       FW  
                  | |       |F|       FW     F  = finger_width
                  |_|       |_|       FW     FW = finger_width
                                      BG     BG = bridge_gap
             _______________________  BG
            |__________   __________| TPS    TPS = t_pad_size
            <-->       | |       <-->        TPE = t_pad_extra
             TPE       |S|        TPE
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
        * t_pad_size='0.3um'
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
                           t_pad_size='0.3um',
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
        struct_width = p.t_pad_size + p.bridge_gap + p.finger_length + p.prong_length + p.fork_pad_size
        len_comp = np.sqrt((p.end_x-p.pos_x)**2+(p.end_y-p.pos_y)**2)
        len_stem = (len_comp - struct_width)/2

        #The T-Section and Stem
        pad_T = [
                (0, p.stem_width*0.5),
                (-len_stem, p.stem_width*0.5),
                (-len_stem, -p.stem_width*0.5),
                (0, -p.stem_width*0.5),
                (0, -p.squid_width*0.5-p.t_pad_extra),
                (p.t_pad_size, -p.squid_width*0.5-p.t_pad_extra),
                (p.t_pad_size, p.squid_width*0.5+p.t_pad_extra),
                (0, p.squid_width*0.5+p.t_pad_extra)]
    
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
        pad_T[:,0] += len_stem
        pad_Fork = np.array(pad_Fork)
        pad_Fork[:,0] += p.t_pad_size + p.bridge_gap + len_stem

        pin1 = shapely.LineString(pad_T[1:3])
        pin2 = shapely.LineString([pad_Fork[-1], pad_Fork[0]])

        sim_JJ = shapely.LineString([
            np.mean(pad_T[1:3],axis=0),
            np.mean([pad_Fork[-1], pad_Fork[0]],axis=0)
        ])

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
        * t_pad_size  - Thickness of the T-Section pad connecting to the leads
        * fork_pad_size  - Thickness of the Fork Section pad connecting to the leads
        * t_pad_extra - Extra length added to the sides of the T-Section to ensure proper overlap with fingers
        * reverse    - Default False. If True, the direction of the prongs is flipped by 180 degrees (that is, the prongs will
                       be near the target pin)


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
                       |S|                   SW  = squid_width
                 ______| |______             S   = stem_width 
                |      ___      |     FPS    FPS = fork_pad_size
                |     |   |PW.PW|     PL     PL  = prong_length
                |_   _|   |_   _|     PL     PW  = prong_width
                  | |       | |       FW  
                  | |       |F|       FW     F  = finger_width
                  |_|       |_|       FW     FW = finger_width
                                      BG     BG = bridge_gap
             _______________________  BG
            |__________   __________| TPS    TPS = t_pad_size
            <-->       | |       <-->        TPE = t_pad_extra
             TPE       |S|        TPE        x   = Location where target poin is attached.
                       |x|

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
        * t_pad_size='0.3um'
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
                           t_pad_size='0.3um',
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
