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
        |_________   TPW   __________| TFW    TFW = t_finger_width  
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




