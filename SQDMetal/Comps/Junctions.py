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

    As usual, the positioning can be done dynamically as a vector given by the supplied parameters: (pos_x,pos_y) to (end_x,end_y)
        
    Pins:
        There are no pins given that overlap and precise positioning is usually the concern...

    Sketch:
        Below is a sketch of the Josephson Junction Shadow Evaporation masking template (there is no ground cut-out)
        ::

            SW.............SW
                   | |
                   |S|              SW  = squid_width
             ______| |______        S   = stem_width 
            |      ___      | FPS   FPS = fork_pad_size
            |     |   |PW.PW| PL    PL  = prong_length
            |_   _|   |_   _| PL    PW  = prong_width
              | |       | |   FW 
              | |       |F|   FW    F  = finger_width
              |_|       |_|   FW    FW = finger_width
                              BG    BG = bridge_gap
             _______________  BG
            |______   ______| TPS   TPS = t_pad_size
                   | |
                   |S|        
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
                           fork_pad_size='0.5um')

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################

        struct_width = p.t_pad_size + p.bridge_gap + p.finger_length + p.prong_length + p.fork_pad_size
        len_comp = np.sqrt((p.end_x-p.pos_x)**2+(p.end_y-p.pos_y)**2)
        len_stem = (len_comp - struct_width)/2

        #The T-Section and Stem
        pad_T = [
                (0, p.stem_width*0.5),
                (-len_stem, p.stem_width*0.5),
                (-len_stem, -p.stem_width*0.5),
                (0, -p.stem_width*0.5),
                (0, -p.squid_width*0.5),
                (p.t_pad_size, -p.squid_width*0.5),
                (p.t_pad_size, p.squid_width*0.5),
                (0, p.squid_width*0.5)]
    
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

        pad_T = shapely.Polygon(pad_T)
        pad_Fork = shapely.Polygon(pad_Fork)

        polys = [pad_T, pad_Fork]
        polys = draw.rotate(polys, np.arctan2(p.end_y-p.pos_y,p.end_x-p.pos_x), origin=(0, 0), use_radians=True)
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [pad_T, pad_Fork] = polys


        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(pad1=pad_T, pad_Fork=pad_Fork),
                           layer=p.layer)

        #subtracts out ground plane on the layer it's on
        # self.add_qgeometry('poly',
        #                    dict(padGap=padGap),
        #                    subtract=True,
        #                    layer=p.layer)

        # Generates its own pins
        # self.add_pin('a', pin1.coords[::-1], width=p.cpw_width)
        # self.add_pin('b', pin2.coords[::-1], width=p.cpw_width)
