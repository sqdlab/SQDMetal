from qiskit_metal.qlibrary.core import QComponent, QRoute
from qiskit_metal import draw
import numpy as np
import shapely
from qiskit_metal.toolbox_python.attr_dict import Dict
from SQDMetal.Utilities.QUtilities import QUtilities

class HallBar(QComponent):
    """
    A Hall Bar component consisting of 8 probe pads for resistance measurements.

    Inherits QComponent class.

    Probe Pad geometry:
        * pad_width     - Width (x) of probe pads
        * pad_height    - Height (y) of probe pads
        * pocket_width  - Width (x) of ground plane cutout
        * pocket_height - Height (x) of ground plane cutout

    Hall Bar geometry:
        * trace_width   - Width of the main conductor
        * T_height      - Height of the T-bars (along y)
        * num_sq_short  - Number of squares along the conductor between the two inner probe pads
        * num_sq_long   - Number of squares along the conductor between the two outer probe pads

    Positioning:
        * pos_x         - X position
        * pos_y         - Y position
        * layer         - layer number

    Pins:
        There are no pins.

    Sketch:
        Below is a sketch of the Hall bar

        Total width determined by num_sq_long (default 100).

        <-----NSL*TW----->          NSS : num_sq_short
         _   _      _   _           NSL : num_sq_short
        |_| |_|    |_| |_|     TH   TH  : T_height
        |     |    |     |     TH   TW  : trace_width       
        |_____|____|_____|     TH   PW  : pad_width   
        |     |    |     |     TH   PH  : pad_height  
        |_   _|    |_   _|     TH
        |_| |_|    |_| |_| PH  TH   y
               <-->     PW          |
               NSS                  --- x
        ::
            
    """
    default_options = Dict(pos_x="0mm",
                           pos_y="0mm",
                           trace_width='10um',
                           num_sq_short=10,
                           num_sq_long=100,
                           pad_width="150um",
                           pad_height="100um",
                           pocket_width="1.4mm",
                           pocket_height="1.0mm",
                           T_height="500um",
                           fillet_radius="20um",
                           layer=0,
                           pocket=True
                        )
    
    def __init__(self, design,
                    name: str = None,
                    options: Dict = None,
                    type: str = "CPW",
                    **kwargs):
        super().__init__(design, name, options, **kwargs)
    
    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        trace_short = p.num_sq_short * p.trace_width
        trace_long = p.num_sq_long * p.trace_width

        # draw
        main_trace = draw.rectangle(trace_long, p.trace_width)
        pocket = draw.rectangle(p.pocket_width, p.pocket_height)
        far_T_L = draw.rectangle(p.trace_width, p.T_height,  xoff=-0.5*trace_long)
        far_T_R = draw.rectangle(p.trace_width, p.T_height,  xoff=0.5*trace_long)
        near_T_L = draw.rectangle(p.trace_width, p.T_height, xoff=-0.5*trace_short)
        near_T_R = draw.rectangle(p.trace_width, p.T_height, xoff=0.5*trace_short)
        pad_U_L = draw.rectangle(p.pad_width, p.pad_height,  xoff=-0.5*trace_long+0.5*p.pad_width-0.5*p.trace_width,  yoff=0.5*p.T_height)
        pad_U_ML = draw.rectangle(p.pad_width, p.pad_height, xoff=-0.5*trace_short-0.5*p.pad_width+0.5*p.trace_width, yoff=0.5*p.T_height)
        pad_U_MR = draw.rectangle(p.pad_width, p.pad_height, xoff=0.5*trace_short+0.5*p.pad_width-0.5*p.trace_width,  yoff=0.5*p.T_height)
        pad_U_R = draw.rectangle(p.pad_width, p.pad_height,  xoff=0.5*trace_long-0.5*p.pad_width+0.5*p.trace_width,  yoff=0.5*p.T_height)
        pad_L_L = draw.rectangle(p.pad_width, p.pad_height,  xoff=-0.5*trace_long+0.5*p.pad_width-0.5*p.trace_width,  yoff=-0.5*p.T_height)
        pad_L_ML = draw.rectangle(p.pad_width, p.pad_height, xoff=-0.5*trace_short-0.5*p.pad_width+0.5*p.trace_width, yoff=-0.5*p.T_height)
        pad_L_MR = draw.rectangle(p.pad_width, p.pad_height, xoff=0.5*trace_short+0.5*p.pad_width-0.5*p.trace_width,  yoff=-0.5*p.T_height)
        pad_L_R = draw.rectangle(p.pad_width, p.pad_height,  xoff=0.5*trace_long-0.5*p.pad_width+0.5*p.trace_width,  yoff=-0.5*p.T_height)

        # merge 
        T_ends = shapely.ops.unary_union([far_T_L, far_T_R])
        T_mids = shapely.ops.unary_union([near_T_L, near_T_R])
        pads_L = shapely.ops.unary_union([pad_L_L, pad_L_ML, pad_L_MR, pad_L_R])
        pads_U = shapely.ops.unary_union([pad_U_L, pad_U_ML, pad_U_MR, pad_U_R])
        # merge all pads and T's
        pads_and_Ts = shapely.ops.unary_union([T_ends, T_mids, pads_L, pads_U])
        pads_and_Ts = pads_and_Ts.buffer(p.fillet_radius).buffer(-p.fillet_radius)
        pads_and_Ts = pads_and_Ts.buffer(-0.48*p.trace_width).buffer(0.48*p.trace_width)

        # rotate and translate
        polys = [main_trace, pocket, pads_and_Ts]
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [main_trace, pocket, pads_and_Ts] = polys

        # subtract pocket
        if p.pocket:
            self.add_qgeometry('poly', 
                            dict(pocket=pocket), 
                            subtract=True,
                            layer=p.layer)

        # add metals
        self.add_qgeometry('poly', 
                           dict(trace=main_trace, 
                                pads=pads_and_Ts
                                ), 
                           subtract=False,
                           layer=p.layer)
        
class FourProbeJJ(QComponent):

    # default options
    default_options = Dict(
        with_labels=False,
        pos_x="0um",
        pos_y="0um",
        pocket_size="800um",
        pad_size="150um",
        pad_dist="200um",
        pin_gap="40um",
        wire_width="8um",
        component_prefix="",
        sample_label="",
        orientation=0,
        cap_style=3,
        pocket=True
    )

    # metadata
    component_metadata = Dict(short_name='fourProbePins',
                              _qgeometry_table_path='True',
                              _qgeometry_table_poly='True',
                              _qgeometry_table_junction='True')

    # def make(self):
    def __init__(self, design,
                    name: str = None,
                    options: Dict = None,
                    **kwargs):
        
        super().__init__(design, name, options, **kwargs)

    def make(self):
        # parse options
        p = self.p
        x_off = p.pos_x
        y_off = p.pos_y
        pad_w = p.pad_size
        pad_h = pad_w
        pad_off = p.pad_dist
        pin_gap = p.pin_gap
        w = p.wire_width
        p_w = p.pocket_size
        c = p.cap_style
        pin_len = 1e-6

        # pocket
        pocket = draw.rectangle(p_w, p_w)

        # probe pads
        pad1 = draw.rectangle(pad_w, pad_h, 0, pad_off)
        pad2 = draw.rectangle(pad_w, pad_h, pad_off, 0)
        pad3 = draw.rectangle(pad_w, pad_h, 0, -pad_off)
        pad4 = draw.rectangle(pad_w, pad_h, -pad_off, 0)

        # probe wires
        wire1 = draw.LineString([[0, pad_off], [0, pin_gap/2]]).buffer(w/2, cap_style=c)
        wire2 = draw.LineString([[pad_off, 0], [pin_gap, 0], [pin_gap/2, pin_gap/2], [0, pin_gap/2]]).buffer(w/2, cap_style=c)
        wire3 = draw.LineString([[0, -pad_off], [0, -pin_gap/2]]).buffer(w/2, cap_style=c)
        wire4 = draw.LineString([[-pad_off, 0], [-pin_gap, 0], [-pin_gap/2, -pin_gap/2], [0, -pin_gap/2]]).buffer(w/2, cap_style=c)

        # join pads and wires
        pad1 = draw.union([pad1, wire1])
        pad2 = draw.union([pad2, wire2])
        pad3 = draw.union([pad3, wire3])
        pad4 = draw.union([pad4, wire4])

        # define pins (pointing from wire ends to origin)
        pin_top = draw.LineString([[0, pin_gap/2 + pin_len], [0, pin_gap/2]])
        pin_bottom = draw.LineString([[0, -pin_gap/2 - pin_len], [0, -pin_gap/2]])

        # bulk translate and rotate
        polys = [pad1, pad2, pad3, pad4, pin_top, pin_bottom, pocket]
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, x_off, y_off)
        [pad1, pad2, pad3, pad4, pin_top, pin_bottom, pocket] = polys

        # add QGeometries - pocket
        if p.pocket:
            self.add_qgeometry('poly',
                            dict(pocket=pocket),
                            subtract=True)

        # add QGeometries - pads and wires
        self.add_qgeometry('poly', 
                           dict(
                               pad_top=pad1, 
                               pad_right=pad2,
                               pad_bottom=pad3,
                               pad_left=pad4
                            ))
        
        # add pins to connect junction component
        self.add_pin("pin_top", 
                     points=list(pin_top.coords), 
                     width=w, 
                     input_as_norm=True)
        self.add_pin("pin_bottom", 
                     points=list(pin_bottom.coords), 
                     width=w, 
                     input_as_norm=True)