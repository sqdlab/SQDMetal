from qiskit_metal.qlibrary.core import QComponent, QRoute
from qiskit_metal import draw
import numpy as np
import shapely
from qiskit_metal.toolbox_python.attr_dict import Dict
from SQDMetal.Utilities.QUtilities import QUtilities

class BacklessLaunchpad(QComponent):
    r"""Launch pad to feed/read signals to/from the chip. Almost identical to the one built into quantum metal except
    the pocket at the opposite end of the pin is intentionally small. Intended for when the pad is at the edge of the coupon
    and a majority of the pocket will be diced off anyways

    .. meta::
        :description: Launchpad Wirebond

    Creates a 50 ohm launch pad with a ground pocket cutout.
    Limited but expandable parameters to control the launchpad polygons.
    The (0,0) point is the center of the necking of the launch tip.
    The pin attaches directly to the built in lead length at its midpoint

    Pocket and pad:
        Pocket and launch pad geometries are currently fixed.
        (0,0) point is the midpoint of the necking of the launch tip.
        Pocket is a negative shape that is cut out of the ground plane

    Values (unless noted) are strings with units included, (e.g., '30um')

    Sketch:
        Below is a sketch of the launch
        ::

            -----------
            |          \
            |      ---------\\
            |      |    0    |    (0,0) pin at midpoint of necking, before the lead
            |      ---------//
            |          /
            -----------

            y
            ^
            |
            |------> x

    .. image::
        LaunchpadWirebond.png

    Default Options:
        * trace_width: 'cpw_width' -- Width of the transmission line attached to the launch pad
        * trace_gap: 'cpw_gap' -- Gap of the transmission line
        * lead_length: '25um' -- Length of the transmission line attached to the launch pad
        * pad_width: '80um' -- Width of the launch pad
        * pad_height: '80um' -- Height of the launch pad
        * pad_gap: '58um' -- Gap of the launch pad
        * taper_height: '122um' -- Height of the taper from the launch pad to the transmission line
        * back_gap: '5um' -- The gap at the back of the launchpad
    """

    default_options = Dict(
        trace_width="cpw_width",
        trace_gap="cpw_gap",
        lead_length="25um",
        pad_width="80um",
        pad_height="80um",
        pad_gap="58um",
        taper_height="122um",
        back_gap="5um"
    )
    """Default options"""

    TOOLTIP = """Launch pad to feed/read signals to/from the chip."""

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""

        p = self.p

        pad_width = p.pad_width
        pad_height = p.pad_height
        pad_gap = p.pad_gap
        trace_width = p.trace_width
        trace_width_half = trace_width / 2.0
        pad_width_half = pad_width / 2.0
        lead_length = p.lead_length
        taper_height = p.taper_height
        trace_gap = p.trace_gap
        back_gap=p.back_gap

        pad_gap = p.pad_gap
        #########################################################

        # Geometry of main launch structure
        # The shape is a polygon and we prepare this point as orientation is 0 degree
        launch_pad = draw.Polygon(
            [
                (0, trace_width_half),
                (-taper_height, pad_width_half),
                (-(pad_height + taper_height), pad_width_half),
                (-(pad_height + taper_height), -pad_width_half),
                (-taper_height, -pad_width_half),
                (0, -trace_width_half),
                (lead_length, -trace_width_half),
                (lead_length, trace_width_half),
                (0, trace_width_half),
            ]
        )

        # Geometry pocket (gap)
        # Same way applied for pocket
        pocket = draw.Polygon(
            [
                (0, trace_width_half + trace_gap),
                (-taper_height, pad_width_half + pad_gap),
                (-(pad_height + taper_height + back_gap), pad_width_half + pad_gap),
                (-(pad_height + taper_height + back_gap), -(pad_width_half + pad_gap)),
                (-taper_height, -(pad_width_half + pad_gap)),
                (0, -(trace_width_half + trace_gap)),
                (lead_length, -(trace_width_half + trace_gap)),
                (lead_length, trace_width_half + trace_gap),
                (0, trace_width_half + trace_gap),
            ]
        )

        # These variables are used to graphically locate the pin locations
        main_pin_line = draw.LineString(
            [(lead_length, trace_width_half), (lead_length, -trace_width_half)]
        )

        # Create polygon object list
        polys1 = [main_pin_line, launch_pad, pocket]

        # Rotates and translates all the objects as requested. Uses package functions in
        # 'draw_utility' for easy rotation/translation
        polys1 = draw.rotate(polys1, p.orientation, origin=(0, 0))
        polys1 = draw.translate(polys1, xoff=p.pos_x, yoff=p.pos_y)
        [main_pin_line, launch_pad, pocket] = polys1

        # Adds the object to the qgeometry table
        self.add_qgeometry("poly", dict(launch_pad=launch_pad), layer=p.layer)

        # Subtracts out ground plane on the layer its on
        self.add_qgeometry("poly", dict(pocket=pocket), subtract=True, layer=p.layer)

        # Generates the pins
        self.add_pin("tie", main_pin_line.coords, trace_width)
        
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