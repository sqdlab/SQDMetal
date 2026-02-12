
### This file is part of the SQDMetal project: Test to create a new Capacitor class
### Author: Léon Thinus

import numpy as np
from qiskit_metal import draw, Dict
from qiskit_metal.qlibrary.core import BaseQubit
from shapely.geometry.base import CAP_STYLE
from shapely.geometry import Polygon
from SQDMetal.Utilities.QUtilities import QUtilities 
from shapely.geometry import LineString
from qiskit_metal.toolbox_python.attr_dict import Dict
from qiskit_metal.qlibrary.core import QComponent
import shapely
from SQDMetal.Utilities.CpwParams import CpwParams
from SQDMetal.Comps.Resonators import ResonatorMeander
from SQDMetal.Utilities.QiskitShapelyRenderer import QiskitShapelyRenderer
from SQDMetal.Utilities.ShapelyEx import ShapelyEx

class TransmonTaperedCapacitor(QComponent):
    """Test of a capacitor design for the Tapered Pocket Transmon. With smooth edges

    Inherits `QComponent` class
     The qubit-resonator coupling has been modified to minimise participation ratios
        Insets (optional) on the pads for additional couplings

    Description:
        Create a standard capacitor design  that looks like a synapse
   
    Options:
        Convention: Values (unless noted) are strings with units included,
        (e.g., '30um')

    Creates two pads upon which one pad is a fork that wraps around the target lead. The ground plane bisects the two pads.

    Inherits QComponent class.

    Capacitor Metal Geometry consisting of a pad wrapped around by a fork (2 prongs):
        * trace_width  - Width of the first pad
        * trace_length - Length of the first pad
        * prong_trace_gap - Gap between the prong and the first pad
        * pad_trace_gap   - Gap between the two pads
        * gnd_prong_trace - Thickness of the ground plane that bisects the two pads
        * prong_width   - Width of the two prongs
        * prong_length  - Length of the two prongs
        * pad_thickness - Length of the Capacitor fork pad
        * gnd_prong_trace_dist - If less than zero, the ground bisecting the prong and the trace is in the centre. Otherwise,
                                 this value sets the gap from the prongs.
        * gnd_pad_trace_dist   - If less than zero, the ground bisecting the pad and the trace is in the centre. Otherwise,
                                 this value sets the gap from the pads.
        * swap_direction - If True, the fork is the first pad and the trace is the second pad instead.
        
    Spacing around the structure (i.e. cuts into the ground plane) can be controlled via:
        * gap_side  - Spacing to ground plane on the outer sides of the fork
        * gap_front - Spacing to ground plane from the prongs in the direction of the first pad
        * gap_back - Spacing to ground plane behind the forked pad

    The positioning can be done dynamically via:
        * pin_inputs=Dict(start_pin=Dict(component=f'...',pin='...')) - Specifying the target component pin

    Pins:
        There is one pin 'a' to link to the forked pad. It has a width equal to trace_width.

    Sketch:
        Below is a sketch of the capacitor
        ::

         @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@     @ = Ground Plane
         @                                 B  @     t = trace_width
         @    ____________aaaa____________ B  @     l = trace_length
         @   |                            |T  @     S = prong_trace_gap             x = Target pin
         @   |   ______________________   |T  @     P = pad_trace_gap               a = Connecting pin generated on this part
         @   |  |          P      H    |  |L  @     g = gnd_prong_trace             G = gnd_prong_trace_dist
         @   |  |   @@@@@@@P@@@@@@@@f  |  |L  @     f = gnd_pad_trace               H = gnd_pad_trace_dist
         @sss|  |GGG@@@@@@@P@@@@@@@@f  |WW|L  @     W = prong_width
         @   |  |   @@    _P__    @@   |  |L  @     L = prong_length
         @   |__|SSSSSSSS|    |l  @@   |__|L  @     T = pad_thickness
         @          @@   |tttt|l  @@       f  @     s = gap_side
         @@@@@@@@@@@@@   |    |l  @@@@@@@@@@@@@     f = gap_front
                           xx     gg                B = gap_back

    .. image::
        Cap3Interdigital.png

    .. meta::
        Cap 3 Interdigital

    Default Options:
        * trace_width='20um'
        * trace_length='30um'
        * prong_trace_gap='5um'
        * pad_trace_gap='5um'
        * gnd_prong_trace='2um'
        * gnd_pad_trace='2um'
        * prong_width='4um'
        * prong_length='10um'
        * pad_thickness='10um'
        * gnd_prong_trace_dist='2um'
        * gnd_pad_trace_dist='2um'
        * swap_direction=False
        * gap_side='5um'
        * gap_front='5um'
        * gap_back='5um'
        * orientation = 0

    Edge options:
        
        ### fillet options for the main body of the capacitor
        * fillet_radius = '5um' - Radius of edge fillets. If 0, no fillets are created.
        * fillet_resolution = 16 - Resolution of the fillets in number of points. Higher number means smoother fillets but longer simulation time.

        ### Fillet options for the prong part of the capacitor

        ### Fillet options for the pocket, ground plane cuts: 
        * gap_fillet_radius = '5um' - Fillet radius for the cuts in the ground plane around the capacitor. If 0, no fillets are created.


    """

    #  Define structure functions

    default_options = Dict(trace_width='20um',
                           trace_length='30um',
                           prong_trace_gap='5um',
                           pad_trace_gap='10um',
                           gnd_prong_trace='2um',
                           gnd_pad_trace='2um',
                           prong_width='4um',
                           prong_length='10um',
                           pad_thickness='10um',
                           gnd_prong_trace_dist='2um',
                           gnd_pad_trace_dist='2um',
                           swap_direction=False,
                           gap_side='5um',
                           gap_front='5um',
                           gap_back='5um',
                           fillet_radius = '5um',
                           fillet_resolution = 16,
                           gap_fillet_radius = '5um',
                           orientation = 0)
    """Default drawing options"""

    TOOLTIP = "Create a three finger planar capacitor with a ground pocket cuttout."

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################

        start_point = self.design.components[self.options.pin_inputs.start_pin.component].pins[self.options.pin_inputs.start_pin.pin]
        startPt = start_point['middle']
        startPtNorm = start_point['normal']

        wid = start_point['width']  # noqa: F841 # abhishekchak52: unused variable wid
         #define the polygons for the basic capacitor structure (before filleting and rotation) 
        pad1 = [                                                                                 # 1st pad - the trace
               (0, p.trace_width*0.5),
               (0, -p.trace_width*0.5),
               (p.trace_length, -p.trace_width*0.5),
               (p.trace_length, p.trace_width*0.5)]
        
        if p.gnd_prong_trace_dist > 0:
            gap_to_inner_gnd_side = p.prong_trace_gap - p.gnd_prong_trace - p.gnd_prong_trace_dist
        else:
            gap_to_inner_gnd_side = (p.prong_trace_gap - p.gnd_prong_trace)*0.5
        if p.gnd_pad_trace_dist > 0:
            gap_to_inner_gnd_front = p.pad_trace_gap - p.gnd_pad_trace - p.gnd_pad_trace_dist
        else:
            gap_to_inner_gnd_front = (p.pad_trace_gap - p.gnd_pad_trace)*0.5
        gap1 = [                                                                                 # gap around 1st pad to create ground plane cutout
               (0, p.trace_width*0.5+gap_to_inner_gnd_side),
               (0, -p.trace_width*0.5-gap_to_inner_gnd_side),
               (p.trace_length+gap_to_inner_gnd_front, -p.trace_width*0.5-gap_to_inner_gnd_side),
               (p.trace_length+gap_to_inner_gnd_front, p.trace_width*0.5+gap_to_inner_gnd_side)]

        pad2 = [                                                                                 # 2nd pad - the fork that wraps around the target pin with principal pad
               (p.trace_length+p.pad_trace_gap+p.pad_thickness, p.trace_width*0.5+p.prong_trace_gap+p.prong_width),
               (p.trace_length+p.pad_trace_gap-p.prong_length, p.trace_width*0.5+p.prong_trace_gap+p.prong_width),
               (p.trace_length+p.pad_trace_gap-p.prong_length, p.trace_width*0.5+p.prong_trace_gap),
               (p.trace_length+p.pad_trace_gap, p.trace_width*0.5+p.prong_trace_gap),
               (p.trace_length+p.pad_trace_gap, -p.trace_width*0.5-p.prong_trace_gap),
               (p.trace_length+p.pad_trace_gap-p.prong_length, -p.trace_width*0.5-p.prong_trace_gap),
               (p.trace_length+p.pad_trace_gap-p.prong_length, -p.trace_width*0.5-p.prong_trace_gap-p.prong_width),
               (p.trace_length+p.pad_trace_gap+p.pad_thickness, -p.trace_width*0.5-p.prong_trace_gap-p.prong_width)]

        if p.gnd_prong_trace_dist > 0:
            gap_to_inner_gnd_side = p.gnd_prong_trace_dist
        if p.gnd_pad_trace_dist > 0:
            gap_to_inner_gnd_front = p.gnd_pad_trace_dist
        gap2 = [                                                                                # gap around 2nd pad to create ground plane cutout
               (p.trace_length+p.pad_trace_gap+p.pad_thickness+p.gap_back, p.trace_width*0.5+p.prong_trace_gap+p.prong_width+p.gap_side),
               (p.trace_length+p.pad_trace_gap-p.prong_length-p.gap_front, p.trace_width*0.5+p.prong_trace_gap+p.prong_width+p.gap_side),
               (p.trace_length+p.pad_trace_gap-p.prong_length-p.gap_front, p.trace_width*0.5+p.prong_trace_gap-gap_to_inner_gnd_side),
               (p.trace_length+p.pad_trace_gap-gap_to_inner_gnd_front, p.trace_width*0.5+p.prong_trace_gap-gap_to_inner_gnd_side),
               (p.trace_length+p.pad_trace_gap-gap_to_inner_gnd_front, -p.trace_width*0.5-p.prong_trace_gap+gap_to_inner_gnd_side),
               (p.trace_length+p.pad_trace_gap-p.prong_length-p.gap_front, -p.trace_width*0.5-p.prong_trace_gap+gap_to_inner_gnd_side),
               (p.trace_length+p.pad_trace_gap-p.prong_length-p.gap_front, -p.trace_width*0.5-p.prong_trace_gap-p.prong_width-p.gap_side),
               (p.trace_length+p.pad_trace_gap+p.pad_thickness+p.gap_back, -p.trace_width*0.5-p.prong_trace_gap-p.prong_width-p.gap_side)]

        # convert to shapely polygons for further processing (filleting, rotation, translation)
        pad1 = shapely.Polygon(pad1) 
        gap1 = shapely.Polygon(gap1)
        pad2 = shapely.Polygon(pad2)
        gap2 = shapely.Polygon(gap2)
        pin = shapely.LineString([(p.trace_length+p.pad_trace_gap+p.pad_thickness, p.trace_width*0.5),
                                  (p.trace_length+p.pad_trace_gap+p.pad_thickness, -p.trace_width*0.5)])
        
        def smooth_polygon(gap_poly, fillet_radius):
         # Apply create_trapezoid-style filleting to ANY polygon 
         raw_data = np.array(gap_poly.exterior.coords)  # Get raw coordinates of the polygon exterior (including duplicate last point)
    
         # Add start/end points like trapezoid (fixes first/last corner!)
         end_x, end_y = raw_data[1]        # Use polygon end as 2nd coordinate to overlap the first corner
         end = [[end_x, end_y]]
    
         # Stack: polygon + end (to ensure fillet is applied to first and last corner)
         raw_data_2 = np.vstack((raw_data,end))
    
         # Fillet ALL corners perfectly
         data = QUtilities.calc_filleted_path(raw_data_2, fillet_radius, p.fillet_resolution)
         data = data[2:-1]  # Remove the artificially added start/end points after filleting
         return Polygon(data)

        #  =========Filleting and smoothing operations=========


        if (p.fillet_radius == 0 and  # If no fillet radius is specified, skip all fillet operations to save time
           p.gap_fillet_radius == 0):
           # No smoothing requested - skip all fillet ops and proceed directly
           pass  # Polygons already ready
        else:
             if p.fillet_radius > 0:
                 pad1 = smooth_polygon(pad1, p.fillet_radius)
                 pad2 = smooth_polygon(pad2, p.fillet_radius)
             if p.gap_fillet_radius > 0:
                 gap1 = smooth_polygon(gap1, p.gap_fillet_radius)
                 gap2 = smooth_polygon(gap2, p.gap_fillet_radius)


        if p.swap_direction:
            polys = [pad1, pad2, gap1, gap2]
            polys = draw.translate(polys, -p.trace_length-p.pad_trace_gap-p.pad_thickness)
            polys = draw.rotate(polys, 180, origin=(0,0))
            [pad1, pad2, gap1, gap2] = polys

        polys = [pad1, pad2, gap1, gap2, pin]
        user_angle_deg = float(p.orientation) # Get user-specified orientation angle in degrees
        polys = draw.rotate(polys, user_angle_deg, origin=(0, 0))
        polys = draw.rotate(polys, np.arctan2(startPtNorm[1], startPtNorm[0]), origin=(0, 0), use_radians=True)
        polys = draw.translate(polys, *startPt)
        [pad1, pad2, gap1, gap2, pin] = polys

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(pad1=pad1, pad2=pad2),
                           layer=p.layer)

        #subtracts out ground plane on the layer it's on
        self.add_qgeometry('poly',
                           dict(padGap1=gap1, padGap2=gap2),
                           subtract=True,
                           layer=p.layer)
    
       # Generates its own pin
        self.add_pin('a', pin.coords, width=p.trace_width)
        

    '''












    # """Default drawing options"""
    component_metadata = Dict(
        short_name="Pocket",
        _qgeometry_table_path="True",
        _qgeometry_table_poly="True",
        _qgeometry_table_junction="True",
    )
    # """Component metadata"""

    # TOOLTIP = """Transmon pocket with tapering."""
    def __init__(self, design, name: str = None, options: Dict = None, **kwargs):
        super().__init__(design, name, options, **kwargs)

    def make(self):
        self.make_pocket()
        self.make_connection_pads()

    # Function to create trapezoid
    def create_trapezoid(
        self, center_x, top_width, base_width, height, y_offset, rfillet, startx
    ):
        half_top = float(top_width) / 2
        half_base = float(base_width) / 2
        coordinates = [
            (center_x - half_base, y_offset),
            (center_x - half_top, y_offset + height),
            (center_x + half_top, y_offset + height),
            (center_x + half_base, y_offset),
        ]

        # Create a Polygon object from the coordinates
        trapezoid_polygon = Polygon(coordinates)

        # Convert the coordinates to numpy array
        raw_data = np.array(trapezoid_polygon.exterior.coords)

        start_x = center_x - startx
        start_y = y_offset

        end_x = center_x + startx
        end_y = y_offset

        start = [[start_x, start_y]]
        end = [[end_x, end_y]]

        # Concatenate arrays after converting start and end points to 2D arrays
        raw_data_2 = np.vstack((start, raw_data, end))

        raw_data_1 = np.delete(raw_data_2, 5, 0)

        # Curve the edges of the polygon
        data = QUtilities.calc_filleted_path(raw_data_1, rfillet, 9)

        # Create a new Polygon object from the curved coordinates
        trapezoid_with_curved_edges = Polygon(data)

        return trapezoid_with_curved_edges

    # """Makes standard transmon in a pocket."""
    def make_pocket(self):

        # self.p allows us to directly access parsed values (string -> numbers) from the user option
        p = self.p
        #  pcop = self.p.coupled_pads[name]  # parser on connector options

        # since we will reuse these options, parse them once and define them as variables
        pad_width = p.pad_width
        pad_height = p.pad_height
        pad_gap = p.pad_gap
        coupled_pad_height = p.coupled_pad_height
        coupled_pad_width = p.coupled_pad_width
        coupled_pad_gap = p.coupled_pad_gap

        # TOP PAD
        # make the pads as rectangles (shapely polygons)
        pad = draw.rectangle(pad_width, pad_height)
        pad_top = draw.translate(pad, 0, +(pad_height + pad_gap) / 2.0)
        # Here, you make your pads round. Not sharp shape on the left and right sides and also this should be the same for the bottom pad as the top pad.
        circ_left_top = draw.Point(
            -pad_width / 2.0, +(pad_height + pad_gap) / 2.0
        ).buffer(pad_height / 2, resolution=16, cap_style=CAP_STYLE.round)
        circ_right_top = draw.Point(
            pad_width / 2.0, +(pad_height + pad_gap) / 2.0
        ).buffer(pad_height / 2, resolution=16, cap_style=CAP_STYLE.round)
        pad_top_tmp1 = draw.union([circ_left_top, pad_top, circ_right_top])
        # TAPER - Add trapezoid to the top pad
        if p.taper_height > 0:
            trapezoid_top = self.create_trapezoid(
                0,
                p.taper_width_top,
                p.taper_width_base,
                p.taper_height,
                (pad_gap / 2) - float(p.taper_height),
                p.taper_fillet_radius,
                p.pad_width / 2,
            )
            # pad_top_tmp = draw.union(pad_top_tmp1, trapezoid_top)
            trapezoid_top_rotated = draw.rotate(
                trapezoid_top,
                180,
                origin=(0, pad_gap / 4 + (pad_gap / 2 - float(p.taper_height)) / 2),
            )  # Rotate trapezoid by 180 degrees
            # create union
            pad_top_tmp = draw.union(
                pad_top_tmp1.buffer(0), trapezoid_top_rotated.buffer(0)
            )
        else:
            pad_top_tmp = pad_top_tmp1

        # COUPLER REGION - TOP PAD
        # In here you create the teeth part and then you union them as one with the pad. Teeth only belong to top pad.
        coupled_pad = draw.rectangle(coupled_pad_width, coupled_pad_height + pad_height)
        coupler_pad_round = draw.Point(
            0.0, (coupled_pad_height + pad_height) / 2
        ).buffer(coupled_pad_width / 2, resolution=16, cap_style=CAP_STYLE.round)
        coupled_pad = draw.union(coupled_pad, coupler_pad_round)
        coupled_pad_left = draw.translate(
            coupled_pad,
            -(coupled_pad_width / 2.0 + coupled_pad_gap / 2.0),
            +coupled_pad_height / 2.0 + pad_height + pad_gap / 2.0 - pad_height / 2,
        )
        coupled_pad_right = draw.translate(
            coupled_pad,
            (coupled_pad_width / 2.0 + coupled_pad_gap / 2.0),
            +coupled_pad_height / 2.0 + pad_height + pad_gap / 2.0 - pad_height / 2,
        )
        # The coupler pads are only created if low_W=0 and low_H=+1
        for name in self.options.connection_pads:
            if (
                self.options.connection_pads[name]["loc_W"] == 0
                and self.options.connection_pads[name]["loc_H"] == +1
            ):
                pad_top_tmp = draw.union(
                    [
                        circ_left_top,
                        coupled_pad_left,
                        pad_top_tmp,
                        coupled_pad_right,
                        circ_right_top,
                    ]
                )
        pad_top = pad_top_tmp
        # Curving the edges of the teeth where it joins the pad
        # outer corners
        pad_top = pad_top.buffer(p.fillet_radius, join_style=2, cap_style=3).buffer(
            -p.fillet_radius,
            cap_style=1,
            join_style=1,
            mitre_limit=2.0,
            quad_segs=p.fillet_resolution,
        )
        # inner corners
        pad_top = pad_top.buffer(-p.fillet_radius, join_style=2, cap_style=3).buffer(
            p.fillet_radius,
            cap_style=1,
            join_style=1,
            mitre_limit=2.0,
            quad_segs=p.fillet_resolution,
        )

        # BOTTOM PAD
        # Round part for the bottom pad. And again you should unite all of them.
        pad_bot = draw.translate(pad, 0, -(pad_height + pad_gap) / 2.0)
        circ_left_bot = draw.Point(
            -pad_width / 2, -(pad_height + pad_gap) / 2.0
        ).buffer(pad_height / 2, resolution=16, cap_style=CAP_STYLE.round)
        circ_right_bot = draw.Point(
            pad_width / 2, -(pad_height + pad_gap) / 2.0
        ).buffer(pad_height / 2, resolution=16, cap_style=CAP_STYLE.round)
        pad_bot = draw.union([pad_bot, circ_left_bot, circ_right_bot])
        # TAPER - Add trapezoid to the bottom pad
        if p.taper_height > 0:
            trapezoid_bot = self.create_trapezoid(
                center_x=0,
                top_width=p.taper_width_top,
                base_width=p.taper_width_base,
                height=p.taper_height,
                y_offset=-(pad_gap / 2),
                rfillet=p.taper_fillet_radius,
                startx=p.pad_width / 2,
            )
            pad_bot = draw.union(pad_bot, trapezoid_bot)

        # outer corners
        pad_bot = pad_bot.buffer(p.fillet_radius, join_style=2, cap_style=3).buffer(
            -p.fillet_radius,
            cap_style=1,
            join_style=1,
            mitre_limit=2.0,
            quad_segs=p.fillet_resolution,
        )
        # inner corners
        pad_bot = pad_bot.buffer(-p.fillet_radius, join_style=2, cap_style=3).buffer(
            p.fillet_radius,
            cap_style=1,
            join_style=1,
            mitre_limit=2.0,
            quad_segs=p.fillet_resolution,
        )
        # pad_bot = pad_bot.buffer(p.fillet_radius, join_style=1).buffer(-p.fillet_radius, join_style=1)

        ###############################################################################################
        # INSETS - cut into the pads for additional coupling
        # define trapezoid
        if p.inset_depth != 0 and p.inset_width != 0:
            inset = draw.rectangle(
                p.inset_depth - (2 * p.inset_fillet_radius),
                p.inset_width - (2 * p.inset_fillet_radius),
            )
            inset_rounded = inset.buffer(
                p.inset_fillet_radius,
                cap_style=1,
                join_style=1,
                mitre_limit=2.0,
                quad_segs=9,
            )
            inset = inset_rounded
            # define x and y locations for insets
            inset_y_top = (pad_gap / 2) + (p.pad_height / 2)
            inset_y_bot = -inset_y_top
            inset_x_left = (
                -(pad_width / 2)
                - (pad_height / 2)
                + (p.inset_depth / 2)
                - p.inset_fillet_radius
            )
            inset_x_right = -inset_x_left
            # position trapezoids
            inset_left_top = draw.translate(inset, inset_x_left, inset_y_top)
            inset_right_top = draw.translate(inset, inset_x_right, inset_y_top)
            inset_left_bot = draw.translate(inset, inset_x_left, inset_y_bot)
            inset_right_bot = draw.translate(inset, inset_x_right, inset_y_bot)
            # Subtract the insets from the pads
            pad_top_tmp1 = draw.subtract(pad_top, inset_left_top)
            pad_top_tmp2 = draw.subtract(pad_top_tmp1, inset_right_top)
            pad_top = pad_top_tmp2
            pad_bot_tmp1 = draw.subtract(pad_bot, inset_left_bot)
            pad_bot_tmp2 = draw.subtract(pad_bot_tmp1, inset_right_bot)
            pad_bot = pad_bot_tmp2
            # re-round the new pads (inner corners only)
            pad_bot = pad_bot.buffer(
                -p.inset_fillet_radius, join_style=2, cap_style=3
            ).buffer(
                p.inset_fillet_radius,
                cap_style=1,
                join_style=1,
                mitre_limit=2.0,
                quad_segs=p.fillet_resolution * 3,
            )
            pad_top = pad_top.buffer(
                -p.inset_fillet_radius, join_style=2, cap_style=3
            ).buffer(
                p.inset_fillet_radius,
                cap_style=1,
                join_style=1,
                mitre_limit=2.0,
                quad_segs=p.fillet_resolution * 3,
            )
            pad_bot = pad_bot.buffer(0)
            pad_top = pad_top.buffer(0)

        ###############################################################################################
        # Pins from the center of the qubit pads
        # Coordinates for the center of the pin
        taper_width_top = p.taper_width_top
        if p.taper_height == 0:
            pin_route_l = 1e-6
        else:
            pin_route_l = p.taper_height
        top_pin = LineString(
            [
                [0, pad_gap / 2 + p.junction_pin_inset],
                [0, pad_gap / 2 + p.junction_pin_inset - pin_route_l / 2],
                # Extend outward
            ]
        )
        bottom_pin = LineString(
            [
                [0, -pad_gap / 2 - p.junction_pin_inset],
                [0, -pad_gap / 2 - p.junction_pin_inset + pin_route_l / 2],  # Start from pad bottom edge
                # Extend outward
            ]
        )

        ###############################################################################################
        # Pins outside the qubit pocket
        # Define the pin positions relative to (0,0)
        top_pocket_pin = LineString(
            [
                [p.pocket_width / 2 + p.chrgln_pin_y_offset, p.chrgln_pin_x_offset],
                [
                    p.pocket_width / 2 + p.chrgln_pin_y_offset,
                    0,
                ],  # Start point (near pocket)
                # Extend outward
            ]
        )

        bottom_pocket_pin = LineString(
            [
                [-p.pocket_width / 2 - p.chrgln_pin_y_offset, p.chrgln_pin_x_offset],
                [
                    -p.pocket_width / 2 - p.chrgln_pin_y_offset,
                    0,
                ],  # Start point (near pocket)
                # Extend outward
            ]
        )

        ###############################################################################################
        # Josephson Junction
        rect_jj = draw.LineString(
            [(0, -p.inductor_height / 2), (0, p.inductor_height / 2)]
        )

        # Pocket
        rect_pk = draw.rectangle(
            p.pocket_width - 2 * p.fillet_radius_gap,
            p.pocket_height - p.pocket_lower_tighten - 2 * p.fillet_radius_gap,
            yoff=p.pocket_lower_tighten / 2,
        )

        # to curve the edges of the qubit pocket
        rect_pk = rect_pk.buffer(
            p.fillet_radius_gap,
            cap_style=1,
            join_style=1,
            mitre_limit=2.0,
            quad_segs=p.fillet_resolution,
        )

        # Rotate and translate all qgeometry as needed.
        polys = [
            rect_jj,
            pad_top,
            pad_bot,
            rect_pk,
            top_pin,
            bottom_pin,
            top_pocket_pin,
            bottom_pocket_pin,
        ]
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [
            rect_jj,
            pad_top,
            pad_bot,
            rect_pk,
            top_pin,
            bottom_pin,
            top_pocket_pin,
            bottom_pocket_pin,
        ] = polys

        # Add shapes as qiskit geometries
        self.add_qgeometry("poly", dict(pad_top=pad_top, pad_bot=pad_bot))
        self.add_qgeometry("poly", dict(rect_pk=rect_pk), subtract=True)
        self.add_qgeometry("junction", dict(rect_jj=rect_jj), width=p.inductor_width)

        # Pins from the center of the qubit pads
        self.add_pin(
            "pin_island",
            points=list(top_pin.coords),
            width=taper_width_top,
            input_as_norm=True,
        )
        self.add_pin(
            "pin_reservior",
            points=list(bottom_pin.coords),
            width=taper_width_top,
            input_as_norm=True,
        )

        # Add pins to the component
        self.add_pin(
            "bottom_pin", points=list(top_pocket_pin.coords), width=taper_width_top
        )
        self.add_pin(
            "top_pin",
            points=list(bottom_pocket_pin.coords),
            width=taper_width_top,
            input_as_norm=True,
        )

    def make_connection_pads(self):
        # """Makes standard transmon in a pocket."""
        for name in self.options.connection_pads:
            self.make_connection_pad(name)

    def make_connection_pad(self, name: str):
        # """Makes n individual connector.

        # Args:
        #     name (str) : Name of the connector
        # """

        # self.p allows us to directly access parsed values (string -> numbers) form the user option
        p = self.p
        pc = self.p.connection_pads[name]  # parser on connector options

        # define commonly used variables once
        r = pc.fillet_radius_inner
        r_outer = pc.fillet_radius_outer
        cpw_width = pc.cpw_width
        cpw_extend = pc.cpw_extend
        pad_width_fillet = pc.pad_width - 2 * r
        pad_height_fillet = pc.pad_height - 2 * r
        pad_width = pc.pad_width
        pad_height = pc.pad_height
        pad_cpw_shift = pc.pad_cpw_shift
        pocket_rise = pc.pocket_rise
        pocket_extent = pc.pocket_extent

        assert (
            pad_width_fillet >= 0
        ), f"Error: pad_width {pad_width} is too small for the fillet radius {r}. Either increase the width or decrease the fillet radius of the connection pads."
        assert (
            pad_height_fillet >= 0
        ), f"Error: pad_height {pad_height} is too small for the fillet radius {r}. Either increase the height or decrease the fillet radius of the connection pads."
        assert (
            pad_width / 2 - cpw_width / 2 - r >= r_outer
        ), f"Error: fillet_radius_outer {r_outer} is too large for the pad_width {pad_width} and cpw_width {cpw_width}. Either decrease the fillet_radius_outer or increase the pad_width."

        loc_W = float(pc.loc_W)
        loc_W, loc_H = float(pc.loc_W), float(pc.loc_H)
        if float(loc_W) not in [-1.0, +1.0, 0] or float(loc_H) not in [-1.0, +1.0]:
            self.logger.info(
                "Warning: Did you mean to define a transmon qubit with loc_W and"
                " loc_H that are not +1, -1, or 0? Are you sure you want to do this?"
            )

        # Define the geometry
        # Connector pad
        if float(loc_W) != 0:
            connector_pad = draw.rectangle(
                pad_width_fillet, pad_height_fillet, -pad_width / 2, pad_height / 2
            )
            connector_pad = connector_pad.buffer(
                r,
                cap_style=1,
                join_style=1,
                mitre_limit=2.0,
                quad_segs=pc.fillet_resolution,
            )
            # add corners
            if r_outer > 0:
                print("Skipping: connector pad rounding not yet implemented for W != 0")

            # Connector CPW wire
            connector_wire_path = draw.wkt.loads(
                f"""LINESTRING (\
                0 {pad_cpw_shift+cpw_width/2}, \
                {pc.pad_cpw_extent}                           {pad_cpw_shift+cpw_width/2}, \
                {(p.pocket_width-p.pad_width)/2-pocket_extent} {pad_cpw_shift+cpw_width/2+pocket_rise}, \
                {(p.pocket_width-p.pad_width)/2+cpw_extend}    {pad_cpw_shift+cpw_width/2+pocket_rise}\
                                            )"""
            )
        else:
            connector_pad = draw.rectangle(
                pad_width_fillet, pad_height_fillet, 0, pad_height / 2
            )
            connector_pad = connector_pad.buffer(
                r,
                cap_style=1,
                join_style=1,
                mitre_limit=2.0,
                quad_segs=pc.fillet_resolution,
            )
            # add corners
            if r_outer > 0:
                # --- Add the small corner box first ---
                connector_pad_corners_sq = draw.rectangle(
                    r_outer,
                    r_outer,
                    -(cpw_width / 2) - (r_outer / 2),
                    pad_height + (r_outer / 2),
                )
                # --- Create a quarter circle to cut the outer corner ---
                corner_circle = draw.Point(
                    -(cpw_width / 2) - (r_outer), pad_height + r_outer
                ).buffer(r_outer, resolution=32)
                corner_subtract_l = draw.subtract(
                    connector_pad_corners_sq, corner_circle
                )
                # repeat on right side
                connector_pad_corners_sq = draw.rectangle(
                    r_outer,
                    r_outer,
                    (cpw_width / 2) + (r_outer / 2),
                    pad_height + (r_outer / 2),
                )
                corner_circle = draw.Point(
                    (cpw_width / 2) + (r_outer), pad_height + r_outer
                ).buffer(r_outer, resolution=32)
                corner_subtract_r = draw.subtract(
                    connector_pad_corners_sq, corner_circle
                )
                # Union the corners
                connector_pad = draw.union(
                    [connector_pad, corner_subtract_l, corner_subtract_r]
                )
            # CPW path
            connector_wire_path = draw.LineString(
                [
                    [0, pad_height],
                    [
                        0,
                        (p.pocket_width / 2 - p.pad_height - p.pad_gap / 2 - pc.pad_gap)
                        + cpw_extend,
                    ],
                ]
            )

        # Position the connector, rotate and translate
        objects = [connector_pad, connector_wire_path]

        if loc_W == 0:
            loc_Woff = 1
        else:
            loc_Woff = loc_W

        objects = draw.scale(objects, loc_Woff, loc_H, origin=(0, 0))
        objects = draw.translate(
            objects,
            loc_W * (p.pad_width) / 2.0,
            loc_H * (p.pad_height + p.pad_gap / 2 + pc.pad_gap),
        )
        objects = draw.rotate_position(objects, p.orientation, [p.pos_x, p.pos_y])
        [connector_pad, connector_wire_path] = objects

        self.add_qgeometry("poly", {f"{name}_connector_pad": connector_pad})
        self.add_qgeometry(
            "path", {f"{name}_wire": connector_wire_path}, width=cpw_width
        )
        self.add_qgeometry(
            "path",
            {f"{name}_wire_sub": connector_wire_path},
            width=cpw_width + 2 * pc.cpw_gap,
            subtract=True,
        )

        ############################################################

        # add pins
        points = np.array(connector_wire_path.coords)
        self.add_pin(name, points=points[-2:], width=cpw_width, input_as_norm=True)

    def get_resonator_length_mm(self, connection_pad="readout"):
        """
        Gives total length of resonator segment.
        """
        return QUtilities.calc_points_on_path([0], self.design, component_name=self.name, trace_name=f"{connection_pad}_wire")[-1]'''