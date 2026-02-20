
### This file is part of the SQDMetal project: Test to create a new Capacitor class
### Author: Léon Thinus

import numpy as np
from qiskit_metal import draw, Dict
from qiskit_metal.qlibrary.core import BaseQubit
from shapely.geometry.base import CAP_STYLE
from shapely.geometry import Polygon
from SQDMetal.Utilities.QUtilities import QUtilities 
from shapely.geometry import LineString
from shapely.geometry import Point
from qiskit_metal.toolbox_python.attr_dict import Dict
from qiskit_metal.qlibrary.core import QComponent
import shapely
from SQDMetal.Utilities.CpwParams import CpwParams
from SQDMetal.Comps.Resonators import ResonatorMeander
from SQDMetal.Utilities.QiskitShapelyRenderer import QiskitShapelyRenderer
from SQDMetal.Utilities.ShapelyEx import ShapelyEx

class Smooth_CapacitorUcapGroundPin(QComponent):
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
        

class Smooth_CapacitorProngPin(QComponent):
    """Creates a forked gap capacitor that wraps around the target lead. Take already existing code and add fillets to the edges of the capacitor and the ground plane cutouts.

    Inherits QComponent class.

    Capacitor Metal Geometry and Ground Cutout Pocket:
        * prong_width  - Width of the two prongs
        * prong_length - Length of the two prongs
        * pin_gap      - Capacitor gap from the fork to the target pin
        * pin_gap_side - Capacitor gap of side prongs to the sides of the target pin
        * pad_thickness - Length of the Capacitor fork pad

    Spacing around the structure (i.e. cuts into the ground plane) can be controlled via:
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
        * orientation = 0

        ## Fillet options for the prongs and ground plane cutouts:
        * fillet_radius = '5um' - Radius of edge fillets. If 0, no fillets are created.
        * fillet_resolution = 16 - Resolution of the fillets in number of points. Higher number means smoother fillets but longer simulation time.
        * gap_fillet_radius = '5um' - Fillet radius for the cuts in the ground plane around the capacitor. If 0, no fillets are created.
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
                           trace_width='10um',
                            orientation = 0,
                            fillet_radius = '5um',
                            fillet_resolution = 16,
                            gap_fillet_radius = '5um')
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
                 pad = smooth_polygon(pad, p.fillet_radius)
                 
             if p.gap_fillet_radius > 0:
                 gap = smooth_polygon(gap, p.gap_fillet_radius)

        polys = [pad, gap, pin]
        user_angle_deg = float(p.orientation) # Get user-specified orientation angle in degrees
        polys = draw.rotate(polys, user_angle_deg, origin=(0, 0))
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

class Smooth_Capacitor_Semicircle(QComponent):
    """Creates a capacitor design with a semicircular pad and a prong that wraps around the target lead.

    Inherits QComponent class.

    Original rectangle capacitor design:
        * Rect_width  - Width of the rectangular pad
        * Rect_length - Length of the rectangular pad
        * orientation - Angle of the capacitor in degrees

    Spacing around the structure (i.e. cuts into the ground plane) can be controlled via:
        * gap_side  - Spacing to ground plane on the outer sides 
        * gap_front - Spacing to ground plane in the direction of the target pin
        * gap_back - Spacing to ground plane away from the direction of the target pin

    Semicircle geometry:
        * semi_radius - Radius of the semicircular pad
        * circle_offset - Offset of the center of the semicircle from the end of the prong in the direction of the target pin.


    The positioning can be done dynamically via:
        * pin_inputs=Dict(start_pin=Dict(component=f'...',pin='...')) - Specifying the target component pin

    Pins:
        There is one pin 'a' to link to the fork (there is no trace drawn in this structure intrinsically). Only needs a width defined
        * trace_width - center trace width of the trace that attaches to this capacitor.

    Sketch:
        Below is a sketch of the capacitor
        ::

            @@@@@@@@@@@@@@@@@@@     @  = Ground Plane
            @  W W W W W W GB @     #  = Target pin
            @   _____x_____GB @     x = pin location with width trace_width
            @ L|           |  @     W = rect_width
            @ L|     _     |  @     L = rect_length
            @ L|   /   \   |GS@     r = semi_radius
            @ L| /       \ |  @    
            @ L|/     rrrr\|  @     
            @       ###   GF  @     BG = gap_back
            @       ###   GF  @     GF = gap_front
            @@@@@@  ###  @@@@@@     GS = gap_side

    .. image::
        Cap3Interdigital.png

    .. meta::
        Cap 3 Interdigital

    Default Options:
        * rect_width='50um'
        * rect_length='50um'
        * gap_side='5um'
        * gap_front='5um'
        * gap_back='5um'
        * semi_radius='25um'
        * circle_offset='0um'
        * orientation = 0

        ## Fillet options for the prongs and ground plane cutouts:
        * fillet_radius = '5um' - Radius of edge fillets. If 0, no fillets are created.
        * fillet_resolution = 16 - Resolution of the fillets in number of points. Higher number means smoother fillets but longer simulation time.
        * gap_fillet_radius = '5um' - Fillet radius for the cuts in the ground plane around the capacitor. If 0, no fillets are created.
    """

    #  Define structure functions

    default_options = Dict(rect_width='50um',
                           rect_length='50um',
                           semi_radius='25um',
                           circle_offset='0um',
                           gap_side='5um',
                           gap_front='5um',
                           gap_back='10um',
                           trace_width='10um',
                            orientation = 0,
                            fillet_radius = '5um',
                            fillet_resolution = 16,
                            gap_fillet_radius = '5um')
    """Default drawing options"""

    TOOLTIP = """Create a rectangular capacitor pad with a ground pocket cutout."""

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

        #wid = start_point['width']
            #define the polygons for the basic capacitor structure (before filleting and rotation)
            # Rectangular pad 
        pad = [(0,p.rect_width*0.5),
                (0,-p.rect_width*0.5),
                (p.rect_length, -p.rect_width*0.5),
                (p.rect_length, p.rect_width*0.5)
              ]

        gap = [(-p.gap_back,p.rect_width*0.5+p.gap_side),
                (-p.gap_back,-p.rect_width*0.5-p.gap_side),
                (p.rect_length+p.gap_front, -p.rect_width*0.5-p.gap_side),
                (p.rect_length+p.gap_front, p.rect_width*0.5+p.gap_side)
              ]
        
        circle_center_x = p.rect_length + p.circle_offset
        circle_center_y = 0
        circle = Point(circle_center_x, circle_center_y).buffer(p.semi_radius)

        pad = shapely.Polygon(pad[::-1])
        pad = draw.subtract(pad, circle)
        gap = shapely.Polygon(gap)
        pin = shapely.LineString([(p.gap_front, p.rect_width*0.5),
                                  (p.gap_front, -p.rect_width*0.5)])
        
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
                 pad = smooth_polygon(pad, p.fillet_radius)
                 
             if p.gap_fillet_radius > 0:
                 gap = smooth_polygon(gap, p.gap_fillet_radius)

        polys = [pad, gap, pin]
        user_angle_deg = float(p.orientation) # Get user-specified orientation angle in degrees
        polys = draw.rotate(polys, user_angle_deg, origin=(0, 0))
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