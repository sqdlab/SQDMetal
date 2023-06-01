from qiskit_metal.qlibrary.core import QComponent, QRoute
from qiskit_metal.qlibrary.tlines.anchored_path import RouteAnchors
from qiskit_metal import draw
import numpy as np
import shapely
from qiskit_metal.toolbox_python.attr_dict import Dict
from SQDMetal.Utilities.QUtilities import QUtilities

class WirePinStretch(QRoute):
    default_options = Dict(dist_extend='10um')

    def __init__(self, design,
                        name: str = None,
                        options: Dict = None,
                        type: str = "CPW",
                        **kwargs):
        #QRoute forces an end-pin to exist... So make it artificial...
        options.pin_inputs['end_pin'] = Dict(component=options.pin_inputs.start_pin.component, pin=options.pin_inputs.start_pin.pin)
        options.fillet='0um'
        super().__init__(design, name, options, **kwargs)

    def make(self):
        p = self.p

        self.set_pin("start")
        
        start_point = self.design.components[self.options.pin_inputs.start_pin.component].pins[self.options.pin_inputs.start_pin.pin]
        startPt = start_point['middle']
        norm = start_point['normal']
        endPt = startPt + norm*p.dist_extend

        self.add_pin('end', [endPt-norm, endPt], width=p.trace_width, input_as_norm=True)
        tang = np.array([-norm[1],norm[0]])
        pt = endPt+tang*p.trace_width*0.5-norm*p.trace_width*0.5
        self.add_pin('endL', [pt-tang, pt], width=p.trace_width, input_as_norm=True)
        pt = endPt-tang*p.trace_width*0.5-norm*p.trace_width*0.5
        self.add_pin('endR', [pt+tang, pt], width=p.trace_width, input_as_norm=True)
        
        self.make_elements(np.array([startPt, endPt]))

class WirePins(QRoute):
    default_options = Dict(pathObjPins=[],
                           trace_width='10um', trace_gap='10um', fillet='50um', end_gap='0um',
                           advanced=Dict(avoid_collision='false'))
    
    def __init__(self, design,
                        name: str = None,
                        options: Dict = None,
                        type: str = "CPW",
                        **kwargs):
        #QRoute forces an end-pin to exist... So make it artificial...
        pins = []
        for cur_pin in [options.pathObjPins[0], options.pathObjPins[-1]]:
            if isinstance(cur_pin, list) or isinstance(cur_pin, tuple):
                comp, pin = cur_pin[0], cur_pin[1]
            else:
                comp, pin = cur_pin, 'a'
            pins += [(comp, pin)]
        assert len(pins) > 0, "Must have at least 2 pins/waypoints."
        options.pin_inputs['start_pin'] = Dict(component=pins[0][0], pin=pins[0][1])
        options.pin_inputs['end_pin'] = Dict(component=pins[1][0], pin=pins[1][1])
        super().__init__(design, name, options, **kwargs)

    def make(self):
        p = self.p

        #Parse Path Joints...
        pins = []
        for cur_pin in self.options.pathObjPins:
            if isinstance(cur_pin, list) or isinstance(cur_pin, tuple):
                comp, pin = cur_pin[0], cur_pin[1]
            else:
                comp, pin = cur_pin, 'a'
            pins += [self.design.components[comp].pins[pin]['middle']]
        
        line = shapely.LineString(np.vstack(pins))
        if p.end_gap > 0:
            norm_vec = pins[-1]-pins[-2]
            norm_vec = norm_vec / np.linalg.norm(norm_vec) * p.end_gap
            pins[-1] = pins[-1] + norm_vec
            line_gap = shapely.LineString(np.vstack(pins))
        else:
            line_gap = line

        self.add_qgeometry('path', {'trace': line},
                           width=p.trace_width,
                           fillet=p.fillet,
                           layer=p.layer)
        self.add_qgeometry('path', {'cut': line_gap},
                               width=p.trace_width + 2 * p.trace_gap,
                               fillet=p.fillet,
                               layer=p.layer,
                               subtract=True)

class WireElbowParallelPinPin(QRoute):
    """A rudimentary wire connector that mimics the 'Microsoft Office Autoshapes Elbow Connector with yellow handle' with curved corners.
    The connecting pins must be either parallel or anti-parallel

    Inherits QComponent class.

    Resonator Metal Geometry and Ground Cutout Pocket:
        * frac_pos_elbow - Position of the crossward part of the elbow; only exists when the pins are anti-parallel. It is a fraction where
                           0.0 is next to the start pin while 1.0 is next to the end pin.
        * fillet - Radius of the turns on the elbows
        * pin_pad - Minimum distance to traverse before turning from a given pin (making it non-zero helps ensure it passes the 'bad-fillet'
                    checks in Qiskit-Metal). Ignored in the 2 Elbows, anti-parallel case shown below...
    The parameters trace_width and trace_gap define the CPW dimensions.

    The positioning can be done dynamically via:
        * pin_inputs=Dict(start_pin=Dict(component=f'...',pin='...'), end_pin=Dict(component=f'...',pin='...')) - Specifying start and end
          positions via a component pins

    Sketch:
        Below is a sketch of the wiring configurations available depending on the position of the pins
        ::
                     P  R                     f                             P  R            R = fillet
                   p1>----\            p1>----\ R                       p1>----\            P = pin_pad
                          |                   |                                |            f = frac_pos_elbow
                          |                   |                   /------------/ f
                          |                   |                   |
                 p2>------/                   \---<p2             \---<p2
            Parallel, 2 Elbows     Anti-parallel, 2 Elbows     Anti-parallel, 4 Elbows

        Notice that frac_pos_elbow is only relevant in the anti-parallel pin configurations. In the 2-elbow case, it controls the position of
        the long vertical edge in the direction of the pins, while in the 4-elbow cases, it controls the direction of the long horizontal edge
        in the direction perpendicular to the pins.

    .. image::
        Cap3Interdigital.png

    .. meta::
        Wire with Elbow between colinear pins

    Default Options:
        * frac_pos_elbow=0.5
        * fillet='20um'
        * pin_pad='2um'
    """

    default_options = Dict(frac_pos_elbow=0.5, fillet='20um', pin_pad='2um')

    def make(self):
        p = self.p

        start_point = self.design.components[self.options.pin_inputs.start_pin.component].pins[self.options.pin_inputs.start_pin.pin]
        startPt = start_point['middle']
        start_norm = start_point['normal']

        end_point = self.design.components[self.options.pin_inputs.end_pin.component].pins[self.options.pin_inputs.end_pin.pin]
        endPt = end_point['middle']
        end_norm = end_point['normal']

        assert np.abs(start_norm[0]*end_norm[1] - start_norm[1]*end_norm[0]) < 1e-15, "The two pins must be facing either parallel or anti-parallel."

        matRot = np.array([[start_norm[0], start_norm[1]], [-start_norm[1], start_norm[0]]])
        start_norm = matRot @ start_norm
        end_norm = matRot @ end_norm
        endPt = matRot @ (endPt-startPt)

        if end_norm[0] > 0:
            #2 Elbows - parallel
            if endPt[0] > 0:
                path = [[0,0], [endPt[0]+p.pin_pad+p.fillet, 0], [endPt[0]+p.pin_pad+p.fillet, endPt[1]], endPt]
            else:
                path = [[0,0], [p.pin_pad+p.fillet,0], [p.pin_pad+p.fillet, endPt[1]], endPt]
        elif endPt[0] >= 2*(p.pin_pad + p.fillet):
            #2 Elbows - anti-parallel
            path = [[0,0], [endPt[0]*p.frac_pos_elbow, 0], [endPt[0]*p.frac_pos_elbow, endPt[1]], endPt]
        else:
            #4 Elbows - anti-parallel            
            path = [[0,0], [p.pin_pad+p.fillet, 0],
                    [p.pin_pad+p.fillet, endPt[1]*p.frac_pos_elbow], [endPt[0]-p.pin_pad-p.fillet, endPt[1]*p.frac_pos_elbow],
                    [endPt[0]-p.pin_pad-p.fillet, endPt[1]], endPt]
        
        path = matRot.T @ np.array(path).T
        path = path.T
        path[:,0] += startPt[0]
        path[:,1] += startPt[1]

        self.set_pin("start")
        self.set_pin("end")
        self.make_elements(path)


class WireTaperPin(QComponent):
    #  Define structure functions

    default_options = Dict(trace_width='20um',
                           trace_gap='10um',
                           taper_length='30um',
                           orig_gap='10um')
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
        rot_angle = np.arctan2(norm[1], norm[0])
        width = start_point['width']
        if 'trace_gap' in self.design.components[self.options.pin_inputs.start_pin.component].options:
            gap = self.design.components[self.options.pin_inputs.start_pin.component].options['trace_gap']
            gap = QUtilities.parse_value_length(gap) / QUtilities.get_units(self.design)
        else:
            gap = p.orig_gap

        taper = [
                (0, width*0.5),
                (p.taper_length, p.trace_width*0.5),
                (p.taper_length, -p.trace_width*0.5),
                (0, -width*0.5)]

        taper_gap = [
                    (0, width*0.5+gap),
                    (p.taper_length, p.trace_width*0.5+p.trace_gap),
                    (p.taper_length, -p.trace_width*0.5-p.trace_gap),
                    (0, -width*0.5-gap)]

        pin = shapely.LineString(taper[1:3])
        taper = shapely.Polygon(taper)
        taper_gap = shapely.Polygon(taper_gap)

        polys = [taper, taper_gap, pin]
        polys = draw.rotate(polys, rot_angle, origin=(0, 0), use_radians=True)
        polys = draw.translate(polys, *startPt)
        [taper, taper_gap, pin] = polys
        
        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(taper=taper),
                           layer=p.layer)

        #subtracts out ground plane on the layer it's on
        self.add_qgeometry('poly',
                           dict(taper_gap=taper_gap),
                           subtract=True,
                           layer=p.layer)

        # Generates its own pins
        self.add_pin('a', pin.coords, width=p.trace_width)
