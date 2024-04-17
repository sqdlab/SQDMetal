# -*- coding: utf-8 -*-
# Author: Prasanna Pakkiam
# Creation Date: 01/05/2023
# Description: Collection of classes to wires - from routing objects to tapers.

from qiskit_metal.qlibrary.core import QComponent, QRoute
from qiskit_metal.qlibrary.tlines.anchored_path import RouteAnchors
from qiskit_metal import draw
import numpy as np
import shapely
from qiskit_metal.toolbox_python.attr_dict import Dict
from SQDMetal.Utilities.QUtilities import QUtilities

class WirePinStretch(QRoute):
    """Creates a wire that extends some distance away from some target component's pin. It also has provisions to provide
    gaps in the start and/or end to the ground plane (e.g. for open-ended terminations).

    Inherits QComponent class.

    The positioning can be done dynamically via:
        * pin_inputs=Dict(start_pin=Dict(component=f'...',pin='...')) - Specifying start position via a component pin
        * dist_extend - Distance upon to stretch away from the start pin.
    The resulting wire is a straight section with width equal to that specified by the start pin and extends dist_extend
    from the start pin (in the direction given by the start pin).

    The part additionally has options to add ground cuts to the start and end portions (e.g. open-ended terminations):
        * start_gap - Distance to cut into ground plane at start.
        * end_gap   - Distance to cut into ground plane at end.
    The values are interpreted as follows:
        - start_gap, end_gap > 0, the gap to the ground plane is the supplied value
        - start_gap, end_gap == 0, the gap to the ground plane is taken as trace_gap
        - start_gap, end_gap < 0, no gap is given.        

    As usual trace_width and trace_gap specify the CPW dimensions.

    Pins:
        There are three pins in the end section:
            * end
            * end_L
            * end_R
        where end is flush in the direction of the wire, while end_L and end_R are on the sides:
        
                        <-w->             x = start pin location 
         _________________l__             e = 'end' pin
        |                    |  /|\       l = 'end_L' pin
        x                    e   w        r = 'end_R' pin
        |_________________r__|  \|/       d = dist_extend
        <----------d--------->            w = width of start pin

    .. image::
        Cap3Interdigital.png

    .. meta::
        Wire that extends off a pin.

    Default Options:
        * dist_extend='10um'
        * start_gap='-1um'
        * end_gap='-1um'
    """

    default_options = Dict(dist_extend='10um',
                           start_gap='-1um',
                           end_gap='-1um')

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
        
        line = shapely.LineString(np.array([startPt, endPt]))
        line_gap = np.array([startPt, endPt])
        if p.start_gap > 0:
            line_gap[0] -= norm*p.start_gap
        elif p.start_gap == 0:
            line_gap[0] -= norm*p.trace_gap
        if p.end_gap > 0:
            line_gap[1] += norm*p.end_gap
        elif p.end_gap == 0:
            line_gap[1] += norm*p.trace_gap
        line_gap = shapely.LineString(line_gap)

        self.add_qgeometry('path', {'trace': line},
                           width=p.trace_width,
                           fillet=p.fillet,
                           layer=p.layer)
        self.add_qgeometry('path', {'cut': line_gap},
                               width=p.trace_width + 2 * p.trace_gap,
                               fillet=p.fillet,
                               layer=p.layer,
                               subtract=True)
        # self.make_elements(np.array([startPt, endPt]))

class WirePins(QRoute):
    """Rudimentary manual routing wire that passes through all the listed pins. It supports CPW cutouts and an additional
    gap in the end for open-ended terminations.

    Inherits QComponent class.

    The CPW wire parameters are given via:
        * trace_width - width of CPW
        * trace_gap - gap to ground plane around CPW
        * fillet    - radius of fillet/turns on every corner
        * end_gap   - gap to ground plane at the end of the wire (useful for open-ended terminations)

    The waypoints for the positioning/routing is given via a list of pins in 'pathObjPins'. Each pin can be either:
        * 2-tuple signifying the component name and the pin name
        * Single string - in this case, it is presumed that the pin name is 'a' (as is the case with most Joints)
    For example:
        pathObjPins=[('Launcher1', 'tie'), 'jnt_1', 'jnt_2', ('Launcher2', 'tie')]

    Pins:
        There are no inherent pins in this object, but they can be added manually via constructs like RouteJoint

    .. image::
        Cap3Interdigital.png

    .. meta::
        Wire that hits manually selected pins.

    Default Options:
        * pathObjPins=[]    <--- Needs to have at least 2 pins to form a path!
        * trace_width='10um'
        * trace_gap='10um'
        * fillet='50um'
        * end_gap='0um'
    """

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

        self.set_pin("start")
        self.set_pin("end")
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

    Inherits QRoute class.

    Wire Metal Geometry and Ground Cutout Pocket:
        * frac_pos_elbow - Position of the crossward part of the elbow; only exists when the pins are anti-parallel. It is a fraction where
                           0.0 is next to the start pin while 1.0 is next to the end pin.
        * fillet - Radius of the turns on the elbows
        * pin_pad - Minimum distance to traverse before turning from a given pin (making it non-zero helps ensure it passes the 'bad-fillet'
                    checks in Qiskit-Metal). Ignored in the 2 Elbows, anti-parallel case shown below...
        * end_gap   - gap to ground plane at the end of the wire (useful for open-ended terminations)
    The parameters trace_width and trace_gap define the CPW dimensions.

    The positioning can be done dynamically via:
        * pin_inputs=Dict(start_pin=Dict(component=f'...',pin='...'), end_pin=Dict(component=f'...',pin='...')) - Specifying start and end
          positions via a component pins

    Sketch:
        Below is a sketch of the wiring configurations available depending on the position of the pins
        ::
                       P  R                   f                             P  R            R = fillet
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
        * end_gap='0um'
    """

    default_options = Dict(frac_pos_elbow=0.5, fillet='20um', pin_pad='2um', end_gap='0um')

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

        line = shapely.LineString(np.vstack(path))
        if p.end_gap > 0:
            norm_vec = path[-1]-path[-2]
            norm_vec = norm_vec / np.linalg.norm(norm_vec) * p.end_gap
            path[-1] = path[-1] + norm_vec
            line_gap = shapely.LineString(np.vstack(path))
        else:
            line_gap = line

        self.set_pin("start")
        self.set_pin("end")
        self.add_qgeometry('path', {'trace': line},
                           width=p.trace_width,
                           fillet=p.fillet,
                           layer=p.layer)
        self.add_qgeometry('path', {'cut': line_gap},
                               width=p.trace_width + 2 * p.trace_gap,
                               fillet=p.fillet,
                               layer=p.layer,
                               subtract=True)


class WireElbowSingle(QRoute):
    """A wire connection that has one corner/elbow turn to connect two pins.

    Inherits QRoute class.

    Wire Metal Geometry and Ground Cutout Pocket:
        * fillet - Radius of the turns on the elbows
    The parameters trace_width and trace_gap define the CPW dimensions.

    The positioning can be done dynamically via:
        * pin_inputs=Dict(start_pin=Dict(component=f'...',pin='...'), end_pin=Dict(component=f'...',pin='...')) - Specifying start and end
          positions via a component pins

    Sketch:
        Below is a sketch of the wiring configurations available depending on the position of the pins
        ::
            P1###########
                         #
                          #             # = Wire
                           #            P1 and P2 = the pins between which the wire connects
                            #
                             P2

        Notice that the pin directions point in direction of wire segments. The turning corner is set at the line intersection of the two pin
        vectors. Thus, this shouldn't be used when the pins are facing away from one another.

    .. image::
        Cap3Interdigital.png

    .. meta::
        Wire with single corner between pins

    Default Options:
        * trace_width='10um'
        * trace_gap='10um'
        * fillet='20um'
    """

    default_options = Dict(trace_width='10um',
                           trace_gap='10um',
                           fillet='20um')

    def make(self):
        p = self.p

        start_point = self.design.components[self.options.pin_inputs.start_pin.component].pins[self.options.pin_inputs.start_pin.pin]
        startPt = start_point['middle']
        start_norm = start_point['normal']

        end_point = self.design.components[self.options.pin_inputs.end_pin.component].pins[self.options.pin_inputs.end_pin.pin]
        endPt = end_point['middle']
        end_norm = end_point['normal']
        
        try:
            t,u = np.linalg.solve([[start_norm[0], end_norm[0]], [start_norm[1], end_norm[1]]], [endPt[0]-startPt[0], endPt[1]-startPt[1]])
        except:
            assert False, "The pins are parallel and this construct should not be used. Try for example, WireElbowParallelPinPin or WirePins instead."
        assert t > 0, "The pins are facing away and will cause the turn to occur at some point behind the components. Check direction of pins."
        midPt = startPt+start_norm*t

        path = [startPt, midPt, endPt]

        self.set_pin("start")
        self.set_pin("end")
        self.make_elements(path)

class WireTaperPinStretch(QComponent):
    """A taper to change width of wire at some pin.

    Inherits QComponent class.

    Wire Metal Geometry and Ground Cutout Pocket:
        * dist_extend - Length of taper
        * orig_gap    - If the target pin is a CPW, this argument will be ignored and the trace_gap parameter from the pin's CPW component
                        will be taken. Otherwise, the initial CPW gap will be taken from orig_gap.
    The parameters trace_width and trace_gap define the final CPW dimensions.

    The positioning can be done dynamically via:
        * pin_inputs=Dict(start_pin=Dict(component=f'...',pin='...')) - Specifying start and end
          positions via a component pins

    Pins:
        There is a pin called 'a' to continue the taper onto a wire of different width.

    Sketch:
        Below is a sketch of the wiring configurations available depending on the position of the pins
        ::
              @@@@@@G
              @@@@  G           # = Taper
              @@    G           @ = Ground plane
             g    ## /|\        P = Pin to attach taper
             g  ####  |         W = trace_width
              ######  |         G = trace_gap
            P ######  W         D = dist_extend
              ######  |         g = orig_gap
                ####  |
                  ## \|/
              <..D.>
              @@@@
              @@@@@@

        Notice that the pin directions point in direction of wire segments. The turning corner is set at the line intersection of the two pin
        vectors. Thus, this shouldn't be used when the pins are facing away from one another.

    .. image::
        Cap3Interdigital.png

    .. meta::
        Taper to change width of wire at some pin.

    Default Options:
        * trace_width='20um'
        * trace_gap='10um'
        * dist_extend='30um'
        * orig_gap='10um'
    """
    default_options = Dict(trace_width='20um',
                           trace_gap='10um',
                           dist_extend='30um',
                           orig_gap='10um')

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
                (p.dist_extend, p.trace_width*0.5),
                (p.dist_extend, -p.trace_width*0.5),
                (0, -width*0.5)]

        taper_gap = [
                    (0, width*0.5+gap),
                    (p.dist_extend, p.trace_width*0.5+p.trace_gap),
                    (p.dist_extend, -p.trace_width*0.5-p.trace_gap),
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

class WireTaperProbePinStretch(QComponent):
    """A taper that splits into two probes.

    Inherits QComponent class.

    Wire Metal Geometry and Ground Cutout Pocket:
        * dist_extend - Length of taper
        * orig_gap    - If the target pin is a CPW, this argument will be ignored and the trace_gap parameter from the pin's CPW component
                        will be taken. Otherwise, the initial CPW gap will be taken from orig_gap.
        * probe_gap1  - Gap between the two probes in the taper at origin
        * probe_gap2  - Gap between the two probes in the taper at the end of the extension
        
    The parameters trace_width and trace_gap define the final CPW dimensions.

    The positioning can be done dynamically via:
        * pin_inputs=Dict(start_pin=Dict(component=f'...',pin='...')) - Specifying start and end
          positions via a component pins

    Pins:
        There are pins called 'L' and 'R' on the ends of the two probes

    Sketch:
        Below is a sketch of the wiring configurations available depending on the position of the pins
        ::
              @@@@@@G
              @@@@  G           # = Taper
              @@    G           @ = Ground plane
             g     #  W         
             g    ##L W         P = Pin to attach taper
             g  ####  W         W = trace_width
              ###  2            G = trace_gap
              1    2            g = orig_gap
            P 1    2            D = dist_extend
              1    2            1 = probe_gap1
              ###  2            2 = probe_gap2
                ####  W         L = New pin to attach a probe wire
                  ##R W         R = New pin to attach a probe wire
                   #  W
              <..D.>
              @@@@
              @@@@@@

        Notice that the pin directions point in direction of wire segments. The turning corner is set at the line intersection of the two pin
        vectors. Thus, this shouldn't be used when the pins are facing away from one another.

    .. image::
        Cap3Interdigital.png

    .. meta::
        Taper to change width of wire at some pin.

    Default Options:
        * trace_width='20um'
        * trace_gap='10um'
        * dist_extend='30um'
        * orig_gap='10um'
        * probe_gap1='2um'
        * probe_gap2='4um'
    """
    default_options = Dict(trace_width='20um',
                           trace_gap='10um',
                           dist_extend='30um',
                           orig_gap='10um',
                           probe_gap1='2um',
                           probe_gap2='4um')

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

        taper1 = [
                (0, width*0.5),
                (0, p.probe_gap1*0.5),
                (p.dist_extend, p.probe_gap2*0.5),
                (p.dist_extend, p.trace_width + p.probe_gap2*0.5)]
        taper2 = [
                (0, -p.probe_gap1*0.5),
                (0, -width*0.5),
                (p.dist_extend, -p.trace_width - p.probe_gap2*0.5),
                (p.dist_extend, -p.probe_gap2*0.5)]

        taper_gap = [
                    (0, width*0.5+gap),
                    (p.dist_extend, p.trace_width + p.probe_gap2*0.5+p.trace_gap),
                    (p.dist_extend, -p.trace_width - p.probe_gap2*0.5-p.trace_gap),
                    (0, -width*0.5-gap)]

        pinL = shapely.LineString(taper1[2:])
        pinR = shapely.LineString(taper2[2:])
        taper1 = shapely.Polygon(taper1)
        taper2 = shapely.Polygon(taper2)
        taper_gap = shapely.Polygon(taper_gap)

        polys = [taper1, taper2, taper_gap, pinL, pinR]
        polys = draw.rotate(polys, rot_angle, origin=(0, 0), use_radians=True)
        polys = draw.translate(polys, *startPt)
        [taper1, taper2, taper_gap, pinL, pinR] = polys
        
        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(taper1=taper1, taper2=taper2),
                           layer=p.layer)

        #subtracts out ground plane on the layer it's on
        self.add_qgeometry('poly',
                           dict(taper_gap=taper_gap),
                           subtract=True,
                           layer=p.layer)

        # Generates its own pins
        self.add_pin('L', pinL.coords[::-1], width=p.trace_width)
        self.add_pin('R', pinR.coords[::-1], width=p.trace_width)

