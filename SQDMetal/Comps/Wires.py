# -*- coding: utf-8 -*-

# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0
# Author: Prasanna Pakkiam
# Creation Date: 01/05/2023
# Description: Collection of classes to wires - from routing objects to tapers.

from qiskit_metal.qlibrary.core import QComponent, QRoute
from qiskit_metal import draw
import itertools
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
        * pin_inputs=Dict(start_pin=Dict(component='...',pin='...'), end_pin=Dict(component='...',pin='...')) - Specifying start and end
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
        p = self.p  # noqa: F841 # abhishekchak52: unused variable p

        start_point = self.design.components[self.options.pin_inputs.start_pin.component].pins[self.options.pin_inputs.start_pin.pin]
        startPt = start_point['middle']
        start_norm = start_point['normal']

        end_point = self.design.components[self.options.pin_inputs.end_pin.component].pins[self.options.pin_inputs.end_pin.pin]
        endPt = end_point['middle']
        end_norm = end_point['normal']
        
        try:
            t,u = np.linalg.solve([[start_norm[0], end_norm[0]], [start_norm[1], end_norm[1]]], [endPt[0]-startPt[0], endPt[1]-startPt[1]])
        except Exception:
            raise Exception("The pins are parallel and this construct should not be used. Try for example, WireElbowParallelPinPin or WirePins instead.")
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
            #In case it is using a parameter - e.g. cpw_width or cpw_gap...
            if isinstance(gap, str):
                if gap in self.design.variables:
                    gap = QUtilities.parse_value_length(self.design.variables[gap])
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
        * layer_probe2 - If greater than -1, it will place the second probe into this layer.
        
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
         w    ###  2            G = trace_gap
         w    1    2            g = orig_gap
         w  P 1    2            D = dist_extend
         w    1    2            1 = probe_gap1
         w    ###  2            2 = probe_gap2
                ####  W         L = New pin to attach a probe wire
                  ##R W         R = New pin to attach a probe wire
                   #  W         w = The width parameter of pin P
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
        * layer_probe2=-1
    """
    default_options = Dict(trace_width='20um',
                           trace_gap='10um',
                           dist_extend='30um',
                           orig_gap='10um',
                           probe_gap1='2um',
                           probe_gap2='4um',
                           layer_probe2=-1)

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
            #In case it is using a parameter - e.g. cpw_width or cpw_gap...
            if isinstance(gap, str):
                gap = QUtilities.parse_value_length(self.design.variables[gap])
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
                           dict(taper1=taper1),
                           layer=p.layer)

        if p.layer_probe2 < 0:
            p.layer_probe2 = p.layer
        self.add_qgeometry('poly',
                           dict(taper2=taper2),
                           layer=p.layer_probe2)

        #subtracts out ground plane on the layer it's on
        self.add_qgeometry('poly',
                           dict(taper_gap=taper_gap),
                           subtract=True,
                           layer=p.layer)

        # Generates its own pins
        self.add_pin('L', pinL.coords[::-1], width=p.trace_width)
        self.add_pin('R', pinR.coords[::-1], width=p.trace_width)

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
    
    Patches:
        * add_patches   - add patches to probe pads

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
                           patches=True,
                           patch_layer=3
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

        # patches
        if p.patches:
            patch_U_L = draw.rectangle(p.pad_width, p.pad_height,  xoff=-0.5*trace_long+0.5*p.pad_width-0.5*p.trace_width,  yoff=0.5*p.T_height)
            patch_U_ML = draw.rectangle(p.pad_width, p.pad_height, xoff=-0.5*trace_short-0.5*p.pad_width+0.5*p.trace_width, yoff=0.5*p.T_height)
            patch_U_MR = draw.rectangle(p.pad_width, p.pad_height, xoff=0.5*trace_short+0.5*p.pad_width-0.5*p.trace_width,  yoff=0.5*p.T_height)
            patch_U_R = draw.rectangle(p.pad_width, p.pad_height,  xoff=0.5*trace_long-0.5*p.pad_width+0.5*p.trace_width,  yoff=0.5*p.T_height)
            patch_L_L = draw.rectangle(p.pad_width, p.pad_height,  xoff=-0.5*trace_long+0.5*p.pad_width-0.5*p.trace_width,  yoff=-0.5*p.T_height)
            patch_L_ML = draw.rectangle(p.pad_width, p.pad_height, xoff=-0.5*trace_short-0.5*p.pad_width+0.5*p.trace_width, yoff=-0.5*p.T_height)
            patch_L_MR = draw.rectangle(p.pad_width, p.pad_height, xoff=0.5*trace_short+0.5*p.pad_width-0.5*p.trace_width,  yoff=-0.5*p.T_height)
            patch_L_R = draw.rectangle(p.pad_width, p.pad_height,  xoff=0.5*trace_long-0.5*p.pad_width+0.5*p.trace_width,  yoff=-0.5*p.T_height)
            patches = shapely.ops.unary_union([patch_L_L, patch_L_ML, patch_L_MR, patch_L_R, patch_U_L, patch_U_ML, patch_U_MR, patch_U_R])
            patches = patches.buffer(p.fillet_radius).buffer(-p.fillet_radius)
            patches = patches.buffer(-0.48*p.trace_width).buffer(0.48*p.trace_width)
        else: 
            patches = []

        # rotate and translate
        polys = [main_trace, pocket, pads_and_Ts, patches]
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        [main_trace, pocket, pads_and_Ts, patches] = polys

        # subtract pocket
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
    
        # add patches
        if patches:
            self.add_qgeometry('poly', 
                           dict(patches=patches), 
                           subtract=False,
                           layer=p.patch_layer)
        
        

class Lollipop(QComponent):
    """
    A circle connected to a straigh wire with a pin (named "pin") at the end. Intended to be used as a charge line.
    Mostly copied from _FluxoniumPocket originally made by Figen YILMAZ and Christian Kraglund Andersen.

    Created by Alex Nguyen
    
               |                               |           ^
               |                               |           |
               |                               |           |
               |              | |   L          |           |----->  X
               |              | |   L          |    
               |              | |   L          |           x = (pox_x, pos_y)
               |              | |   L          |           L = length
               |              |_|   L          |           P = Pin 
               |             /   \             |           ^   
               |            |  x  |            |           |
               |             \   /             |
               |               -               |
               |               P               |
    
    
    """

    default_options=Dict(pos_x='0um', pos_y='0um', orientation=0,
                         length = '100um', cpw_width='cpw_width', cpw_gap='cpw_gap', radius='cpw_width', make_pocket=True)
    
    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p

        pos_x=p.pos_x
        pos_y=p.pos_y
        length=p.length
        cpw_width=p.cpw_width
        cpw_gap=p.cpw_gap
        radius=p.radius
        make_pocket=p.make_pocket


    # Define the geometry
        # Charge Line
        charge_line = draw.rectangle(cpw_width, length, 0, 0)
        charge_line_gap = draw.rectangle(cpw_width+2*cpw_gap, length+cpw_gap/2, 0, -cpw_gap/4)

        # Making the charge_line and charge_line_gap circle and union them to charge_line and it's gap
        charge_line_round = draw.Point(0., (-length)/2).buffer(radius)
        charge_line_gap_round = draw.Point(0., -(length+radius)/2.2).buffer((radius+2*cpw_gap))
        charge_line = draw.union(charge_line, charge_line_round)
        charge_line_gap = draw.union(charge_line_gap, charge_line_gap_round)

        # Charge Line CPW wire
        port_line = draw.LineString([(-cpw_width/2, length/2),
                                     (cpw_width/2, length/2)])
        
        # Position the charge_line, rotate and translate
        objects = [charge_line, charge_line_gap, port_line]
        objects = draw.scale(objects, 1, -1, origin=(0, 0))#Figen Yilmaz made an upside lollipop
        objects = draw.rotate(objects, p.orientation, origin=(0,0))
        objects = draw.translate(
            objects,
            pos_x,
            pos_y)
        [charge_line, charge_line_gap, port_line] = objects

        self.add_qgeometry('poly', {'charge_line': charge_line})
        if make_pocket:
            self.add_qgeometry('poly', {'charge_line_gap': charge_line_gap}, subtract=True)

        ############################################################

        # add pins
        port_line_cords = list(draw.shapely.geometry.shape(port_line).coords)
        port_line_cords = port_line_cords[::-1]# if loc_H==1 else port_line_cords[::-1]
        points = list(port_line_cords)
        self.add_pin('pin', 
                    points, cpw_width)
        
class Wire(QComponent):
    """
    A straight CPW with pins on both ends. 

    Created by Alex Nguyen
                               P1
                               _
               |              | |   L          |           |----->  X
               |              | |   L          |    
               |              | |   L          |           x = (pox_x, pos_y)
               |              | |   L          |           L = length
               |              |_|   L          |           P1 = Pin1
               |               P2              |           P2 = Pin2
    
    
    """

    default_options=Dict(pos_x='0um', pos_y='0um', orientation=0,
                         length = '100um', cpw_width='cpw_width', cpw_gap='cpw_gap')
    
    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p

        pos_x=p.pos_x
        pos_y=p.pos_y
        length=p.length
        cpw_width=p.cpw_width
        cpw_gap=p.cpw_gap


    # Define the geometry
        # Charge Line
        charge_line = draw.rectangle(cpw_width, length, 0, 0)
        charge_line_gap = draw.rectangle(cpw_width+2*cpw_gap, length, 0, 0)

        

        # For the pins
        port_line1 = draw.LineString([(-cpw_width/2, length/2),
                                     (cpw_width/2, length/2)])
        
        port_line2 = draw.LineString([(-cpw_width/2, -length/2),
                                     (cpw_width/2, -length/2)])
        
        # Position the charge_line, rotate and translate
        objects = [charge_line, charge_line_gap, port_line1, port_line2]
        objects = draw.rotate(objects, p.orientation, origin=(0, 0))
        objects = draw.translate(
            objects,
            pos_x,
            pos_y)
        [charge_line, charge_line_gap, port_line1, port_line2] = objects

        self.add_qgeometry('poly', {'charge_line': charge_line})
        self.add_qgeometry('poly', {'charge_line_gap': charge_line_gap}, subtract=True)

        ############################################################

        # add pins
        port_line_cords1 = list(draw.shapely.geometry.shape(port_line1).coords)
        port_line_cords1 = port_line_cords1# if loc_H==1 else port_line_cords[::-1]
        points1 = list(port_line_cords1)
        self.add_pin('pin1', 
                    points1, cpw_width)
        
        port_line_cords2 = list(draw.shapely.geometry.shape(port_line2).coords)
        port_line_cords2 = port_line_cords2[::-1]# if loc_H==1 else port_line_cords[::-1]
        points2 = list(port_line_cords2)
        self.add_pin('pin2', 
                    points2, cpw_width)


class WireMeander(QRoute):
    """Meandered CPW between two pins that solves for an exact (filleted) length and searches for a
    layout that fits. Similar in use to Qiskit-Metal's RouteMeander, but more robust:

        * Fillets are computed by this component (QUtilities.calc_lines_and_fillets_on_path) and
          drawn into the geometry, so the renderer never re-fillets it and the qiskit-metal
          "short segments ... fillet" warning cannot occur. Where corners are too close for the
          full radius, the radius is reduced locally instead of being left sharp.
        * total_length is the true length of the drawn (filleted) centre line.
        * The meander amplitude is solved numerically for the requested length.
        * If the first layout breaks a rule (fillet radius reduced below fillet_min, the path
          crossing itself, or overlapping a component in avoid_components), it tries other
          layouts: mirrored, biased to one side, fewer legs, tighter spacing. The first layout
          that passes is drawn.
        * If no layout can reach total_length (pins too far apart, or not enough room under
          max_amplitude), the build fails with a message giving the reachable length range.
        * Optional anchors (waypoints) the route passes through.
        * Optional bonded crossings over other CPWs (e.g. a feedline): the trace is broken either
          side of the crossed line with open ends (optionally bond pads) and a WireBond is placed
          across.

    Inherits QRoute class.

    Layout:
        A lead of length start_straight leaves the start pin and one of end_straight leaves the end pin
        (along each pin's normal). Between the two lead tips (S and E), the meander has legs
        perpendicular to the line S->E, stacked along it every meander.spacing:

                     ______        ______
                    |      |      |      |                spacing  = centre-to-centre leg distance
            S ___   |      |      |      |   ___ E        amplitude = solved for total_length
                 |  |      |      |      |  |
                 |__|      |______|      |__|

    The CPW is given via:
        * trace_width, trace_gap - CPW dimensions
        * fillet            - Requested corner radius (use <= spacing/2 for semicircular U-turns)
        * fillet_min        - Smallest acceptable radius after local reduction ('' -> 0.99*fillet)
        * fillet_resolution - Points per fillet arc in the drawn geometry

    Length and leads:
        * total_length - Length of the filleted centre line from start pin to end pin. Leave
                         empty ('') for the shortest route instead (no meander): with snap, the
                         first of an L or Z shape along x or y that passes the checks. Useful
                         for control lines that only need to get from A to B.
        * lead.start_straight, lead.end_straight - Straight sections off each pin (use >= fillet)
        * end_gap - If set, the end pin is an open circuit: the ground cut continues this far past
                    the end of the trace (e.g. a charge line ending near a qubit). '' -> no gap.

    Anchors:
        * anchors - List of (x, y) points the route passes through, in order (e.g.
                    [('-3mm', '1.2mm'), ('-1mm', '1.2mm')]). The route is split into legs
                    S -> anchor 1 -> ... -> E. Every leg takes the shortest shape (straight, L or
                    Z) that passes the checks; with total_length, one leg also carries the meander.
        * meander.leg - Index of the leg that carries the meander ('' -> the longest leg)

    Meander:
        * meander.spacing       - Distance between adjacent legs (centre-to-centre). Should be at
                                  least trace_width + 2*trace_gap (+ some ground in between).
        * meander.asymmetry     - Offset of the meander centre line from the line S->E
        * meander.max_amplitude - Maximum leg half-height ('' -> unlimited)
        * meander.search        - If True, also tries mirrored/one-sided layouts, fewer legs and
                                  tighter spacing (x0.75, x0.5; never below trace_width +
                                  2*trace_gap + clearance or 2*fillet_min)

    Snapping:
        * snap - If True (default), the meander runs along whichever of x or y has the longest
                 S->E projection (most room for legs) and a straight jog joins it to E, keeping
                 the route Manhattan. If False, it runs straight along S->E.

    Obstacles:
        * avoid_components - List of component names to keep clear of, or 'all' (every other
                             component); 'all' can be mixed with explicit entries, which
                             override it. An entry may also be (name, 'metal') to avoid only that
                             component's metal (e.g. cross a qubit pocket to reach a flux inlet
                             while keeping clear of the pads). The straight leads are exempt from
                             the two components the route connects to (so those can be listed
                             too); make the leads long enough to clear them.
        * clearance        - Minimum ground between this CPW's gap and any avoided geometry

    Bonded crossings:
        * cross_components - Names of CPW components this route may cross (e.g. a feedline).
                             They are not obstacles; instead, each crossing must be at 90 degrees
                             on a straight stretch. Crossing anything else still counts as
                             "too close".
        * crossing.bond_gap   - Ground left between the crossed CPW's gap and this route's open
                                end (the open end itself has a trace_gap cut beyond it)
        * crossing.pad_width, crossing.pad_length - Bond pad at each open end ('' or 0 -> no
                                pad: the plain trace end is the bond landing)
        * crossing.add_wirebonds - If True, places a WireBond (SQDMetal.Comps.WireBonds) across
                                each crossing, named <this name>_bond<k>
        * crossing.wirebond   - Extra options passed to each WireBond (e.g. width, material)
        At each crossing the trace is broken: going outwards from the crossed CPW's centre there
        is its own gap, then bond_gap of ground, then this route's end gap (trace_gap), then the
        open end (and pad). The bond wire is not part of the drawn length, so a route with
        both total_length and crossings logs a warning.

    After building, the chosen layout is stored in self.route_info (legs, amplitude, bias, mirror,
    length, smallest fillet radius, issues, layouts tried, shape, crossings).

    Pins:
        'start' and 'end' (as with all QRoutes).

    Default Options:
        * trace_width='10um', trace_gap='10um'
        * fillet='30um', fillet_min='', fillet_resolution=16
        * total_length='5mm'
        * lead=Dict(start_straight='100um', end_straight='50um')
        * anchors=[]
        * meander=Dict(spacing='200um', asymmetry='0um', max_amplitude='', search=True, leg='')
        * avoid_components='all', clearance='10um'
        * cross_components=[]
        * crossing=Dict(bond_gap='20um', pad_width='60um', pad_length='80um', add_wirebonds=True,
                        wirebond=Dict())
        * end_gap=''
        * snap=True
    """

    default_options = Dict(trace_width='10um',
                           trace_gap='10um',
                           fillet='30um',
                           fillet_min='',
                           fillet_resolution=16,
                           total_length='5mm',
                           lead=Dict(start_straight='100um', end_straight='50um'),
                           anchors=[],
                           meander=Dict(spacing='200um', asymmetry='0um', max_amplitude='', search=True,
                                        leg=''),
                           avoid_components='all',
                           clearance='10um',
                           cross_components=[],
                           crossing=Dict(bond_gap='20um', pad_width='60um', pad_length='80um',
                                         add_wirebonds=True, wirebond=Dict()),
                           end_gap='',
                           snap=True)

    # Layouts to try, in order: (bias, mirror). bias shifts legs to one side of the S->E line
    # (0 = centred, +-0.8 = mostly on one side); mirror flips which side the first leg goes to.
    _SEARCH_LAYOUTS = [(0.0, 1), (0.0, -1), (0.5, 1), (0.5, -1), (-0.5, 1), (-0.5, -1),
                       (0.8, 1), (0.8, -1), (-0.8, 1), (-0.8, -1)]
    # Most leg-shape combinations evaluated by the shortest-route search (with many anchors, only
    # the first two shapes of each leg are combined)
    _MAX_LEG_COMBOS = 256

    @staticmethod
    def _optional_length(val):
        """Parsed optional length -> positive float, or None when unset. Anything that is not a
        positive number counts as unset: '', None, or an empty Dict (which is what a qiskit Dict
        returns for a missing key)."""
        try:
            val = float(val)
        except (TypeError, ValueError):
            return None
        return val if val > 0 else None

    @staticmethod
    def _is_true(val):
        return str(val).strip().lower() in ('true', '1')

    @staticmethod
    def _snap_axis(delta, start_dir):
        """x/y unit vector to meander along: the axis with the largest positive projection of
        delta (most room for legs). Near-ties (within 5%) go to the start lead's direction."""
        axes = [np.array(a, dtype=float) for a in ((1, 0), (-1, 0), (0, 1), (0, -1))]
        proj = [float(np.dot(delta, a)) for a in axes]
        for a, pr in zip(axes, proj):
            if np.dot(a, start_dir) > 1 - 1e-9 and pr >= 0.95 * max(proj):
                return a
        return axes[int(np.argmax(proj))]

    # ------------------------------------------------------------------ #
    def make(self):
        p = self.p
        start = self.set_pin("start")
        end = self.set_pin("end")
        precision = self.design.template_options.PRECISION
        tol = 10 ** (-precision)

        r = p.fillet
        r_min = self._optional_length(p.fillet_min)
        if r_min is None:
            r_min = 0.99 * r
        max_amp = self._optional_length(p.meander.max_amplitude)
        target = self._optional_length(p.total_length)   # None -> shortest route

        # Lead tips S, E and the anchors between them: the route's legs run node to node
        S = start.position + start.direction * p.lead.start_straight
        E = end.position + end.direction * p.lead.end_straight
        anchors = [np.array([float(a[0]), float(a[1])]) for a in (p.anchors or [])]
        nodes = [S, *anchors, E]
        for k in range(len(nodes) - 1):
            if np.linalg.norm(nodes[k + 1] - nodes[k]) < tol:
                raise ValueError(f"{self.name}: consecutive route points {k} and {k + 1} (lead tips / "
                                 f"anchors) coincide; move an anchor or change the lead lengths.")

        obstacles = self._get_obstacles(p)

        # Shortest route through all the legs (this is also the base for the meander)
        short, tried = self._shortest_route(p, start, end, nodes, r, r_min, precision, obstacles)
        if target is None:
            self._finish(p, short, tried, r, precision, obstacles, target)
            return

        # --- Meander on one leg; the other legs keep their shortest-route shapes ---
        legs = short.legs
        n_route_legs = len(legs)
        m = self._optional_int(p.meander.leg)
        if m is None:
            m = int(np.argmax([np.linalg.norm(nodes[k + 1] - nodes[k]) for k in range(n_route_legs)]))
        if not 0 <= m < n_route_legs:
            raise ValueError(f"{self.name}: meander.leg={m} but the route only has legs 0-{n_route_legs - 1}.")
        P, Q = nodes[m], nodes[m + 1]
        prefix = [start.position] + [pt for leg in legs[:m] for pt in leg[:-1]]
        suffix = [pt for leg in legs[m + 1:] for pt in leg[1:]] + [end.position]
        dir_in = start.direction if m == 0 else self._unit(legs[m - 1][-1] - legs[m - 1][-2])
        dir_out = -end.direction if m == n_route_legs - 1 else self._unit(legs[m + 1][1] - legs[m + 1][0])

        # Meander frame: u is the meander axis, v is perpendicular to it. With snap, u is the x/y
        # axis with the longest P->Q projection. In this frame P=(0, 0) and Q=(D, h): the
        # meander's last leg ends at height h, so no extra jog is needed.
        delta = Q - P
        if self._is_true(p.snap):
            u = self._snap_axis(delta, dir_in)
        else:
            u = delta / np.linalg.norm(delta)
        v = np.array([-u[1], u[0]])
        D = float(np.dot(delta, u))
        h = float(np.dot(delta, v))
        if abs(h) < tol:
            h = 0.0

        # Leg margins at each end of the axis: room for 1 fillet, or 2 if there is a corner there too
        corner_S = np.dot(dir_in, u) < 1 - 1e-9
        corner_E = np.dot(dir_out, u) < 1 - 1e-9
        margin_min = r * (1 + max(corner_S, corner_E))

        # A U-turn between legs can be at most a semicircle of radius spacing/2. If the requested
        # fillet is bigger, accept that reduction once here rather than failing every layout on it.
        base_spacing = p.meander.spacing
        if r_min > base_spacing / 2:
            self.logger.warning(f"{self.name}: fillet ({r:.4g}) > spacing/2 ({base_spacing / 2:.4g}); "
                                f"U-turns will use radius spacing/2. Reduce fillet or increase spacing.")
            r_min = 0.99 * base_spacing / 2

        # Spacings to try: the requested one, then tighter ones (search only) that still leave
        # ground between adjacent legs and room for full-radius U-turns
        min_spacing = max(p.trace_width + 2 * p.trace_gap + p.clearance, 2 * r_min)
        spacings = [base_spacing]
        if self._is_true(p.meander.search):
            spacings += [base_spacing * s for s in (0.75, 0.5) if base_spacing * s >= min_spacing]
        layouts = self._SEARCH_LAYOUTS if self._is_true(p.meander.search) else self._SEARCH_LAYOUTS[:1]

        def assemble(local):
            local = np.asarray(local, dtype=float)
            world = P + np.outer(local[:, 0], u) + np.outer(local[:, 1], v)
            return np.vstack([*prefix, world, *suffix])

        def build(n_legs, amp, bias, mirror, spacing):
            # Built mirrored, then flipped back, so the end height h is the same in both cases
            local = self._meander_local_pts(D, n_legs, spacing, amp, bias, p.meander.asymmetry,
                                            mirror * h)
            local[:, 1] *= mirror
            return assemble(local)

        def length_of(pts):
            segs = QUtilities.calc_lines_and_fillets_on_path(pts, r, precision)
            return sum(s['dist'] for s in segs), segs

        best, solved_any = None, False
        reach = [short.length, short.length]    # [shortest, longest] lengths seen, for errors

        for spacing in spacings:
            n_max = int(np.floor((D - 2 * margin_min) / spacing)) + 1 if D > 2 * margin_min else 0
            if n_max < 2 and spacing == base_spacing and D - 2 * margin_min >= 2 * r:
                spacing = D - 2 * margin_min    # squeeze in a single U-turn
                n_max = 2
            # Most legs first (most compact), trying every layout at each leg count, so a compact
            # one-sided layout wins over a tall loop with few legs
            too_short = set()   # layouts that cannot reach the length even at max amplitude
            for n_legs in range(n_max, 1, -1):
                for bias, mirror in layouts:
                    if (bias, mirror) in too_short:
                        continue    # fewer legs would only be shorter still
                    amp_lo = r
                    amp_hi = max_amp if max_amp is not None else max(target, r)
                    len_lo = length_of(build(n_legs, amp_lo, bias, mirror, spacing))[0]
                    len_hi = length_of(build(n_legs, amp_hi, bias, mirror, spacing))[0]
                    reach = [min(reach[0], len_lo), max(reach[1], len_hi)]
                    if len_lo > target:
                        continue    # too long even at minimum amplitude: needs fewer legs
                    if len_hi < target:
                        too_short.add((bias, mirror))
                        continue
                    amp = self._solve_amplitude(
                        lambda a: length_of(build(n_legs, a, bias, mirror, spacing))[0],
                        target, amp_lo, amp_hi, len_lo, len_hi)
                    pts = build(n_legs, amp, bias, mirror, spacing)
                    length, segs = length_of(pts)
                    issues = self._check_layout(p, pts, segs, r, r_min, spacing, obstacles, precision)
                    tried += 1
                    solved_any = True
                    info = Dict(n_legs=n_legs, amplitude=amp, bias=bias, mirror=mirror, spacing=spacing,
                                length=length, issues=issues, segs=segs, pts=pts, meander_leg=m)
                    if best is None or len(issues) < len(best.issues):
                        best = info
                    if not issues:
                        break
                if best is not None and not best.issues:
                    break
            if best is not None and not best.issues:
                break

        # The shortest route is fine if it already has the right length
        if (best is None or best.issues) and abs(short.length - target) < 1e-6:
            best, solved_any = short, True

        if not solved_any:
            units = self.design.get_units()
            raise ValueError(
                f"{self.name}: total_length={target:.4f}{units} cannot be reached. Reachable lengths "
                f"with these leads/anchors/spacing{'/max_amplitude' if max_amp is not None else ''} are "
                f"{reach[0]:.4f}-{reach[1]:.4f}{units} (meander leg {m} spans "
                f"{np.linalg.norm(delta):.4f}{units}). Move the pins or anchors, shorten the leads, or "
                f"change spacing/max_amplitude.")

        self._finish(p, best, tried, r, precision, obstacles, target)

    # ------------------------------------------------------------------ #
    @staticmethod
    def _unit(vec):
        return vec / np.linalg.norm(vec)

    @staticmethod
    def _optional_int(val):
        try:
            return int(val)
        except (TypeError, ValueError):
            return None

    def _leg_shapes(self, P, Q, dir_in, snap, tol):
        """Candidate corner lists (P ... Q inclusive) for one leg: straight if P and Q line up;
        otherwise, with snap, an L with its corner at either end and a Z through the middle,
        along both axes."""
        delta = Q - P
        if not snap:
            return [("straight", np.array([P, Q]))]
        first = self._snap_axis(delta, dir_in)
        other = np.array([first[1], -first[0]])
        if np.dot(delta, other) < 0:
            other = -other
        frames = [first] + ([other] if np.dot(delta, other) > tol else [])
        shapes = []
        for u in frames:
            v = np.array([-u[1], u[0]])
            D, h = float(np.dot(delta, u)), float(np.dot(delta, v))
            if abs(h) <= tol:
                local = [("straight", [(0.0, 0.0), (D, 0.0)])]
            else:
                local = [("L (corner at end)", [(0.0, 0.0), (D, 0.0), (D, h)]),
                         ("L (corner at start)", [(0.0, 0.0), (0.0, h), (D, h)]),
                         ("Z", [(0.0, 0.0), (D / 2, 0.0), (D / 2, h), (D, h)])]
            for name, pts in local:
                pts = np.asarray(pts, dtype=float)
                world = P + np.outer(pts[:, 0], u) + np.outer(pts[:, 1], v)
                world[-1] = Q       # exact end point (no rounding drift)
                shapes.append((name, world))
        return shapes

    def _shortest_route(self, p, start, end, nodes, r, r_min, precision, obstacles):
        """The shortest route through the nodes (lead tips and anchors): combinations of each
        leg's shapes, keeping the first that passes the checks (otherwise the one with the fewest
        issues, then the shortest)."""
        tol = 10 ** (-precision)
        snap = self._is_true(p.snap)
        options = []
        for k in range(len(nodes) - 1):
            dir_in = start.direction if k == 0 else self._unit(nodes[k] - nodes[k - 1])
            options.append(self._leg_shapes(nodes[k], nodes[k + 1], dir_in, snap, tol))
        n_combos = int(np.prod([len(o) for o in options]))
        if n_combos > self._MAX_LEG_COMBOS:
            options = [o[:2] for o in options]

        best, best_key, tried = None, None, 0
        for combo in itertools.product(*options):
            legs = [shape[1] for shape in combo]
            pts = np.vstack([start.position] + [leg if k == 0 else leg[1:] for k, leg in enumerate(legs)]
                            + [end.position])
            segs = QUtilities.calc_lines_and_fillets_on_path(pts, r, precision)
            length = sum(s['dist'] for s in segs)
            issues = self._check_layout(p, pts, segs, r, r_min, None, obstacles, precision)
            tried += 1
            key = (len(issues), length)
            if best is None or key < best_key:
                best_key = key
                best = Dict(n_legs=0, amplitude=0.0, bias=0.0, mirror=1, spacing=None, length=length,
                            issues=issues, segs=segs, pts=pts, legs=legs,
                            shape=" + ".join(shape[0] for shape in combo))
            if not issues:
                break
        return best, tried

    # ------------------------------------------------------------------ #
    def _finish(self, p, best, tried, r, precision, obstacles, target):
        """Draws the chosen layout (with any bonded crossings) and records route_info."""
        if best.issues:
            self.logger.warning(f"{self.name}: no layout passed all checks after {tried} "
                                f"attempt(s); drawing the best one. Issues: {'; '.join(best.issues)}")

        # Filleted centre line (no renderer filleting: fillet=0 in the geometry tables)
        fpts = QUtilities.calc_filleted_path(best.pts, r, precision, resolution=p.fillet_resolution)
        keep = np.r_[True, np.linalg.norm(np.diff(fpts, axis=0), axis=1) > 10 ** (-precision)]
        line = shapely.LineString(fpts[keep])

        crossings = self._find_crossings(p, line, obstacles[2])
        self._remove_old_bonds()
        if not crossings:
            self.add_qgeometry('path', {'trace': line}, width=p.trace_width, fillet=0, layer=p.layer)
            self.add_qgeometry('path', {'cut': line}, width=p.trace_width + 2 * p.trace_gap, fillet=0,
                               layer=p.layer, subtract=True)
        else:
            self._draw_with_crossings(p, line, crossings)
            if target is not None:
                self.logger.warning(f"{self.name}: total_length includes {len(crossings)} bonded "
                                    f"crossing(s); the bond wires are not part of the drawn length.")

        # Open-circuit end: the ground cut continues end_gap past the end of the trace
        end_gap = self._optional_length(p.end_gap)
        if end_gap is not None:
            coords = np.array(line.coords)
            tip = coords[-1]
            direction = self._unit(tip - coords[-2])
            self.add_qgeometry('path', {'cut_open_end': shapely.LineString([tip, tip + direction * end_gap])},
                               width=p.trace_width + 2 * p.trace_gap, fillet=0, layer=p.layer,
                               subtract=True)

        radii = [s['radius'] for s in best.segs if 'radius' in s]
        self._filleted_length = best.length
        self.options._actual_length = f"{best.length} {self.design.get_units()}"
        self.intermediate_pts = best.pts[2:-2]      # keeps QRoute.get_points() meaningful
        self.route_info = Dict(n_legs=best.n_legs, amplitude=best.amplitude, bias=best.bias,
                               mirror=best.mirror, spacing=best.spacing, length=best.length,
                               min_fillet=min(radii) if radii else None, issues=best.issues,
                               layouts_tried=tried, shape=best.get('shape', 'meander'),
                               meander_leg=best.get('meander_leg'),
                               crossings=[Dict(crossed=c.name, position=tuple(c.pos), bond=c.get('bond'))
                                          for c in crossings])

    # ------------------------------------------------------------------ #
    # Bonded crossings
    def _crossing_dims(self, p, cut_half):
        """Distances (from the crossed CPW's centre line, along this route) of this route's open
        trace end (t_end) and how far the straight stretch must reach on each side (reach)."""
        pad_w = self._optional_length(p.crossing.pad_width)
        pad_l = self._optional_length(p.crossing.pad_length)
        has_pad = pad_w is not None and pad_l is not None
        t_end = cut_half + p.crossing.bond_gap + p.trace_gap
        reach = t_end + (pad_l + p.trace_gap if has_pad else p.trace_gap)
        return t_end, reach, (pad_w, pad_l) if has_pad else None

    def _find_crossings(self, p, centre, crossers):
        """Where the centre line crosses the traces of the cross_components: list of Dicts with
        the crossing position, arclength s, this route's direction u, the crossed component and
        its half cut width, plus 'problems' (not at 90 degrees, on a bend, overlapping)."""
        from shapely.ops import substring
        found = []
        ds = min(1e-4, centre.length / 1e3)

        def direction(line, s):
            a = np.array(line.interpolate(max(s - ds, 0)).coords[0])
            b = np.array(line.interpolate(min(s + ds, line.length)).coords[0])
            return (b - a) / np.linalg.norm(b - a)

        for name, traces, cut_half in crossers:
            _, reach, _ = self._crossing_dims(p, cut_half)
            for trace in traces:
                inter = centre.intersection(trace)
                for g in [g for g in getattr(inter, 'geoms', [inter]) if not g.is_empty]:
                    pt = g if g.geom_type == 'Point' else g.representative_point()
                    s = centre.project(pt)
                    u, w = direction(centre, s), direction(trace, trace.project(pt))
                    problems = []
                    if g.geom_type != 'Point':
                        problems.append(f"runs along {name} instead of crossing it")
                    elif abs(np.dot(u, w)) > 1e-3:
                        problems.append(f"crosses {name} at {np.degrees(np.arccos(abs(np.dot(u, w)))):.1f} "
                                        f"degrees, not 90")
                    s0, s1 = s - reach, s + reach
                    if s0 < 0 or s1 > centre.length:
                        problems.append(f"crossing of {name} too close to a route end")
                    else:
                        piece = substring(centre, s0, s1)
                        chord = np.linalg.norm(np.array(piece.coords[-1]) - np.array(piece.coords[0]))
                        if piece.length - chord > 1e-6:
                            problems.append(f"crossing of {name} is on or too near a bend")
                    found.append(Dict(name=name, pos=np.array(pt.coords[0]), s=s, u=u,
                                      cut_half=cut_half, problems=problems))
        found.sort(key=lambda c: c.s)
        for a, b in zip(found, found[1:]):
            if b.s - a.s < self._crossing_dims(p, a.cut_half)[1] + self._crossing_dims(p, b.cut_half)[1]:
                b.problems.append("crossings too close together")
        return found

    def _crossing_pads(self, p, crossing):
        """(metal, cut) pad polygons for both sides of a crossing (empty if pads are off)."""
        t_end, _, pad = self._crossing_dims(p, crossing.cut_half)
        if pad is None:
            return [], []
        pad_w, pad_l = pad
        u = crossing.u
        v = np.array([-u[1], u[0]])
        metal, cut = [], []
        for sgn in (-1, 1):
            def rect(t0, t1, half_w):
                c0, c1 = crossing.pos + sgn * u * t0, crossing.pos + sgn * u * t1
                return shapely.Polygon([c0 + v * half_w, c1 + v * half_w, c1 - v * half_w, c0 - v * half_w])
            metal.append(rect(t_end, t_end + pad_l, pad_w / 2))
            cut.append(rect(t_end - p.trace_gap, t_end + pad_l + p.trace_gap, pad_w / 2 + p.trace_gap))
        return metal, cut

    def _draw_with_crossings(self, p, line, crossings):
        """Trace broken at each crossing (open ends, optional pads) plus a WireBond across each."""
        from shapely.ops import substring
        L = line.length
        trace_breaks, cut_breaks = [], []
        for c in crossings:
            t_end, _, _ = self._crossing_dims(p, c.cut_half)
            trace_breaks.append((c.s - t_end, c.s + t_end))
            cut_breaks.append((c.s - (t_end - p.trace_gap), c.s + (t_end - p.trace_gap)))

        def pieces(breaks):
            out, s = [], 0.0
            for b0, b1 in breaks:
                if b0 > s:
                    out.append(substring(line, s, b0))
                s = max(s, b1)
            if s < L:
                out.append(substring(line, s, L))
            return [piece for piece in out if piece.length > 0]

        for k, piece in enumerate(pieces(trace_breaks)):
            self.add_qgeometry('path', {'trace' if k == 0 else f'trace_{k}': piece}, width=p.trace_width,
                               fillet=0, layer=p.layer)
        for k, piece in enumerate(pieces(cut_breaks)):
            self.add_qgeometry('path', {'cut' if k == 0 else f'cut_{k}': piece},
                               width=p.trace_width + 2 * p.trace_gap, fillet=0, layer=p.layer,
                               subtract=True)

        from SQDMetal.Comps.WireBonds import WireBond
        self._bond_names = []
        for k, c in enumerate(crossings):
            metal, cut = self._crossing_pads(p, c)
            for j, (pm, pc) in enumerate(zip(metal, cut)):
                self.add_qgeometry('poly', {f'bond_pad_{k}_{j}': pm}, layer=p.layer)
                self.add_qgeometry('poly', {f'bond_pad_gap_{k}_{j}': pc}, layer=p.layer, subtract=True)
            if self._is_true(p.crossing.add_wirebonds):
                t_end, _, pad = self._crossing_dims(p, c.cut_half)
                t_foot = t_end + (pad[1] / 2 if pad is not None else p.trace_width)
                a, b = c.pos - c.u * t_foot, c.pos + c.u * t_foot
                name = f"{self.name}_bond{k}"
                WireBond(self.design, name, options=Dict(start_x=float(a[0]), start_y=float(a[1]),
                                                         end_x=float(b[0]), end_y=float(b[1]),
                                                         **dict(self.options.crossing.wirebond or {})))
                self._bond_names.append(name)
                c.bond = name

    def _own_bond_names(self):
        """WireBonds made by this route in earlier builds (named <this name>_bond<k>), including
        ones left by a previous component of the same name (e.g. after a parent rebuild)."""
        prefix = f"{self.name}_bond"
        return [n for n in self.design.components.keys()
                if n.startswith(prefix) and n[len(prefix):].isdigit()]

    def _remove_old_bonds(self):
        """Deletes the WireBonds this route made in a previous build."""
        for name in self._own_bond_names():
            self.design.delete_component(name, force=True)
        self._bond_names = []

    # ------------------------------------------------------------------ #
    @property
    def length(self) -> float:
        """True length of the drawn (filleted) centre line."""
        return getattr(self, '_filleted_length', 0.0)

    @staticmethod
    def _meander_local_pts(D, n_legs, spacing, amp, bias, asymmetry, end_level=0.0):
        """Corner points in the (u, v) frame, from S=(0,0) to E=(D,end_level), with n_legs >= 2
        legs. Legs run between v=top and v=bot; the first leg starts on v=0 and the last leg
        ends on v=end_level."""
        top = asymmetry + amp * (1 + bias)
        bot = asymmetry - amp * (1 - bias)
        margin = (D - (n_legs - 1) * spacing) / 2
        xs = margin + spacing * np.arange(n_legs)

        pts = [(0.0, 0.0), (xs[0], 0.0), (xs[0], top)]
        level = top
        for k in range(1, n_legs - 1):
            pts.append((xs[k], level))
            level = bot if level == top else top
            pts.append((xs[k], level))
        pts += [(xs[-1], level), (xs[-1], end_level), (D, end_level)]
        return np.array(pts, dtype=float)

    @staticmethod
    def _solve_amplitude(length_fn, target, lo, hi, f_lo, f_hi, tol=1e-10, max_iter=100):
        """Solves length_fn(amp) == target on [lo, hi] (length grows with amplitude) using the
        Illinois false-position method. Length is piecewise linear in amplitude, so this usually
        converges in a few evaluations; the bracket guarantees convergence regardless."""
        g_lo, g_hi = f_lo - target, f_hi - target
        side = 0
        amp = lo
        for _ in range(max_iter):
            amp = (lo * g_hi - hi * g_lo) / (g_hi - g_lo) if g_hi != g_lo else 0.5 * (lo + hi)
            g = length_fn(amp) - target
            if abs(g) < tol or hi - lo < tol:
                break
            if g < 0:
                lo, g_lo = amp, g
                if side == -1:
                    g_hi *= 0.5
                side = -1
            else:
                hi, g_hi = amp, g
                if side == 1:
                    g_lo *= 0.5
                side = 1
        return amp

    # ------------------------------------------------------------------ #
    # Obstacles and checks
    def _get_obstacles(self, p):
        """Returns (ends, others, crossers):
            * ends   - prepared footprint of the two components this route connects to (or None)
            * others - prepared footprint of every other component to avoid (or None)
            * crossers - [(name, [trace LineStrings], half cut width)] for cross_components
        An avoid_components entry can be a name (all its geometry, including ground cut-outs) or
        (name, 'metal') (only its metal, e.g. to cross a qubit's pocket while keeping clear of its
        pads). Components in cross_components are never obstacles."""
        cross_names = [c for c in (p.cross_components or [])]
        entries = p.avoid_components
        if isinstance(entries, str):
            entries = [entries]
        # name -> metal_only. 'all' adds every other component; explicit entries (e.g.
        # (name, 'metal')) override it, so ['all', ('Q1', 'metal')] avoids everything but only
        # Q1's metal.
        avoid = {}
        if any(isinstance(e, str) and e.strip().lower() == 'all' for e in entries):
            avoid = {n: False for n in self.design.components.keys() if n != self.name}
        for entry in entries:
            if isinstance(entry, (list, tuple)):
                avoid[entry[0]] = str(entry[1]).lower() == 'metal'
            elif not (isinstance(entry, str) and entry.strip().lower() == 'all'):
                avoid[entry] = False
        bond_names = set(self._own_bond_names())
        end_comps = {self.options.pin_inputs.start_pin.component,
                     self.options.pin_inputs.end_pin.component}
        shapes = {True: [], False: []}      # keyed by "is an end component"
        for name, metal_only in avoid.items():
            if name in cross_names or name in bond_names:
                continue
            if name not in self.design.components:
                self.logger.warning(f"{self.name}: avoid_components entry '{name}' is not in the design.")
                continue
            for kind in ('path', 'poly'):
                df = self.design.qgeometry.get_component(name, kind)
                if df is None:
                    continue
                for _, row in df.iterrows():
                    # (row['subtract'], not row.subtract: the latter is pandas' Series.subtract)
                    if metal_only and bool(row['subtract']):
                        continue
                    geom = row['geometry']
                    if kind == 'path':
                        geom = geom.buffer(row['width'] / 2, cap_style=2)
                    shapes[name in end_comps].append(geom)

        crossers = []
        for name in cross_names:
            if name not in self.design.components:
                self.logger.warning(f"{self.name}: cross_components entry '{name}' is not in the design.")
                continue
            df = self.design.qgeometry.get_component(name, 'path')
            if df is None or df.empty:
                self.logger.warning(f"{self.name}: cross_components entry '{name}' has no CPW paths.")
                continue
            traces = [row['geometry'] for _, row in df.iterrows() if not bool(row['subtract'])]
            widths = [row['width'] for _, row in df.iterrows()]
            crossers.append((name, traces, max(widths) / 2))

        def union(geoms):
            if not geoms:
                return None
            u = shapely.unary_union(geoms)
            shapely.prepare(u)
            return u
        return union(shapes[True]), union(shapes[False]), crossers

    def _check_layout(self, p, pts, segs, r, r_min, spacing, obstacles, precision):
        """Human-readable list of problems with a candidate layout (empty = OK)."""
        issues = []
        radii = [s['radius'] for s in segs if 'radius' in s]
        if radii and min(radii) < r_min - 10 ** (-precision):
            issues.append(f"fillet reduced to {min(radii):.4g} (< fillet_min {r_min:.4g})")
        if spacing is not None and spacing < p.trace_width + 2 * p.trace_gap:
            issues.append(f"spacing {spacing:.4g} < trace_width + 2*trace_gap: adjacent gaps merge")

        fpts = QUtilities.calc_filleted_path(pts, r, precision, resolution=4)
        # Drop repeated points (e.g. a fillet that starts exactly on a pin): shapely counts them
        # as self-intersections
        fpts = fpts[np.r_[True, np.linalg.norm(np.diff(fpts, axis=0), axis=1) > 10 ** (-precision)]]
        centre = shapely.LineString(fpts)
        if not centre.is_simple:
            issues.append("path crosses itself")

        ends, others, crossers = obstacles
        crossings = self._find_crossings(p, centre, crossers) if crossers else []
        for c in crossings:
            issues += c.problems
        if ends is None and others is None:
            return issues
        from shapely.ops import substring
        half = p.trace_width / 2 + p.trace_gap + p.clearance

        def too_close(geom, skip_start, skip_end):
            # Footprint of the centre line between the two skipped lengths (1nm less than the full
            # width, so that merely touching an avoided shape does not count)
            s0, s1 = skip_start, centre.length - skip_end
            if geom is None or s1 <= s0:
                return False
            return geom.intersects(substring(centre, s0, s1).buffer(half - 1e-6, cap_style=2))

        # Other components: the whole route counts, except the footprint half-width at each pin
        if too_close(others, half, half):
            issues.append("too close to an avoided component")
        # The components at either end: the straight leads (minus the fillet that starts on them)
        # are exempt, since they have to touch them
        elif too_close(ends, max(p.lead.start_straight - r, half), max(p.lead.end_straight - r, half)):
            issues.append("too close to the component at a route end")
        # Bond pads (wider than the trace) must keep clear too
        if others is not None:
            for c in crossings:
                _, pad_cuts = self._crossing_pads(p, c)
                if any(others.intersects(pc.buffer(p.clearance - 1e-6)) for pc in pad_cuts):
                    issues.append(f"bond pad at the {c.name} crossing too close to an avoided component")
                    break
        return issues
