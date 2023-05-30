from qiskit_metal.qlibrary.core import QComponent, QRoute
from qiskit_metal.qlibrary.tlines.anchored_path import RouteAnchors
import numpy as np
import shapely
from qiskit_metal.toolbox_python.attr_dict import Dict

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
                           trace_width='10um', trace_gap='10um', fillet='50um',
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
        pins = np.vstack(pins)
        
        line = shapely.LineString(pins)

        self.add_qgeometry('path', {'trace': line},
                           width=p.trace_width,
                           fillet=p.fillet,
                           layer=p.layer)
        self.add_qgeometry('path', {'cut': line},
                               width=p.trace_width + 2 * p.trace_gap,
                               fillet=p.fillet,
                               layer=p.layer,
                               subtract=True)
