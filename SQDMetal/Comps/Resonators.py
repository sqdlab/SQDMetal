# -*- coding: utf-8 -*-
# Author: Prasanna Pakkiam
# Creation Date: 01/06/2023
# Description: Collection of classes to draw resonators (such as meander resonators).

from qiskit_metal.qlibrary.core import QComponent, QRoute
from qiskit_metal.qlibrary.tlines.anchored_path import RouteAnchors
import numpy as np
import shapely
from qiskit_metal.toolbox_python.attr_dict import Dict
from qiskit_metal import draw

class ResonatorMeander(QComponent):
    """Create a CPW Meander Resonator.

    Inherits QComponent class.

    Resonator Metal Geometry and Ground Cutout Pocket:
        * total_length - Total length of meander resonator (as measured when tracing the centre of the stripline)
        * constr_radius - Constraint on the radius of the meander turns
        * constr_width_max - Constraint on the maximum width of the meanders (measured centre to centre of maximal extent of striplines)
        * constr_extend_length - Constraint on maximum distance to extend the resonator
        * trace_width - Width of CPW trace
        * trace_gap - Gap to ground plane around CPW trace
        * fillet_padding - Padding added on every turn to overcome Qiskit-Metal's "bad fillet checks"
        * start_left - If True, the meander turns left first.
    
    The number of turns is automatically calculated. Only 2 out of the three out of the 3 constraints may be non-zero at any one time
    (zero implies that said 'constraint' is a free parameter). Extra padding on the ends of the resonator is added if the geometry cannot
    be solved exactly when constr_extend_length is given. Nonetheless, the overall length will always be satisfied.
        
    As usual, the positioning can be done dynamically as a vector given by the supplied parameters: pos_x,pos_y, and orientation
        
    Pins:
        There are two pins on the resonator on either end (a and b)

    Sketch:
        Below is a sketch of the resonator
        ::

          LLLLLLLLLLLLLLLLLLLLLLLLLLL
                  _____
                r/  P  \ r              W       r = constr_radius
                |       |               W       W = constr_width_max
          _____/r       |        ____   W       L = constr_extend_length
          P             |      r/   P   W
                        |      |        W       P = fillet_padding - locations where the padding occurs along the resonator length
                       r\_____/r        W
                           P

    .. image::
        Cap3Interdigital.png

    .. meta::
        Cap 3 Interdigital

    Default Options:
        * pos_x='0um'
        * pos_y='0um'
        * orientation=0
        * total_length='1mm'
        * constr_radius='30um'
        * constr_width_max='250um'
        * constr_extend_length='0um'
        * trace_width='10um',
        * trace_gap='10um'
        * fillet_padding='2um'
        * start_left=True
    """

    default_options = Dict(pos_x='0um', pos_y='0um', orientation=0,
                           total_length='1mm',
                           constr_radius='30um',
                           constr_width_max='250um',
                           constr_extend_length='0um',
                           trace_width='10um', trace_gap='10um',
                           fillet_padding='2um',
                           start_left=True
                           )

    def make(self):
        p = self.p

        l = p.total_length
        L = p.constr_extend_length
        r = p.constr_radius
        wid_constr = p.constr_width_max
        f = p.fillet_padding

        lePath, pin1, pin2, r = ResonatorMeander.draw_resonator(l, L, r, wid_constr, f, p.start_left, p.trace_width, p.pos_x, p.pos_y, p.orientation)

        self.add_qgeometry('path', {'trace': lePath},
                           width=p.trace_width,
                           fillet=r,
                           layer=p.layer)
        self.add_qgeometry('path', {'trace_gap': lePath},
                           width=p.trace_width + 2*p.trace_gap,
                           fillet=r,
                           subtract=True,
                           layer=p.layer)
        self.add_pin('a', pin1.coords, width=p.trace_width)
        self.add_pin('b', pin2.coords, width=p.trace_width)

    @staticmethod
    def draw_resonator(l, L, r, wid_constr, f, start_left, trace_width, x, y, angle):
        err_msg = "Overconstrained the meander. Only 2 of the 3 constraints may be set: constr_radius, constr_width_max and constr_extend_length. Currently constr_wid_max and constr_extend_length will not work as the chosen constraints - other combinations will."
        if r > 0 and wid_constr > 0:
            assert L == 0, err_msg
            smax = wid_constr - 2*r
            N = np.ceil( (l - 3*f - smax - 2*np.pi*r + 2*r) / (np.pi*r + f + smax) )
            N = int(N)
            assert N > 0, "Constraints on meander cannot be satisfied. Try reducing fillet_padding or constr_radius. Otherwise try increasing constr_width_max."
            s = (l - 3*f - 2*np.pi*r + 2*r - N*(np.pi*r + f)) / (1+N)
            resid=0
        elif r > 0 and L > 0:
            assert wid_constr == 0, err_msg
            N = int( (L-3*f-4*r) / (f+2*r) )
            if N > 0:
                L0 = 3*f + 4*r + N*(f+2*r)
                resid=(L-L0)*0.5
                s = (l-2*resid - 3*f - 2*np.pi*r + 2*r - N*(np.pi*r + f)) / (1+N)
            else:
                assert N==0 and L==l, "Constraints on meander for this given point-to-point distance cannot be satisfied. Try reducing fillet_padding or constr_radius."
                resid=0
        elif wid_constr > 0 and L > 0:
            Nmax = int((L-3*f)/f)
            if Nmax > 0:
                best_width = -1
                cur_resid_s_r_N = None
                for N in range(1,Nmax+1):
                    r = (L - 3*f - N*f)/(4 + 2*N)
                    if r <= trace_width/2:
                        continue
                    L0 = 3*f + 4*r + N*(f+2*r)
                    resid=(L-L0)*0.5
                    s = (l-2*resid - 3*f - 2*np.pi*r + 2*r - N*(np.pi*r + f)) / (1+N)
                    cur_w = s+2*r
                    # print(N, cur_w, r, s)
                    if cur_w > wid_constr or s/2 <= r:
                        continue
                    if cur_w > best_width:
                        best_width = cur_w
                        cur_resid_s_r_N = (resid, s, r, N)
                assert best_width > 0, "Constraints on meander cannot be satisfied for some reason. Try increasing wid_constr or decreasing fillet_padding."
                resid, s, r, N = cur_resid_s_r_N
            else:
                assert Nmax==0 and L==l, "Constraints on meander for this given point-to-point distance cannot be satisfied."
                resid=0


        # print(N, s, r, resid, L, trace_width)
        # print(L, r, wid_constr)

        if N == 0:
            lePath = [[0,0], [l,0]]
        else:
            assert s/2 > r, "Constraints on meander cannot form a curved fillet/radius on the meander. The width and radius combination make it too narrow."

            lePath = [[0,0], [resid+f+r,0]]
            cur_x = resid+f+r
            if start_left:
                dir = 1
            else:
                dir=-1
            for m in range(N+1):
                lePath += [[cur_x,dir*(s*0.5+r)], [cur_x+2*r+f,dir*(s*0.5+r)]]
                cur_x = cur_x+2*r+f
                dir *= -1
            lePath += [[cur_x, 0], [cur_x+r+f+resid, 0]]

        pin1 = [[0, -trace_width*0.5], [0, trace_width*0.5]]
        pin2 = [[lePath[-1][0],trace_width*0.5], [lePath[-1][0],-trace_width*0.5]]
        pin1 = shapely.LineString(pin1)
        pin2 = shapely.LineString(pin2)
        lePath = shapely.LineString(lePath)

        polys = [lePath, pin1, pin2]
        polys = draw.rotate(polys, angle, origin=(0, 0))
        polys = draw.translate(polys, x, y)
        [lePath, pin1, pin2] = polys
        return lePath, pin1, pin2, r

class ResonatorMeanderPin(QComponent):
    """Create a CPW Meander Resonator.

    Inherits QComponent class.

    Resonator Metal Geometry and Ground Cutout Pocket:
        * total_length - Total length of meander resonator (as measured when tracing the centre of the stripline)
        * constr_radius - Constraint on the radius of the meander turns
        * constr_width_max - Constraint on the maximum width of the meanders (measured centre to centre of maximal extent of striplines)
        * constr_extend_length - Constraint on maximum distance to extend the resonator
        * trace_width - Width of CPW trace
        * trace_gap - Gap to ground plane around CPW trace
        * fillet_padding - Padding added on every turn to overcome Qiskit-Metal's "bad fillet checks"
        * start_left - If True, the meander turns left first.
    
    The number of turns is automatically calculated. Only 2 out of the three out of the 3 constraints may be non-zero at any one time
    (zero implies that said 'constraint' is a free parameter). Extra padding on the ends of the resonator is added if the geometry cannot
    be solved exactly when constr_extend_length is given. Nonetheless, the overall length will always be satisfied.
        
    The positioning can be done dynamically via:
        * pin_inputs=Dict(start_pin=Dict(component=f'...',pin='...')) - Specifying start position via a component pin
    The resulting resonator is right in the centre. This class ignores pos_x, pos_y and orientation...
        
    Pins:
        There are two pins on the resonator on either end (a and b)

    Sketch:
        Below is a sketch of the resonator
        ::

          LLLLLLLLLLLLLLLLLLLLLLLLLLL
                  _____
                r/  P  \ r              W       r = constr_radius
                |       |               W       W = constr_width_max
          _____/r       |        ____   W       L = constr_extend_length
          P             |      r/   P   W
                        |      |        W       P = fillet_padding - locations where the padding occurs along the resonator length
                       r\_____/r        W
                           P

    .. image::
        Cap3Interdigital.png

    .. meta::
        Cap 3 Interdigital

    Default Options:
        * pos_x='0um'
        * pos_y='0um'
        * orientation=0
        * total_length='1mm'
        * constr_radius='50um'
        * constr_width_max='250um'
        * constr_extend_length='0um'
        * trace_width='10um',
        * trace_gap='10um'
        * fillet_padding='2um'
        * start_left=True
    """

    default_options = Dict(pos_x='0um', pos_y='0um', orientation=0,
                           total_length='1mm',
                           constr_radius='50um',
                           constr_width_max='250um',
                           constr_extend_length='0um',
                           trace_width='10um', trace_gap='10um',
                           fillet_padding='2um',
                           start_left=True
                           )

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
        p = self.p

        start_point = self.design.components[self.options.pin_inputs.start_pin.component].pins[self.options.pin_inputs.start_pin.pin]
        startPt = start_point['middle']
        norm = start_point['normal']

        l = p.total_length
        L = p.constr_extend_length
        r = p.constr_radius
        wid_constr = p.constr_width_max
        f = p.fillet_padding

        lePath, pin1, pin2, r = ResonatorMeander.draw_resonator(l, L, r, wid_constr, f, p.start_left, p.trace_width, startPt[0], startPt[1], np.arctan2(norm[1],norm[0])/np.pi*180)

        self.add_qgeometry('path', {'trace': lePath},
                           width=p.trace_width,
                           fillet=r,
                           layer=p.layer)
        self.add_qgeometry('path', {'trace_gap': lePath},
                           width=p.trace_width + 2*p.trace_gap,
                           fillet=r,
                           subtract=True,
                           layer=p.layer)
        self.add_pin('a', pin1.coords, width=p.trace_width)
        self.add_pin('b', pin2.coords, width=p.trace_width)

class ResonatorMeanderPinPin(QComponent):
    """Create a CPW Meander Resonator.

    Inherits QComponent class.

    Resonator Metal Geometry and Ground Cutout Pocket:
        * total_length - Total length of meander resonator (as measured when tracing the centre of the stripline)
        * constr_radius - Constraint on the radius of the meander turns
        * constr_width_max - Constraint on the maximum width of the meanders (measured centre to centre of maximal extent of striplines)
        * trace_width - Width of CPW trace
        * trace_gap - Gap to ground plane around CPW trace
        * fillet_padding - Padding added on every turn to overcome Qiskit-Metal's "bad fillet checks"
        * start_left - If True, the meander turns left first.
    
    The number of turns is automatically calculated. Extra padding on the ends of the resonator is added if the geometry cannot
    be solved exactly. Nonetheless, the overall length will always be satisfied.
        
    The positioning can be done dynamically via:
        * pin_inputs=Dict(start_pin=Dict(component=f'...',pin='...'), end_pin=Dict(component=f'...',pin='...')) - Specifying start and end
          positions via a component pins
    This class ignores pos_x, pos_y and orientation...
        
    Pins:
        There are two pins on the resonator on either end (a and b)

    Sketch:
        Below is a sketch of the resonator
        ::

                  _____
                r/  P  \ r              W       r = constr_radius
                |       |               W       W = constr_width_max
          _____/r       |        ____   W
          P             |      r/   P   W
                        |      |        W       P = fillet_padding - locations where the padding occurs along the resonator length
                       r\_____/r        W
                           P

    .. image::
        Cap3Interdigital.png

    .. meta::
        Cap 3 Interdigital

    Default Options:
        * total_length='1mm'
        * constr_radius='50um'
        * constr_width_max='250um'
        * trace_width='10um',
        * trace_gap='10um'
        * fillet_padding='2um'
        * start_left=True
    """

    default_options = Dict(pos_x='0um', pos_y='0um', orientation=0,
                           total_length='1mm',
                           constr_radius='50um',
                           constr_width_max='0um',
                           trace_width='10um', trace_gap='10um',
                           fillet_padding='2um',
                           start_left=True
                           )

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
        p = self.p

        start_point = self.design.components[self.options.pin_inputs.start_pin.component].pins[self.options.pin_inputs.start_pin.pin]
        startPt = start_point['middle']

        end_point = self.design.components[self.options.pin_inputs.end_pin.component].pins[self.options.pin_inputs.end_pin.pin]
        endPt = end_point['middle']

        vec = endPt-startPt
        ori = np.arctan2(vec[1],vec[0])/np.pi*180

        l = p.total_length
        L = np.linalg.norm(vec)
        r = p.constr_radius
        wid_constr = p.constr_width_max
        f = p.fillet_padding

        lePath, pin1, pin2, r = ResonatorMeander.draw_resonator(l, L, r, wid_constr, f, p.start_left, p.trace_width, startPt[0], startPt[1], ori)

        self.add_qgeometry('path', {'trace': lePath},
                           width=p.trace_width,
                           fillet=r,
                           layer=p.layer)
        self.add_qgeometry('path', {'trace_gap': lePath},
                           width=p.trace_width + 2*p.trace_gap,
                           fillet=r,
                           subtract=True,
                           layer=p.layer)
        self.add_pin('a', pin1.coords, width=p.trace_width)
        self.add_pin('b', pin2.coords, width=p.trace_width)
