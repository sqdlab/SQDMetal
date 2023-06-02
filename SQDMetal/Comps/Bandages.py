from qiskit_metal import draw
from qiskit_metal.toolbox_python.attr_dict import Dict
from qiskit_metal.qlibrary.core import QComponent
import numpy as np
import shapely

class BandageRectPin(QComponent):

    default_options = Dict(width='10um',
                           height='10um'
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
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################

        start_point = self.design.components[self.options.pin_inputs.start_pin.component].pins[self.options.pin_inputs.start_pin.pin]
        startPt = start_point['middle']
        normal = -start_point['normal']
        rot_angle = np.arctan2(normal[1], normal[0])

        bandaid = [
                  (-p.width*0.5, -p.height*0.5),
                  (-p.width*0.5, p.height*0.5),
                  (p.width*0.5, p.height*0.5),
                  (p.width*0.5, -p.height*0.5)]

        bandaid = shapely.Polygon(bandaid)

        polys = [bandaid]
        polys = draw.rotate(polys, rot_angle, origin=(0, 0), use_radians=True)
        polys = draw.translate(polys, *startPt)
        [bandaid] = polys

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(bandaid=bandaid),
                           layer=p.layer)


