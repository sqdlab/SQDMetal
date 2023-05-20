# -*- coding: utf-8 -*-
# Author: Prasanna Pakkiam
# Creation Date: 01/05/2023
# Description: Collection of classes to create bounding objects/planes.

from qiskit_metal import draw
from qiskit_metal.toolbox_python.attr_dict import Dict
from qiskit_metal.qlibrary.core import QComponent
import numpy as np
import shapely
from SQDMetal.Utilities.CpwParams import CpwParams
from SQDMetal.Utilities.QUtilities import QUtilities

class BoundRectangle(QComponent):
    """Creates a rectangular bounding box around a given set of components.

    Inherits QComponent class.

    Rectangle created via:
        * bnd_objects - A comma-separated string of objects to which this rectangle is to tightly bound.
        * is_ground_cut - If True, this bounding box becomes a ground cut-out; otherwise, it's a metal piece.

    Default Options:
        * bnd_objects: ''
        * is_ground_cut=True
    """

    #  Define structure functions

    default_options = Dict(bnd_objects='',
                           is_ground_cut=True)
    """Default drawing options"""

    TOOLTIP = """Bounding box object."""

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################

        str_objs = self.options.bnd_objects

        objs = str_objs.split(',')
        objs = [x.strip() for x in objs]
        minX, minY, maxX, maxY = QUtilities.get_comp_bounds(self.design, objs)
        poly = shapely.Polygon([(minX, minY), (maxX, minY), (maxX, maxY), (minX, maxY)])

        #subtracts out ground plane on the layer its on
        self.add_qgeometry('poly',
                           dict(bndRect=poly),
                           subtract=p.is_ground_cut,
                           layer=p.layer)
