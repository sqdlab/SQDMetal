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
        minX = 1000
        maxX = -1000
        minY = 1000
        maxY = -1000
        for cur_obj in objs:
            paths = self.design.components[cur_obj].qgeometry_table('path')
            for _, row in paths.iterrows():
                cur_minX, cur_minY, cur_maxX, cur_maxY = row['geometry'].buffer(row['width'] / 2, cap_style=shapely.geometry.CAP_STYLE.flat).bounds
                minX = min(cur_minX, minX)
                minY = min(cur_minY, minY)
                maxX = max(cur_maxX, maxX)
                maxY = max(cur_maxY, maxY)
            for cur_poly in self.design.components[cur_obj].qgeometry_list('poly'):
                cur_minX, cur_minY, cur_maxX, cur_maxY = cur_poly.bounds
                minX = min(cur_minX, minX)
                minY = min(cur_minY, minY)
                maxX = max(cur_maxX, maxX)
                maxY = max(cur_maxY, maxY)
        poly = shapely.Polygon([(minX, minY), (maxX, minY), (maxX, maxY), (minX, maxY)])

        #subtracts out ground plane on the layer its on
        self.add_qgeometry('poly',
                           dict(bndRect=poly),
                           subtract=p.is_ground_cut,
                           layer=p.layer)
