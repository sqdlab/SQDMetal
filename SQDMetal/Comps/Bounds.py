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
from SQDMetal.Utilities.QiskitShapelyRenderer import QiskitShapelyRenderer
from SQDMetal.Utilities.ShapelyEx import ShapelyEx

class BoundRectangle(QComponent):
    """Creates a rectangular bounding box around a given set of components.

    Inherits QComponent class.

    Rectangle created via:
        * bnd_objects - A list of objects to which this rectangle is to tightly bound.
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

        obj_names = self.options.bnd_objects
        assert isinstance(obj_names, (list, tuple)), "Must give bnd_objects as a list of strings."

        minX, minY, maxX, maxY = QUtilities.get_comp_bounds(self.design, obj_names)
        poly = shapely.Polygon([(minX, minY), (maxX, minY), (maxX, maxY), (minX, maxY)])

        #subtracts out ground plane on the layer it's on
        self.add_qgeometry('poly',
                           dict(bndRect=poly),
                           subtract=p.is_ground_cut,
                           layer=p.layer)

class BoundGroundShield(QComponent):
    """Specialised operation to draw a perimetric envelope ground shield around a set of objects. This is only really
    applicable if one is not going to pattern the default ground plane.

    Inherits QComponent class.

    The CPW wire parameters are given via:
        * gnd_width - width of the ground envelope
        * include_geoms - list of names of the components to which the ground envelope is to be drawn
        * exclude_geoms - list of names of the components whose regions are to be discarded when drawing the ground shield
        * exclude_pins  - list of tuples given as: (component name, component pin name, exclusion_distance). A box of size
                          exclusion_distance is drawn around these pins. These regions are excluded when drawing the ground
                          shield envelope.
        * curve_resolution - Defaults to 32. This is the number of corners per curved surface used when pre-rendering the
                             geometric structures (i.e. using QiskitShapelyRenderer) when calculating the ground shield
                             envelopes

    Pins:
        There are no inherent pins in this object, but they can be added manually via constructs like RouteJoint

    .. image::
        Cap3Interdigital.png

    .. meta::
        A Ground Shield pattern drawn around a set of objects.

    Default Options:
        * gnd_width='10um'
        * include_geoms=[]
        * exclude_geoms=[]
        * exclude_pins=[]
        * curve_resolution=32
    """

    default_options = Dict(gnd_width='10um',
                           include_geoms=[],
                           exclude_geoms=[],
                           exclude_pins=[],
                           curve_resolution=32)

    def make(self):
        p = self.p

        qmpl = QiskitShapelyRenderer(None, self.design, None)
        gsdf = qmpl.get_net_coordinates(resolution=p.curve_resolution)
        filt = gsdf#.loc[gsdf['subtract'] == True]

        if len(p.exclude_geoms) > 0:
            leGeoms = [self.design.components[x].id for x in self.options.exclude_geoms]
            exfilt = filt[filt['component'].isin(leGeoms)]
        if len(p.include_geoms) > 0:
            leGeoms = [self.design.components[x].id for x in self.options.include_geoms]
            filt = filt[filt['component'].isin(leGeoms)]

        units = QUtilities.get_units(self.design)
        metal_polys = shapely.unary_union(filt['geometry'])
        metal_polys = ShapelyEx.fuse_polygons_threshold(metal_polys, 1e-12/units)

        if len(p.exclude_geoms) > 0:
            ex_metal_polys = shapely.unary_union(exfilt['geometry'])
            ex_metal_polys = ShapelyEx.fuse_polygons_threshold(ex_metal_polys, 1e-12/units)

        polys = metal_polys.buffer(p.gnd_width)
        polys = shapely.difference(polys, metal_polys)

        if len(p.exclude_geoms) > 0:
            polys = shapely.difference(polys, ex_metal_polys)

        exPins = []
        for cur_excl_pin in p.exclude_pins:
            comp_name, pin_name, dist_excl = cur_excl_pin   #The parser already converts units on dist_excl...
            pin_point = self.design.components[comp_name].pins[pin_name]
            startPt = pin_point['middle']        
            vec_norm = pin_point['normal']
            vec_norm /= np.linalg.norm(vec_norm)
            vec_perp = np.array([vec_norm[1], -vec_norm[0]])
            leExclPath = [startPt + vec_perp*dist_excl/2,
                          startPt + vec_perp*dist_excl/2 + vec_norm*dist_excl,
                          startPt - vec_perp*dist_excl/2 + vec_norm*dist_excl,
                          startPt - vec_perp*dist_excl/2]
            exPins.append(shapely.Polygon(leExclPath))
        if len(exPins) > 0:
            exPins = shapely.unary_union(exPins)
            polys = shapely.difference(polys, exPins)

        if isinstance(polys, shapely.geometry.multipolygon.MultiPolygon):
            polys = list(polys.geoms)
        else:
            polys = [polys]
        polyDict = {}
        for m, cur_poly in enumerate(polys):
            polyDict[f'poly{m}'] = cur_poly
        #
        self.add_qgeometry('poly',
                           polyDict,
                           layer=p.layer)

        # polys = metal_polys.buffer(p.gnd_width/2)
        # lines = []
        # if isinstance(polys, shapely.geometry.multipolygon.MultiPolygon):
        #     for cur_poly in list(polys.geoms):
        #         lines += [cur_poly.exterior.coords[:]]
        #         for cur_int in cur_poly.interiors:
        #             lines += [cur_int.coords[:]]
        # else:
        #     lines += [polys.exterior.coords[:]]
        #     for cur_int in polys.interiors:
        #         lines += [cur_int.coords[:]]
        
        # lines = [shapely.LineString(x) for x in lines]

        # for cur_excl_pin in p.exclude_pins:
        #     comp_name, pin_name, dist_excl = cur_excl_pin.comp_name, cur_excl_pin.pin_name, cur_excl_pin.dist_excl
        #     pin_point = self.design.components[comp_name].pins[pin_name]
        #     startPt = pin_point['middle']        
        #     start_norm = pin_point['normal']
        #     startPt

        # for m,cur_line in enumerate(lines):
        #     cur_line = shapely.difference(cur_line, ex_metal_polys)
        #     self.add_qgeometry('path', {'trace': cur_line},
        #                 width=p.gnd_width,
        #                 fillet=0,
        #                 layer=p.layer)



