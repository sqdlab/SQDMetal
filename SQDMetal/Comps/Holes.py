# -*- coding: utf-8 -*-
# Author: Prasanna Pakkiam
# Creation Date: 12/06/2023
# Description: Collection of classes to draw holes.

from qiskit_metal import draw
from qiskit_metal.toolbox_python.attr_dict import Dict
from qiskit_metal.qlibrary.core import QComponent
import numpy as np
import shapely
from SQDMetal.Utilities.QiskitShapelyRenderer import QiskitShapelyRenderer
from SQDMetal.Utilities.QUtilities import QUtilities
from scipy.spatial import KDTree
from SQDMetal.Utilities.ShapelyEx import ShapelyEx

class HoleBorders(QComponent):
    """Draws holes some distance away from metallic objects - similar to via fencing in PCBs.

    Inherits QComponent class.

    The holes are placed on the given layer (set by 'layer' option) are configured via:
        * hole_radius - Radius of each individual hole drawn as a metallic object on the given layer
        * dist_holes  - Distance spacing between holes - both along a given border and between multiple borders of holes
        * dist_init   - Distance between any metallic object and its border of holes (i.e. the spacing of holes to metals)
        * num_hole_lines  - Number of borders of holes to place around any given metallic object
        * segs_per_circle - Number of segments per circle to use when drawing the holes (i.e. they are just regular polygons)
        * dist_min  - After setting the bordering holes, they are culled if they are too close (e.g. two borders of nearby metals
                      intersecting one another). The threshold distance is set by dist_min. On curved surfaces, it's advisable to
                      set this to be smaller than dist_holes
        * hole_radius_gnd - If above zero, there is a ground cut-out of holes with radius hole_radius_gnd. Otherwise, there is no
                            ground cutout below the holes.
        * border_layers  - Given as a list of integers signifying the layers upon which to consider when patterning the bordering
                           holes.
        * exclude_geoms  - Given a list of component names to exclude when patterning the bordering holes.
        * include_geoms  - Given a list of component names to include when patterning the bordering holes. If it's empty, then all
                           metals in the layers given by border_layers are selected for border-hole patterning.
        * layer_gnd_plane - The layer in which to place the ground cutout. This should be typically layer 1 (hence, the default).
        * exclude_gaps   - (Default False) If True, the holes cannot fall within gaps...
    ***READ CAREFULLY*** the holes are placed as follows:
        1. Select all metals (not cut-outs unless exclude_gaps is True) in all layers supplied by border_layers - restrict to only
           include_geoms if non-empty. Careful when supplying include_geoms as it can pattern holes on non-included metallic surfaces.
        2. Pattern holes on selected metals
        3. Using exclude_geoms, a buffer area is made via dist_init and num_hole_lines. All holes in this area are culled.

    This class ignores pos_x and pos_y...
        
    Pins:
        There are no pins in this object.

    Sketch:
        Below is a sketch of the holes and their placement options.
        ::
            O   O   O   O   O   O   O
                                       N   # = metals to which the holes are bordering 
            O   O   O   O   O   O   O      d = dist_holes
                i                          i = dist_init
                i                          N = num_hole_lines (2 in this example)
            ##############        O   O     
            ##############     O   
            i             iiiii     O
            i       
            O   O   O   O   O   O   O d
                                      d
            O   O   O   O   O   O   O d
                             ddd

    .. image::
        Cap3Interdigital.png

    .. meta::
        Holes bordering on metals

    Default Options:
        * hole_radius='5um'
        * dist_holes='20um'
        * dist_init='30um'
        * num_hole_lines=3
        * segs_per_circle=12
        * dist_min='10um'
        * hole_radius_gnd='6um'
        * border_layers=[1,2,3]
        * exclude_geoms=[]
        * include_geoms=[]
        * layer_gnd_plane=1
        * exclude_gaps=False
    """

    default_options = Dict(hole_radius='5um', 
                           dist_holes='20um', 
                           dist_init='30um', 
                           num_hole_lines=3, 
                           segs_per_circle=12, 
                           dist_min='10um', 
                           hole_radius_gnd='6um', 
                           border_layers=[1,2,3],
                           exclude_geoms=[],
                           include_geoms=[],
                           layer_gnd_plane=1,
                           exclude_gaps=False)

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################

        qmpl = QiskitShapelyRenderer(None, self.design, None)
        gsdf = qmpl.get_net_coordinates()
        gsdf = gsdf[gsdf['layer'].isin(p.border_layers)]
        if p.exclude_gaps:
            filt = gsdf
        else:
            filt = gsdf.loc[gsdf['subtract'] == False]

        if len(p.exclude_geoms) > 0:
            leGeoms = [self.design.components[x].id for x in self.options.exclude_geoms]
            exfilt = filt[filt['component'].isin(leGeoms)]
        if len(p.include_geoms) > 0:
            leGeoms = [self.design.components[x].id for x in self.options.include_geoms]
            filt = filt[filt['component'].isin(leGeoms)]

        units = QUtilities.get_units(self.design)
        metal_polys = shapely.unary_union(filt['geometry'])
        metal_polys = ShapelyEx.fuse_polygons_threshold(metal_polys, 1e-12/units)

        leHoles = []
        for hl in range(int(p.num_hole_lines)):
            polys = metal_polys.buffer(p.dist_init + hl*p.dist_holes)

            lines = []
            if isinstance(polys, shapely.geometry.multipolygon.MultiPolygon):
                for cur_poly in list(polys.geoms):
                    lines += [cur_poly.exterior.coords[:]]
                    for cur_int in cur_poly.interiors:
                        lines += [cur_int.coords[:]]
            else:
                lines += [polys.exterior.coords[:]]
                for cur_int in polys.interiors:
                    lines += [cur_int.coords[:]]

            lines = [shapely.LineString(x) for x in lines]
            for cur_line in lines:
                distances = np.arange(0, cur_line.length, p.dist_holes)
                points = [cur_line.interpolate(distance) for distance in distances]
                leHoles += points

        #Cull holes in excluded geometry
        if len(p.exclude_geoms) > 0:
            ex_metal_polys = shapely.unary_union(exfilt['geometry'])
            ex_metal_polys = ShapelyEx.fuse_polygons_threshold(ex_metal_polys, 1e-12/units)
            polys = ex_metal_polys.buffer(p.dist_init + p.num_hole_lines*p.dist_holes)
            leHoleFilts = []
            for m, cur_hole in enumerate(leHoles):
                if not polys.contains(cur_hole):
                    leHoleFilts += [m]
            leHoles = [leHoles[m] for m in leHoleFilts]

        leHoles = np.array([[h.x,h.y] for h in leHoles])

        #Cull holes that are too close
        #Using 'fast' solution from: https://stackoverflow.com/questions/58114650/np-where-to-eliminate-data-where-coordinates-are-too-close-to-each-other
        X_tree = KDTree(leHoles)
        in_radius = np.array(list(X_tree.query_pairs(p.dist_min))).flatten()
        leHoles = leHoles[np.where(~np.in1d(np.arange(leHoles.shape[0]), in_radius))[0]]

        #Convert holes into n-polygons
        polys = []
        for cur_hole in leHoles:
            polys.append((
                tuple( [(cur_hole[0]+p.hole_radius*np.cos(angl), cur_hole[1]+p.hole_radius*np.sin(angl)) for angl in np.linspace(0,2*np.pi, p.segs_per_circle)] ),
                []))

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(holes=shapely.MultiPolygon(polys)),
                           layer=p.layer)

        if p.hole_radius_gnd > 0:
            polys = []
            for cur_hole in leHoles:
                polys.append((
                    tuple( [(cur_hole[0]+p.hole_radius_gnd*np.cos(angl), cur_hole[1]+p.hole_radius_gnd*np.sin(angl)) for angl in np.linspace(0,2*np.pi, p.segs_per_circle)] ),
                    []))
            self.add_qgeometry('poly',
                               dict(holesGND=shapely.MultiPolygon(polys)),
                               subtract=True,
                               layer=p.layer_gnd_plane)

