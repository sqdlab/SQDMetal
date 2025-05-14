# -*- coding: utf-8 -*-
# Author: Prasanna Pakkiam
# Creation Date: 29/05/2023
# Description: Collection of classes to draw Xmon qubits.

from qiskit_metal import draw
from qiskit_metal.toolbox_python.attr_dict import Dict
from qiskit_metal.qlibrary.core import QComponent
import numpy as np
import shapely

from matplotlib.path import Path
import matplotlib.path as mpath
from matplotlib.textpath import TextPath, TextToPath
from matplotlib.font_manager import FontProperties
import matplotlib as mpl


class LabelText(QComponent):
    """Create an Text Label

    Inherits QComponent class.

    The text label has the following attributes:
        * text  - The actual text to display
        * font_size - Vertical height of the text's font (in the length units used in qiskit-metal)
        * is_latex - (default False) If True, the LaTeX parser in matplotlib is used to process the text
        * style - Styles the text via the matplotlib options: must be 'normal', 'italic', 'oblique'
        * font - Font family (irrelevant if is_latex is True)
        * is_gap - (Default False) If True, the text will be a ground plane cut rather than metal
        * centre_text - (Default False) If True, then the positional coordinate is at the centre of the text rather than the bottom-left-hand corner

    The positioning can be done via position (pos_x,pos_y) and angle given by 'orientation'.

    Pins:
        There are no pins present in this label

    .. image::
        Cap3Interdigital.png

    .. meta::
        Xmon cross

    Default Options:
        * pos_x='0um'
        * pos_y='0um'
        * orientation=0
        * text='text'
        * font_size='100um'
        * is_latex=False
        * style="normal"
        * font="sans-serif"
        * is_gap=False
        * centre_text=False
    """

    default_options = Dict(pos_x='0um',pos_y='0um',
                           orientation=0,
                           text='text',
                           font_size='100um',
                           is_latex=False,
                           style="normal",
                           font="sans-serif",
                           is_gap=False,
                           centre_text=False
                           )

    def _matplotlib_path_to_polygons(self, path):
        codes_old = path.codes[path.codes != mpath.Path.CLOSEPOLY]
        verts_old = path.vertices[path.codes != mpath.Path.CLOSEPOLY]
        
        # Find indices of MOVETO commands to split subpaths
        subpath_start_indices = np.append(np.flatnonzero(codes_old == mpath.Path.MOVETO), len(codes_old))
        
        polygons_ext = []
        polygons_int = []
        for i in range(len(subpath_start_indices) - 1):
            start = subpath_start_indices[i]
            end = subpath_start_indices[i + 1]
            
            # Extract vertices for the current subpath
            subpath_verts = verts_old[start:end]
            
            # Create a Shapely Polygon from the subpath vertices
            if len(subpath_verts) > 2:  # Ensure enough points to form a polygon
                lePoly = shapely.Polygon(subpath_verts[:-1])
                if lePoly.exterior.is_ccw:
                    polygons_int.append(lePoly)
                else:
                    polygons_ext.append(lePoly)

        poly_full = shapely.unary_union(polygons_ext)
        if len(polygons_int) > 0:
            poly_ints = shapely.unary_union(polygons_int)
            poly_full = shapely.difference(poly_full, poly_ints)
        return poly_full

    def make(self):
        """This is executed by the user to generate the qgeometry for the
        component."""
        p = self.p
        #########################################################

        allowed_styles = ['normal', 'italic', 'oblique']
        assert p.style in allowed_styles, "The argument 'style' must be 'normal', 'italic', 'oblique'"

        #TODO: Look up how to smooth the text...
        # prevSimplify = mpl.rcParams['path.simplify']
        # mpl.rcParams['path.simplify_threshold'] = 1.0
        # mpl.rcParams['path.simplify'] = False

        fp = FontProperties(family=p.font, style=p.style)
        leText2PathObj = TextToPath()
        leText2PathObj.DPI = 300
        verts, codes = leText2PathObj.get_text_path(fp, self.options.text, ismath=p.is_latex)
        lePath = Path(verts, codes, closed=False)

        lePoly = self._matplotlib_path_to_polygons(lePath)
        xMin, yMin, xMax, yMax = lePoly.bounds

        polys = [lePoly]
        if p.centre_text:
            polys = draw.translate(polys, -(xMin+xMax)/2, -(yMin+yMax)/2)
        polys = draw.scale(polys, p.font_size/300, p.font_size/300, origin=(0,0))
        polys = draw.rotate(polys, p.orientation, origin=(0, 0))
        polys = draw.translate(polys, p.pos_x, p.pos_y)
        lePoly = polys[0]   #Remove [0] if there are multiple polygons later...

        # Adds the object to the qgeometry table
        self.add_qgeometry('poly',
                           dict(textLabel=lePoly),
                           layer=p.layer,
                           subtract=p.is_gap)

        # #subtracts out ground plane on the layer it's on
        # self.add_qgeometry('poly',
        #                    dict(gap=gap),
        #                    subtract=True,
        #                    layer=p.layer)

        # Generates its own pins
        # self.add_pin('up', pin_up.coords[::-1], width=p.vBar_width)
        # self.add_pin('down', pin_down.coords[::-1], width=p.vBar_width)
        # self.add_pin('left', pin_left.coords[::-1], width=p.hBar_width)
        # self.add_pin('right', pin_right.coords[::-1], width=p.hBar_width)
