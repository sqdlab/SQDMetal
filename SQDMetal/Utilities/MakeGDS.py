import gdspy
import shapely
import numpy as np
from SQDMetal.Utilities.QUtilities import QUtilities
from SQDMetal.Utilities.QiskitShapelyRenderer import QiskitShapelyRenderer
from SQDMetal.Utilities.ShapelyEx import ShapelyEx

class MakeGDS:
    def __init__(self, design, threshold=1e-9, precision=1e-9, curve_resolution=22):
        '''
        Inputs:
            design      - The QiskitMetal design object...
            threshold   - The smallest feature width that can exist; anything smaller will get culled out. This is to help remove
                          artefacts from floating-point inaccuracies. It is given in metres.
            precision   - The GDS granularity/precision given in metres.
            curve_resolution - Number of vertices to use on a path fillet
        '''
        self._design = design
        self.curve_resolution = curve_resolution
        self.threshold = threshold
        self.precision = precision
        self.refresh()

    def refresh(self):
        '''
        Deletes everything and rebuilds metallic layers from scratch with the ground plane being in layer zero. The negative versions
        of the layers are given in layers starting from highest layer index plus 10.
        '''
        qmpl = QiskitShapelyRenderer(None, self._design, None)
        gsdf = qmpl.get_net_coordinates(self.curve_resolution)

        gdspy.current_library.cells.clear()

        leUnits = QUtilities.get_units(self._design)
        self.lib = gdspy.GdsLibrary(precision=self.precision)

        cell = self.lib.new_cell('Main')

        threshold = self.threshold
        gds_units = 1e-6

        #Assemble geometry for each layer
        all_layers = []
        #Assuming Layer 0 is GND Plane
        gaps = gsdf.loc[(gsdf['subtract'] == True)]
        gaps = ShapelyEx.fuse_polygons_threshold(gaps['geometry'].tolist(), threshold)
        #
        sx = QUtilities.parse_value_length(self._design.chips['main']['size']['size_x'])/leUnits
        sy = QUtilities.parse_value_length(self._design.chips['main']['size']['size_y'])/leUnits
        cx, cy = [QUtilities.parse_value_length(self._design.chips['main']['size'][x])/leUnits for x in ['center_x', 'center_y']]
        #
        rect = shapely.Polygon([(cx-0.5*sx, cy-0.5*sy), (cx+0.5*sx, cy-0.5*sy), (cx+0.5*sx, cy+0.5*sy), (cx-0.5*sx, cy+0.5*sy), (cx-0.5*sx, cy-0.5*sy)])
        chip = rect.difference(gaps)
        chip = shapely.affinity.scale(chip, xfact=leUnits/gds_units, yfact=leUnits/gds_units, origin=(0,0))
        all_layers.append((0, chip))
        #
        #
        leLayers = np.unique(gsdf['layer'])
        for layer in leLayers:
            metals = gsdf.loc[(gsdf['layer'] == layer) & (gsdf['subtract'] == False)]
            metals = ShapelyEx.fuse_polygons_threshold(metals['geometry'].tolist(), threshold)
            metals = shapely.affinity.scale(metals, xfact=leUnits/gds_units, yfact=leUnits/gds_units, origin=(0,0))
            all_layers.append((layer, metals))

        neg_layer_offset = np.max(leLayers)+10
        for cur_layer in all_layers:
            layer, metals = cur_layer

            if isinstance(metals, shapely.geometry.multipolygon.MultiPolygon):
                temp_cur_metals = [x for x in metals.geoms]
            else:
                temp_cur_metals = [metals]
            gds_metal = None
            for cur_poly in temp_cur_metals:
                if gds_metal:
                    gds_metal = gdspy.boolean(gdspy.Polygon(cur_poly.exterior.coords[:]), gds_metal, "or", layer=layer)
                else:
                    gds_metal = gdspy.Polygon(cur_poly.exterior.coords[:], layer=layer)
                for cur_int in cur_poly.interiors:
                    gds_metal = gdspy.boolean(gds_metal,gdspy.Polygon(cur_int.coords[:]), "not", layer=layer)

            gds_rect_neg = gdspy.boolean(gdspy.Rectangle(((cx-0.5*sx)*leUnits/gds_units, (cy-0.5*sy)*leUnits/gds_units),
                                                        ((cx+0.5*sx)*leUnits/gds_units, (cy+0.5*sy)*leUnits/gds_units)),
                                        gds_metal, "not", layer=layer+neg_layer_offset)

            cell.add(gds_metal)
            if gds_rect_neg:
                cell.add(gds_rect_neg)

    def export(self, file_name):
        self.lib.write_gds(file_name)
