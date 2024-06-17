import gdspy
import shapely
import numpy as np
from SQDMetal.Utilities.QUtilities import QUtilities
from SQDMetal.Utilities.QiskitShapelyRenderer import QiskitShapelyRenderer
from SQDMetal.Utilities.ShapelyEx import ShapelyEx

class MakeGDS:
    def __init__(self, design, threshold=1e-9, precision=1e-9, curve_resolution=22, smooth_radius=0):
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
        self.smooth_radius = smooth_radius
        self.gds_units = 1e-6
        self.refresh()

    def _shapely_to_gds(self, metals, layer):
        if isinstance(metals, shapely.geometry.multipolygon.MultiPolygon):
            temp_cur_metals = [x for x in metals.geoms]
        else:
            temp_cur_metals = [metals]
        gds_metal = None
        for cur_poly in temp_cur_metals:
            if gds_metal:
                gds_metal = gdspy.boolean(gdspy.Polygon(cur_poly.exterior.coords[:]), gds_metal, "or", layer=layer, max_points=0)
            else:
                gds_metal = gdspy.Polygon(cur_poly.exterior.coords[:], layer=layer)
            for cur_int in cur_poly.interiors:
                gds_metal = gdspy.boolean(gds_metal,gdspy.Polygon(cur_int.coords[:]), "not", layer=layer, max_points=0)
        return gds_metal

    def _round_corners(self, gds_metal, output_layer):
        #Perform the usual expand and contract to join fractures... But this fails due to some gdsPy bugs?
        # new_metal = gdspy.boolean(gdspy.offset(new_metal, self.threshold/self.gds_units, max_points=0, join_first=True), None, "or", max_points=0)
        # new_metal = gdspy.boolean(gdspy.offset(new_metal, -self.threshold/self.gds_units, max_points=0, join_first=True), None, "or", max_points=0)
        # new_metal = new_metal.fillet(self.smooth_radius/self.gds_units, max_points=0)

        #So use Shapely instead...
        polyFuse = ShapelyEx.fuse_polygons_threshold([shapely.Polygon(x) for x in gds_metal.polygons], self.threshold/self.gds_units)
        rndDist = self.smooth_radius/self.gds_units
        polyFuse = polyFuse.buffer(rndDist*0.5, join_style=1, quad_segs=self.curve_resolution).buffer(-rndDist, join_style=1).buffer(rndDist*0.5, join_style=1, quad_segs=self.curve_resolution)
        new_metal = self._shapely_to_gds(polyFuse, output_layer)
        return new_metal

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

        self.cell = self.lib.new_cell('Main')

        threshold = self.threshold

        #Assemble geometry for each layer
        all_layers = []
        self._layer_metals = {}
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
        chip = shapely.affinity.scale(chip, xfact=leUnits/self.gds_units, yfact=leUnits/self.gds_units, origin=(0,0))
        all_layers.append((0, chip))
        #
        #
        leLayers = np.unique(gsdf['layer'])
        for layer in leLayers:
            metals = gsdf.loc[(gsdf['layer'] == layer) & (gsdf['subtract'] == False)]
            metals = ShapelyEx.fuse_polygons_threshold(metals['geometry'].tolist(), threshold)
            metals = shapely.affinity.scale(metals, xfact=leUnits/self.gds_units, yfact=leUnits/self.gds_units, origin=(0,0))
            all_layers.append((layer, metals))

        neg_layer_offset = np.max(leLayers)+10
        for cur_layer in all_layers:
            layer, metals = cur_layer

            gds_metal = self._shapely_to_gds(metals, layer)

            gds_rect_neg = gdspy.boolean(gdspy.Rectangle(((cx-0.5*sx)*leUnits/self.gds_units, (cy-0.5*sy)*leUnits/self.gds_units),
                                                        ((cx+0.5*sx)*leUnits/self.gds_units, (cy+0.5*sy)*leUnits/self.gds_units)),
                                        gds_metal, "not", layer=layer+neg_layer_offset, max_points=0)

            if self.smooth_radius > 0:
                gds_metal = self._round_corners(gds_metal, layer)
            self.cell.add(gdspy.boolean(gds_metal, gds_metal, "or", layer=layer))    #max_points = 199 here...
            self._layer_metals[layer] = gds_metal
            if gds_rect_neg:
                self.cell.add(gdspy.boolean(gds_rect_neg, gds_rect_neg, "or", layer=layer+neg_layer_offset))    #max_points = 199 here...
                self._layer_metals[layer+neg_layer_offset] = gds_rect_neg

    def add_boolean_layer(self, layer1_ind, layer2_ind, operation, output_layer=None):
        assert operation in ["and", "or", "xor", "not"], "Operation must be: \"and\", \"or\", \"xor\", \"not\""
        if output_layer == None:
            output_layer = max([x for x in self._layer_metals]) + 1

        new_metal = gdspy.boolean(self._layer_metals[layer1_ind], self._layer_metals[layer2_ind], operation, layer=output_layer, precision=1e-6, max_points=0)
        if self.smooth_radius > 0:
            new_metal = self._round_corners(new_metal, output_layer)

        new_metal = gdspy.boolean(new_metal, new_metal, "or", layer=output_layer)    #max_points = 199 here...

        self.cell.add(new_metal)
        self._layer_metals[output_layer] = new_metal

    def export(self, file_name):
        self.lib.write_gds(file_name)
