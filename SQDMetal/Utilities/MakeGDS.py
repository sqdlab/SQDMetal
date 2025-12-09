# Copyright 2025 Prasanna Pakkiam
# SPDX-License-Identifier: Apache-2.0

# import gdspy
import gdstk
import shapely
import numpy as np
import os
from datetime import datetime
from types import NoneType
from SQDMetal.Utilities.QUtilities import QUtilities
from SQDMetal.Utilities.QiskitShapelyRenderer import QiskitShapelyRenderer
from SQDMetal.Utilities.ShapelyEx import ShapelyEx

class MakeGDS:
    def __init__(self, design, threshold=1e-9, precision=1e-9, curve_resolution=22, smooth_radius=0, export_type='all', export_layers=None, print_statements=False):
        '''
        Inputs:
            design      - The QiskitMetal design object...
            threshold   - The smallest feature width that can exist; anything smaller will get culled out. This is to help remove
                          artefacts from floating-point inaccuracies. It is given in metres.
            precision   - The GDS granularity/precision given in metres.
            curve_resolution - Number of vertices to use on a path fillet
            export_type - 'all' (exports positive and negative patterns and GND plane as seperate layers), 
                          'positive' (exports flattened metals to layer 0), or 
                          'negative' (exports negative flattened metals to layer 0)
            export_layers - layers to export as a list of integers or single integer (WIP)
        '''
        self._design = design
        self.curve_resolution = curve_resolution
        self.threshold = threshold
        self.precision = precision
        self.smooth_radius = smooth_radius
        self.gds_units = 1e-6
        self.export_type = export_type
        self.export_layers = export_layers
        self.print_statements = print_statements
        self.refresh(full_refresh=True)

    def _shapely_to_gds(self, metals, layer):
        layer = int(layer)
        if isinstance(metals, shapely.geometry.multipolygon.MultiPolygon):
            temp_cur_metals = [x for x in metals.geoms]
        else:
            temp_cur_metals = [metals]
        gds_metal = None
        for cur_poly in temp_cur_metals:
            if gds_metal:
                gds_metal = gdstk.boolean(gdstk.Polygon(cur_poly.exterior.coords[:]), gds_metal, "or", layer=layer, precision=self.precision)
            else:
                gds_metal = gdstk.Polygon(cur_poly.exterior.coords[:], layer=layer)
            for cur_int in cur_poly.interiors:
                gds_metal = gdstk.boolean(gds_metal, gdstk.Polygon(cur_int.coords[:]), "not", layer=layer, precision=self.precision)
        return gds_metal

    def _round_corners(self, gds_metal, output_layer):
        rndDist = self.smooth_radius/self.gds_units
        for p in gds_metal:
            p.fillet(rndDist, tolerance=1e-3)
        new_metal = gdstk.boolean(gds_metal, gds_metal, "or", layer=int(output_layer))
        return new_metal

    def refresh(self, full_refresh=False):
        '''
        Deletes everything and rebuilds metallic layers from scratch with the ground plane being in layer zero.
        The negative versions of the layers are given in layers starting from highest layer index plus 10.
        Export all layers (negative plus positive patterns), or only positive or negative patterns.
        User-defined layers (e.g., via add_boolean_layer) are preserved unless their layer index overlaps with a standard layer being recreated.
        '''
        # Preserve any existing user-defined layers
        existing_user_layers = dict(getattr(self, "_layer_metals", {}))

        # Validate export_type and export_layers before rebuilding
        assert self.export_type in ["all", "negative", "positive"], "Export type must be: \"all\", \"negative\", or \"positive\""

        if self.export_layers is not None:
            self.export_type = "all"
            print(f"Exporting layers {self.export_layers}") if self.print_statements else 0

        # Parse geometry and update metadata before rebuilding geometry
        qmpl = QiskitShapelyRenderer(None, self._design, None)
        gsdf = qmpl.get_net_coordinates(self.curve_resolution)
        leUnits = QUtilities.get_units(self._design)
        self.unit_conv = leUnits / self.gds_units
        threshold = self.threshold
        unit_conv = self.unit_conv
        sx = QUtilities.parse_value_length(self._design.chips['main']['size']['size_x']) / leUnits
        sy = QUtilities.parse_value_length(self._design.chips['main']['size']['size_y']) / leUnits
        cx, cy = [QUtilities.parse_value_length(self._design.chips['main']['size'][x]) / leUnits for x in ['center_x', 'center_y']]
        self.sx, self.sy, self.cx, self.cy = sx, sy, cx, cy
        rect = shapely.Polygon([
            (cx - 0.5 * sx, cy - 0.5 * sy),
            (cx + 0.5 * sx, cy - 0.5 * sy),
            (cx + 0.5 * sx, cy + 0.5 * sy),
            (cx - 0.5 * sx, cy + 0.5 * sy),
            (cx - 0.5 * sx, cy - 0.5 * sy)
        ])
        self.rect = rect

        # Now rebuild library and cell, removing only layers that are about to be recreated
        self.lib = gdstk.Library('GDS Export')
        self.cell = self.lib.new_cell('Main')
        self.lib.replace(self.cell)  # over-write existing cell if it exists
        # Rebuild _layer_metals with only standard layers, preserve user layers below
        self._layer_metals = {}

        # Assemble geometry for each layer
        all_layers = []
        # Assuming Layer 0 is GND Plane
        gaps = gsdf.loc[(gsdf['subtract'])]
        gaps = ShapelyEx.fuse_polygons_threshold(gaps['geometry'].tolist(), threshold)
        chip = rect.difference(gaps)
        chip = shapely.affinity.scale(chip, xfact=unit_conv, yfact=unit_conv, origin=(0, 0))
        all_layers.append((0, chip))
        leLayers = np.unique(gsdf['layer'])  # metal layers
        for layer in leLayers:
            metals = gsdf.loc[(gsdf['layer'] == layer) & (~gsdf['subtract'])]
            metals = ShapelyEx.fuse_polygons_threshold(metals['geometry'].tolist(), threshold)
            metals = shapely.affinity.scale(metals, xfact=unit_conv, yfact=unit_conv, origin=(0, 0))
            all_layers.append((layer, metals))

        print(f'\nExport type: {self.export_type}') if self.print_statements else 0
        neg_layer_offset = np.max(leLayers) + 10
        for cur_layer in all_layers:
            layer, metals = cur_layer
            gds_metal = self._shapely_to_gds(metals, layer)
            if self.smooth_radius > 0:
                gds_metal = self._round_corners(gds_metal, layer)
            print(f" Building layer: {layer:3}") if self.print_statements else 0
            gds_rect_neg = gdstk.boolean(
                gdstk.rectangle(
                    ((cx - 0.5 * sx) * unit_conv, (cy - 0.5 * sy) * unit_conv),
                    ((cx + 0.5 * sx) * unit_conv, (cy + 0.5 * sy) * unit_conv)
                ),
                gds_metal,
                "not",
                layer=int(layer + neg_layer_offset)
            )
            self.cell.add(*gdstk.boolean(gds_rect_neg, gds_rect_neg, "or", layer=int(layer + neg_layer_offset)))
            self._layer_metals[layer + neg_layer_offset] = gds_rect_neg
            print(f" Added negative layer {layer + neg_layer_offset}") if self.print_statements else 0
            self.cell.add(*gdstk.boolean(gds_metal, gds_metal, "or", layer=int(layer)))
            self._layer_metals[layer] = gds_metal

        # Merge in user-defined layers, skipping any that would be overwritten
        if not full_refresh:
            for layer, polys in existing_user_layers.items():
                if layer not in self._layer_metals:
                    self._layer_metals[layer] = polys
                    self.cell.add(*polys)

        # fuse layers for a single-layer design (GND, metal)
        if (len(all_layers) == 2) and (self.export_type in ["positive", "negative"]):
            print(f"  Two layers detected - performing boolean operations for {self.export_type} export") if self.print_statements else 0
            self.add_boolean_layer(11, 12, "and", output_layer=0)
            self.merge_polygons_in_layer(0)
            self.cell.remove(*[poly for poly in self.cell.polygons if poly.layer in [11, 12]])
            gds_metal_neg = self._layer_metals[0]
            if self.export_type == "positive":
                gds_gnd_plane = gdstk.rectangle(
                    ((cx - 0.5 * sx) * unit_conv, (cy - 0.5 * sy) * unit_conv),
                    ((cx + 0.5 * sx) * unit_conv, (cy + 0.5 * sy) * unit_conv)
                )
                gds_positive = gdstk.boolean(gds_gnd_plane, gds_metal_neg, "not", layer=0)
                self.cell.filter([(0, 0)], True)
                self.cell.add(*gdstk.boolean(gds_positive, gds_positive, "or", layer=0))
                self._layer_metals[0] = gds_positive
            self.cell.remove(*[p for p in self.cell.polygons if p.layer == 1])
            self.cell.flatten()
        elif (len(all_layers) > 2) and (self.export_type in ["positive", "negative"]):
            "Positive or negative export for multi-layer designs is not supported."

        # layer exports (works best when export_type=="all")
        if self.export_layers is not None:
            assert isinstance(self.export_layers, list), "Please pass a list of integers as 'export_layers'."
            layer_list = self.get_cell_layers()
            assert set(self.export_layers).issubset(layer_list), "Chosen export layers are not present in the design."
            to_delete = list(set(layer_list) - set(self.export_layers))
            if isinstance(to_delete, list):
                for i in list(to_delete):
                    self.cell.remove(*[p for p in self.cell.polygons if p.layer == i])
            else:
                self.cell.remove(*[p for p in self.cell.polygons if p.layer == to_delete])
        
    def merge_polygons_in_layer(self, layer):
        polys_to_merge = [p for p in self.cell.polygons if p.layer == layer]
        merged = gdstk.boolean(polys_to_merge, [], "or", layer=layer)
        self.cell.remove(*polys_to_merge)
        if merged:
            self.cell.add(*merged)

    def add_boolean_layer(self, layer1_ind, layer2_ind, operation, output_layer=None):
        assert operation in ["and", "or", "xor", "not"], "Operation must be: \"and\", \"or\", \"xor\", \"not\""
        if output_layer is None:
            output_layer = max([x for x in self._layer_metals]) + 1
        try:
            new_metal = gdstk.boolean(self._layer_metals[layer1_ind], self._layer_metals[layer2_ind], operation, layer=int(output_layer), precision=self.precision)
        except KeyError:
            print(f"Warning: Boolean operation ({layer1_ind}, {layer2_ind}) failed: one of the input layers ({layer1_ind}, {layer2_ind}) has no polygons.")
            return output_layer
        if self.smooth_radius > 0:
            new_metal = self._round_corners(new_metal, output_layer)
        new_metal = gdstk.boolean(new_metal, new_metal, "or", layer=int(output_layer), precision=self.precision)
        self.cell.filter([(int(output_layer), 0)], True)
        if not new_metal:
            print(f"No geometry resulted from boolean operation '{operation}' between layers {layer1_ind} and {layer2_ind}.")
            return output_layer
        self.cell.add(*new_metal)
        self._layer_metals[output_layer] = new_metal
        return output_layer
    
    # add a text label (default to bottom left)
    def add_text(self, text_label="", layer=0, size=300, position=None):
        """
        Add a text label to the gds export on a given layer at a given position.

        Inputs:
            text_label  - text label to add
            layer       - layer number to export the text to (if none, write to layer 0)
            size        - text size
            position    - tuple containing position (normalised to 1);
                          e.g. (0,0) is bottom left, (1,1) is top right
        """
        assert isinstance(text_label, str), "Please pass a string as the text label."
        assert isinstance(position, (tuple, NoneType)), "Please pass the position as an (x,y) tuple"
        if isinstance(position, tuple):
            for i in position: 
                assert 0 <= i <= 1, "Ensure the positions are normalised to 1 (i.e. (0,0) is bottom left, (1,1) is top right)"

        layer_list = self.get_cell_layers()
        if (layer in layer_list) and (layer!=0):
            print("Warning: text is being added with no action argument to a (non-ground-plane) layer already containing shapes. Consider changing the layer.")

        # default label
        text_gds = f"Made with SQDMetal {datetime.now().strftime('%Y')}." if text_label == "" else text_label

        # set default position (bottom left) or input position (normalised to 1)
        unit_conv = self.unit_conv
        default_pos = 0.48

        if position is None: 
            text_position = ((self.cx-default_pos*self.sx)*unit_conv, (self.cy-default_pos*self.sy)*unit_conv)
        else: 
            text_position = ((self.cx+(position[0] - 0.5)*self.sx)*unit_conv, (self.cy+(position[1] - 0.5)*self.sy)*unit_conv)
        
        # subtract text from ground plane
        if layer == 0:
            text = gdstk.text(text_gds, size, text_position)
            # subtract from ground plane
            if self.export_type in ['all', 'positive']:
                text = gdstk.text(text_gds, size, text_position)
                text = gdstk.boolean(self._layer_metals[0], text, "not", layer=layer)
                self.cell.filter([(0, 0)], True)
                self.cell.add(*text)
            # add to ground plane
            elif self.export_type=='negative':
                text = gdstk.text(text_gds, size, text_position)
                self.cell.add(*text)
        # add as new layer
        elif layer!=0:
            text = gdstk.text(text_gds, size, text_position, layer=layer)
            self.cell.add(*text)
            self._layer_metals[layer] = text

    def get_cell_layers(self):
        layers = set()
        for poly in self.cell.polygons:
            layers.add(poly.layer)
        for path in self.cell.paths:
            layers.add(path.layer)
        for label in self.cell.labels:
            layers.add(label.layer)
        return sorted(layers)
    
    def delete_layers(self, layers_to_remove):
        if isinstance(layers_to_remove, int):
            layers_to_remove = [layers_to_remove]
        self.cell.filter([(layer, 0) for layer in layers_to_remove])
    
    @staticmethod
    def merge_gds_per_layer_inplace(gds_input_file, inplace=True):
        try:
            import pya # type: ignore
            # Load the layout
            layout = pya.Layout()
            layout.read(f"{gds_input_file}.gds") if gds_input_file[-4:] != ".gds" else layout.read(f"{gds_input_file}")
            # Iterate through all layers in the original layout
            for layer_index in layout.layer_indices():
                # Iterate through all cells in the layout
                for cell in layout.each_cell():
                    shapes = cell.shapes(layer_index)
                    # Create a region from all shapes in the layer
                    region = pya.Region(shapes)
                    # Merge all polygons in the region
                    region.merge()
                    # Clear existing shapes and insert merged region
                    shapes.clear()
                    shapes.insert(region)
            # Commit the transaction (finalize the operation)
            layers_in_design = [layout.get_info(layer_index) for layer_index in layout.layer_indices()]
            if inplace:
                layout.write(f"{gds_input_file}")
                print(f"\nSuccessfully saved the merged GDS file inplace as {gds_input_file}.\n Merged layers: {list(layers_in_design)}")
            else:
                layout.write(f"{gds_input_file}_merged")
                print(f"\nSuccessfully saved the merged GDS file as {gds_input_file}_merged.\n Merged layers: {list(layers_in_design)}")
        except ImportError:
            print("KLayout (pya) is not installed. Please install it with 'pip install klayout' to run this function.")
        except AttributeError as e:
            print(f"Attribute error: {e}. Please ensure you're using the correct KLayout API (i.e. NOT 'pip install pya' - this is a different library). Use 'pip install klayout' instead.")

    def perforate_layer(self, layer_ind=0, x_space=20e-6, y_space=20e-6, hole_size=2e-6, keep_hole_layer=False, **kwargs):
        #TODO: Add options to perhaps output into different layer?
        #TODO: Add ability to avoid other layers or rather the edges of the layer (i.e. don't tread close to edge of ground plane etc...)

        x_space /= self.gds_units
        y_space /= self.gds_units
        dx0 = kwargs.get('x_space_edge_proportion', 0.5) * x_space
        dy0 = kwargs.get('x_space_edge_proportion', 0.5) * y_space

        xMin = (self.cx-self.sx/2) * self.unit_conv
        xMax = (self.cx+self.sx/2) * self.unit_conv
        yMin = (self.cy-self.sy/2) * self.unit_conv
        yMax = (self.cy+self.sy/2) * self.unit_conv

        x_pos = np.arange(xMin+dx0, xMax-dx0+x_space/2, x_space)
        y_pos = np.arange(yMin+dy0, yMax-dy0+y_space/2, y_space)
        hole_size /= self.gds_units

        gds_metal = None
        hole_layer = int(np.max([x for x in self._layer_metals]) + 10)

        cur_coords = [[x_pos[0]-hole_size/2, y_pos[0]-hole_size/2], [x_pos[0]+hole_size/2, y_pos[0]-hole_size/2], [x_pos[0]+hole_size/2, y_pos[0]+hole_size/2], [x_pos[0]-hole_size/2, y_pos[0]+hole_size/2]]
        seed_hole = gdstk.Polygon(cur_coords, layer=hole_layer)

        temp_cell = self.lib.new_cell(f'Holes{hole_layer}')
        temp_cell.add(seed_hole)
        # Create a reference to the cell
        ref = gdstk.Reference(temp_cell)
        # Define a rectangular repetition for the reference
        ref.repetition = gdstk.Repetition(columns=x_pos.size, rows=y_pos.size, spacing=(x_space, y_space))
        # ref.repetition = gdstk.Repetition(columns=5, rows=3, spacing=(x_space, y_space))
        temp_cell.flatten()

        holed_metal = gdstk.boolean(self._layer_metals[layer_ind], ref, 'not', layer=layer_ind, precision=self.precision)
        #Delete original layer
        self.cell.filter([(int(layer_ind), 0)], True)
        self.cell.add(*holed_metal)

        if not keep_hole_layer:
            self.cell.filter([(int(hole_layer), 0)], True)
            self.lib.remove(temp_cell)
            return None
        else:
            self._layer_metals[hole_layer] = ref
        return hole_layer

    def export(self, file_name, export_type=None, export_layers=None):
        if export_type is not None:
            self.export_type = export_type
            self.export_layers = None
            # remove any user layers
            if export_type in ["positive", "negative"]:
                self.refresh(full_refresh=True)
            # keep user layers
            elif export_type == "all":
                self.refresh()
        if export_layers is not None:
            self.export_layers = export_layers
            self.export_type = "all"
            self.refresh()
        self.lib.write_gds(file_name)
        abs_path = os.path.abspath(file_name)
        print(f"\nGDS exported at {abs_path}")
