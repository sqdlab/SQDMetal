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
    def __init__(self, design, cell_name='Main', threshold:float=1e-9, precision:float=1e-9, curve_resolution:int=22, smooth_radius:float=0, export_type:str='all', export_layers:list[int]=[], print_statements:bool=False):
        """
        Class to export a Qiskit-Metal design into GDS.

        Parameters
        ----------
        design : QDesign
            The QiskitMetal design object
        fuse_threshold : float 
            The smallest feature crack that can exist; anything smaller will get sealed together. This is to help remove
            artefacts from floating-point inaccuracies. It is given in metres and defaults to 1e-9.
        precision : float
            The GDS granularity/precision given in metres (defaults to 1e-9).
        curve_resolution : int
            Number of vertices to use on a path fillet; that is, number of points over a quarter arc. Default value is 22.
        smooth_radius : float
            The radius to use upon applying a final postprocessing step to round off all corners into smooth bends. The
            default value is 0, in which case no rounding is performed.
        export_type : str
            Select which layers to export (default is ``'all'``):

                - ``'all'`` (exports positive and negative patterns and GND plane as seperate layers), 
                - ``'positive'`` (exports flattened metals to layer 0), or 
                - ``'negative'`` (exports negative flattened metals to layer 0)

        export_layers : list[int]
            Layers to export given as a list of integers or single integer (WIP). Default is an empty list in which case
            all layers are exported.
        print_statements : bool
            Print debug/status messages when exporting. Default value is False.
        """
        self._design = design
        self.cell_name = cell_name
        self.curve_resolution = curve_resolution
        self.threshold = threshold
        self.precision = precision
        self.smooth_radius = smooth_radius
        self.gds_units = 1e-6
        self.export_type = export_type
        self.export_layers = export_layers
        self.print_statements = print_statements
        self._refresh(full_refresh=True)

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

    def _refresh(self, full_refresh=False):
        """
        Refresh the current state of the GDS construction where:

        - Everything is deleted and metallic layers are rebuilt from scratch with the ground plane being in layer zero.
        - Negative versions of the layers are given in layers starting from highest layer index plus 10.
        - The compiled layers are are specified by the user: negative plus positive patterns, or only positive or negative patterns.
        - User-defined layers (e.g., via :func:`~SQDMetal.Utilities.MakeGDS.MakeGDS.add_boolean_layer`) are preserved unless their layer index overlaps with a standard layer being recreated.
        """
        # Preserve any existing user-defined layers
        existing_user_layers = dict(getattr(self, "_layer_metals", {}))

        # Validate export_type and export_layers before rebuilding
        assert self.export_type in ["all", "negative", "positive"], "Export type must be: \"all\", \"negative\", or \"positive\""

        if self.export_layers is not None and len(self.export_layers) > 0:
            self.export_type = "all"
            print(f"Exporting layers {self.export_layers}") if self.print_statements else 0

        # Parse geometry and update metadata before rebuilding geometry
        qmpl = QiskitShapelyRenderer(design=self._design, canvas=None, logger=None)
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
        self.cell = self.lib.new_cell(self.cell_name)
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
        self.neg_layer_offset = int(neg_layer_offset)  # negative of layer L is written to layer L + neg_layer_offset
        self._negative_layers = {int(l) + self.neg_layer_offset for l in [0, *leLayers]}  # auto-generated negatives
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
            self._merge_polygons_in_layer(0)
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
        if self.export_layers is not None and len(self.export_layers) > 0:
            assert isinstance(self.export_layers, list), "Please pass a list of integers as 'export_layers'."
            layer_list = self.get_cell_layers()
            assert set(self.export_layers).issubset(layer_list), "Chosen export layers are not present in the design."
            to_delete = list(set(layer_list) - set(self.export_layers))
            if isinstance(to_delete, list):
                for i in list(to_delete):
                    self.cell.remove(*[p for p in self.cell.polygons if p.layer == i])
            else:
                self.cell.remove(*[p for p in self.cell.polygons if p.layer == to_delete])
        
    def _merge_polygons_in_layer(self, layer):
        """
        Merge all polygons within a given layer.

        Parameters
        ----------
        layer : int
            Index of layer to have all polygons merged.
        """
        polys_to_merge = [p for p in self.cell.polygons if p.layer == layer]
        merged = gdstk.boolean(polys_to_merge, [], "or", layer=layer)
        self.cell.remove(*polys_to_merge)
        if merged:
            self.cell.add(*merged)

    def add_boolean_layer(self, layer1_ind:int, layer2_ind:int, operation:str, output_layer:int=None):
        """
        Perform a boolean operation on the polygonal components across two layers.

        Parameters
        ----------
        layer1_ind : int
            Index of first layer
        layer2_ind : int
            Index of second layer
        operation : str
            Boolean operation to perform over the two layers. It can be: ``"and"``, ``"or"``, ``"xor"``, ``"not"``
        output_layer : int
            If an integer is given here, the result of the boolean operation is written into this layer (creating one if
            it does not exist). Default value is None, in which case, the layer is automatically created by taking an
            integer one above the largest index thus far.
        
        Returns
        -------
        output_layer : int
            Integer index of the layer containing the result of the boolean operation.
        """
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
    def add_text(self, text_label:str="", layer:int=0, size:int=300, position=None):
        """
        Add a text label to the gds export on a given layer at a given position.

        Parameters
        ----------
        text_label : str
            Text label to add
        layer : int
            Layer number to export the text to (defaults to layer 0)
        size : str
            Text size in GDS Units (default value is 300)
        position : tuple
            Tuple containing position (normalised to 1); e.g. (0,0) is bottom left, (1,1) is top right.
        """
        #TODO: @6biscuits, check that it is specified in metres so that it's consistent with all the other interfaces. It looks like GDS units here...
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
        """
        Get all layers across the current cell.

        Returns
        ----------
        layers : list[int]
            List of layers present in current cell.
        """
        layers = set()
        for poly in self.cell.polygons:
            layers.add(poly.layer)
        for path in self.cell.paths:
            layers.add(path.layer)
        for label in self.cell.labels:
            layers.add(label.layer)
        return sorted(layers)
    
    def delete_layers(self, layers_to_remove:list[int]):
        """
        Delete a list of layers.

        Parameters
        ----------
        layers_to_remove : list[int]
            Layers to remove given as a list of integers.
        """
        if isinstance(layers_to_remove, int):
            layers_to_remove = [layers_to_remove]
        self.cell.filter([(layer, 0) for layer in layers_to_remove])
    
    @staticmethod
    def merge_gds_per_layer_inplace(gds_input_file:str, inplace:bool=True):
        """
        Open a GDS file, go through every layer in every cell/region and merge all polygons within, then save the GDS file.

        Parameters
        ----------
        gds_input_file : str
            Path to the GDS file
        inplace : bool
            If True (the default value), the input file is overwritten with the polygons merged within every layer. Otherwise,
            the output is written into the same location as the input file with the suffix *_merged* added to the name.
        """
        try:
            import pya # type: ignore
            gds_input_file = str(gds_input_file)
            if not gds_input_file.lower().endswith(".gds"):
                gds_input_file = f"{gds_input_file}.gds"
            layout = pya.Layout()
            layout.read(gds_input_file)
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
                output_file = gds_input_file
            else:
                base, ext = os.path.splitext(gds_input_file)
                output_file = f"{base}_merged{ext}"
            layout.write(output_file)
            print(f"\nSuccessfully saved the merged GDS file as {output_file}.\n Merged layers: {list(layers_in_design)}")
        except ImportError:
            print("KLayout (pya) is not installed. Please install it with 'pip install klayout' to run this function.")
        except AttributeError as e:
            print(f"Attribute error: {e}. Please ensure you're using the correct KLayout API (i.e. NOT 'pip install pya' - this is a different library). Use 'pip install klayout' instead.")

    def perforate_layer(self, layer_ind:int=0, x_space:float=20e-6, y_space:float=20e-6, hole_size:float=2e-6, keep_hole_layer:bool=False, x_space_edge_proportion:float=0.5, y_space_edge_proportion:float=0.5, clearance:float=2e-6, avoid_geometry:bool=True, avoid_layers:list[int]=None):
        """
        Perforate a given layer with a uniform grid of square holes (a.k.a. cheesing).

        By default, holes are only placed where they fit entirely inside the layer's metal, at least ``clearance``
        away from the metal's edges (e.g. CPW gaps, the chip boundary) and from shapes on every other layer.

        Parameters
        ----------
        layer_ind : int
            Index of layer to perforate
        x_space : float
            Spacing of the holes along the x-axis given in metres. Default value is 20e-6.
        y_space : float
            Spacing of the holes along the y-axis given in metres. Default value is 20e-6.
        hole_size : float
            Size of each square hole in metres. Default value is 2e-6
        keep_hole_layer : bool
            If True, the temporary layer created as a stencil (i.e. the hole polygons) is kept. Defaults to False.
        x_space_edge_proportion : float
            Padding on the left and right edges before placing the holes. Given as a proportion of x_space (default value being 0.5).
        y_space_edge_proportion : float
            Padding on the top and bottom edges before placing the holes. Given as a proportion of y_space (default value being 0.5).
        clearance : float
            Minimum distance in metres between the edge of any hole and (a) the edge of the perforated metal and
            (b) any geometry on the avoided layers. Default value is 2e-6 (holes may not overlap anything, but can
            come right up to it). Placement is slightly conservative near diagonal or curved edges (by up to
            ~0.2 x hole_size), as the check treats each hole as the circle enclosing it.
        avoid_geometry : bool
            If True (default), skip holes that would violate ``clearance``. If False, every grid position is used
            and holes may cut into edges and overlap other layers (the previous behaviour).
        avoid_layers : list[int]
            Layers whose shapes the holes must keep clear of. Default is None, which uses every layer currently in
            the cell except the perforated layer itself and the auto-generated negative layers (these cover the
            whole chip by construction). Pass ``[]`` to only keep clear of the perforated layer's own edges.

        Returns
        -------
        hole_layer : int
            Integer index of the hole layer if ``keep_hole_layer`` is True, otherwise None.
        """
        #TODO: Add options to perhaps output into different layer?

        x_space /= self.gds_units
        y_space /= self.gds_units
        dx0 = x_space_edge_proportion * x_space
        dy0 = y_space_edge_proportion * y_space

        xMin = (self.cx-self.sx/2) * self.unit_conv
        xMax = (self.cx+self.sx/2) * self.unit_conv
        yMin = (self.cy-self.sy/2) * self.unit_conv
        yMax = (self.cy+self.sy/2) * self.unit_conv

        x_pos = np.arange(xMin+dx0, xMax-dx0+x_space/2, x_space)
        y_pos = np.arange(yMin+dy0, yMax-dy0+y_space/2, y_space)
        hole_size /= self.gds_units
        clearance /= self.gds_units
        assert clearance >= 0, "clearance must be non-negative."

        chip_rect = gdstk.rectangle((xMin, yMin), (xMax, yMax))
        negative_pattern = (self.export_type == "negative" and layer_ind == 0)

        # hole centres on the full grid
        xx, yy = np.meshgrid(x_pos, y_pos)
        centres = np.column_stack([xx.ravel(), yy.ravel()])
        n_grid = len(centres)

        if avoid_geometry:
            def to_shapely(polys):
                geoms = [shapely.Polygon(p.points) for p in polys if len(p.points) >= 3]
                # buffer(0) repairs the self-touching "keyhole" polygons that gdstk uses for holes
                return shapely.union_all([g if g.is_valid else g.buffer(0) for g in geoms]) if geoms else shapely.Polygon()

            # metal that can be perforated (in a negative export, layer 0 is the etched pattern instead)
            layer_polys = [p for p in self.cell.polygons if p.layer == layer_ind]
            if negative_pattern:
                metal = shapely.Polygon(chip_rect.points).difference(to_shapely(layer_polys))
            else:
                metal = to_shapely(layer_polys)

            # geometry on other layers to keep clear of
            if avoid_layers is None:
                skip = {int(layer_ind)} | getattr(self, "_negative_layers", set())
                avoid_layers = [l for l in self.get_cell_layers() if l not in skip]
            keepout = to_shapely([p for p in self.cell.polygons if p.layer in set(avoid_layers)])

            # a hole centre is valid if the hole's enclosing circle (radius hole_size/sqrt(2)) plus clearance
            # fits inside the metal and stays clear of the keep-out geometry
            margin = clearance + hole_size / np.sqrt(2)
            allowed = metal.buffer(-margin)
            if not keepout.is_empty:
                allowed = allowed.difference(keepout.buffer(margin))
            shapely.prepare(allowed)
            centres = centres[shapely.contains_xy(allowed, centres[:, 0], centres[:, 1])]

            print(f"Perforating layer {layer_ind}: {len(centres)} of {n_grid} holes placed (avoiding layers {avoid_layers}, clearance {clearance:g} um)") if self.print_statements else 0

        if len(centres) == 0:
            print(f"Warning: no holes could be placed on layer {layer_ind} with the given spacing, size and clearance.")
            return None

        hole_layer = int(np.max([x for x in self._layer_metals]) + 10)

        # one seed hole, placed at the first valid centre and repeated at the others
        # (gdstk counts the reference's own position as the first copy, so offsets are relative to it)
        seed_hole = gdstk.rectangle((-hole_size/2, -hole_size/2), (hole_size/2, hole_size/2), layer=hole_layer)
        temp_cell = self.lib.new_cell(f'Holes{hole_layer}')
        temp_cell.add(seed_hole)
        ref = gdstk.Reference(temp_cell, origin=tuple(centres[0]))
        if len(centres) > 1:
            ref.repetition = gdstk.Repetition(offsets=(centres[1:] - centres[0]).tolist())

        if negative_pattern:
            # layer 0 holds the negative (etched) pattern: holes in the metal are extra etched area
            holes_in_chip = gdstk.boolean(ref, chip_rect, "and", precision=self.precision)
            holed_metal = gdstk.boolean(self._layer_metals[layer_ind], holes_in_chip, "or", layer=layer_ind, precision=self.precision)
        else:
            holed_metal = gdstk.boolean(self._layer_metals[layer_ind], ref, "not", layer=layer_ind, precision=self.precision)

        # Replace the original layer
        self.cell.filter([(int(layer_ind), 0)], True)
        self.cell.add(*holed_metal)
        self._layer_metals[layer_ind] = holed_metal

        # Update the matching negative layer (only present for "all" exports) so the holes appear there too
        neg_layer = int(layer_ind) + getattr(self, "neg_layer_offset", 0)
        if self.export_type == "all" and neg_layer != layer_ind and neg_layer in self._layer_metals:
            holed_neg = gdstk.boolean(chip_rect, holed_metal, "not", layer=neg_layer, precision=self.precision)
            self.cell.filter([(neg_layer, 0)], True)
            self.cell.add(*holed_neg)
            self._layer_metals[neg_layer] = holed_neg

        if not keep_hole_layer:
            self.lib.remove(temp_cell)
            return None
        else:
            self._layer_metals[hole_layer] = ref
        return hole_layer

    def export(self, file_name:str, export_type=None, export_layers=None):
        """
        Perform the final export to write to the GDS file.

        Parameters
        ----------
        file_name : str
            File path to the exported GDS file output.
        export_type : str
            Select which layers to export (default is None and thus, the value given in the initialiser is used):

                - ``'all'`` (exports positive and negative patterns and GND plane as seperate layers), 
                - ``'positive'`` (exports flattened metals to layer 0), or 
                - ``'negative'`` (exports negative flattened metals to layer 0)

        export_layers : list[int]
            Layers to export given as a list of integers or single integer (WIP). Default is None and thus, the
            value given in the initialiser is used.
        """
        if export_type is not None:
            self.export_type = export_type
            self.export_layers = None
            # remove any user layers
            if export_type in ["positive", "negative"]:
                self._refresh(full_refresh=True)
            # keep user layers
            elif export_type == "all":
                self._refresh()
        if export_layers is not None:
            self.export_layers = export_layers
            self.export_type = "all"
            self._refresh()
        self.lib.write_gds(file_name)
        abs_path = os.path.abspath(file_name)
        print(f"\nGDS exported at {abs_path}")