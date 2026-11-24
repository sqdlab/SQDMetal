# Copyright 2026 Prasanna Pakkiam, Zachary Degnan
# SPDX-License-Identifier: Apache-2.0

import gdstk # type: ignore
from SQDMetal.Utilities.QUtilities import QUtilities
from time import gmtime, strftime


class ManipulateGDS:
    def __init__(
        self, import_path_gds, export_path_gds, import_cells=None, origin=(0, 0)
    ):
        """
        Class that will load a GDS file and perform actions such as making an array of the passed design (or a subsection of)

        Inputs:
            import_path_gds - Path to GDS file for import
            export_path_gds - Path to export GDS
            import_cells - (Optional) Cells to import (must be cells present in the passed GDS), defaults to all top-level cells
            origin - (Defaults to (0,0)) Tuple containing relative origin (in um) to import the design
        """

        self.infile = import_path_gds
        self.outfile = export_path_gds
        self.import_cells = [import_cells] if isinstance(import_cells, str) else import_cells
        self.origin = origin
        self.gds_units = 1e-6

        # Working library (all coordinates in um)
        self.lib = gdstk.Library("GDS Export", unit=self.gds_units)
        self.cell_flattened = None
        self.cell_array = None
        self.cell_current = None  # cell that remap/merge_overlaps/remove/report/preview act on

        # Read input, rescaled to this library's units so all coordinates are consistent
        self.gds_in = gdstk.read_gds(import_path_gds, unit=self.gds_units)
        self._cells_in = {c.name: c for c in self.gds_in.cells}

        print(f"Available cells:\n {list(self._cells_in.keys())}\n")
        print(f"Top-level cells:\n {[c.name for c in self.gds_in.top_level()]}\n")
        if self.import_cells is None:
            total_polys = sum(len(c.get_polygons()) for c in self.gds_in.top_level())
            print(f"Imported: All ({total_polys} polygons)\n")
        else:
            self._check_cells(self.import_cells, self._cells_in)
            print(f"Imported:\n {self.import_cells}\n")

    @staticmethod
    def _check_cells(cell_keys, available):
        missing = [k for k in cell_keys if k not in available]
        if missing:
            raise KeyError(
                f"Cell(s) {missing} not found. Available cells: {list(available.keys())}"
            )

    def _reset_lib(self):
        """Remove all cells from the working library (lib.cells returns a copy, so .clear() has no effect)."""
        cells = self.lib.cells
        if cells:
            self.lib.remove(*cells)

    def _working_cell(self):
        """Return the current working cell, flattening the imported design first if needed."""
        if self.cell_current is None:
            self.flatten_cells(cell_keys=self.import_cells)
        return self.cell_current

    @staticmethod
    def _layer_key(spec):
        """Normalise a layer spec: int -> (layer, None) meaning any datatype; tuple -> (layer, datatype)."""
        if isinstance(spec, int):
            return (spec, None)
        layer, datatype = spec
        return (int(layer), int(datatype))

    @classmethod
    def _matches(cls, shape, specs):
        # polygons have .datatype, labels have .texttype
        dt = shape.datatype if hasattr(shape, "datatype") else shape.texttype
        for layer, datatype in specs:
            if shape.layer == layer and (datatype is None or dt == datatype):
                return True
        return False

    def _write_cells(self, cells, path):
        """Write only the given cells (plus any cells they reference) to a GDS file."""
        out = gdstk.Library(self.lib.name, unit=self.lib.unit, precision=self.lib.precision)
        to_write = list(cells)
        for c in cells:
            to_write.extend(d for d in c.dependencies(True) if d not in to_write)
        out.add(*to_write)
        out.write_gds(path)

    # ---------------------------------------------------------------- actions

    def flatten_cells(self, cell_keys=None, export=False, export_path=None):
        """
        Takes `cell_keys` (or all top-level cells by default) from the imported GDS and flattens them,
        including any sub-cell references, into a single cell. Original layers/datatypes are preserved.

        Inputs:
            - cell_keys - (Optional) Cell name or list of cell names to extract from the passed GDS file
            - export - (Defaults to False) Choose whether to export the generated file, or just to operate on the design (default)
            - export_path - (Optional) Export path (if different from self)

        Output:
            - cell_flattened - gdstk cell containing the flattened design
        """
        self._reset_lib()
        self.cell_array = None

        timestamp = strftime("%Y%m%d_%H%M", gmtime())
        cell_name = f"Flattened_{timestamp}"
        self.cell_flattened = self.lib.new_cell(cell_name)
        self.cell_current = self.cell_flattened

        if cell_keys is None:
            # top-level only: flattening every cell would duplicate sub-cell geometry
            cell_keys = [c.name for c in self.gds_in.top_level()]
        elif isinstance(cell_keys, str):
            print(f"Single cell: {cell_keys}")
            cell_keys = [cell_keys]
        else:
            cell_keys = list(cell_keys)
        self._check_cells(cell_keys, self._cells_in)

        for key in cell_keys:
            src = self._cells_in[key]
            # get_polygons() returns flattened copies (sub-cell references included)
            polys = src.get_polygons()
            labels = src.get_labels()
            for p in polys:
                p.translate(self.origin)
            for lbl in labels:
                lbl.origin = (lbl.origin[0] + self.origin[0], lbl.origin[1] + self.origin[1])
            self.cell_flattened.add(*polys, *labels)
            layers = sorted({(p.layer, p.datatype) for p in polys})
            print(f"Copied {key:<13} ({len(polys)} polygons, layers/datatypes {layers})")

        if export:
            path = self.outfile if export_path is None else export_path
            self._write_cells([self.cell_flattened], path)
            print(
                f"\nOutputting cell {cell_name} ({len(self.cell_flattened.polygons)} polygons)\n--> {path}\n"
            )
        else:
            print(
                f"\n{cell_keys} were flattened to a single cell called '{cell_name}' with {len(self.cell_flattened.polygons)} polygons.\n"
            )

        return self.cell_flattened

    def merge_gds(
        self,
        merge_paths,
        offset=(0, 0),
        merge_cells=None,
        export=True,
        export_path=None,
    ):
        """
        Merges one or more GDS files into the current working cell. Each file is flattened and its
        polygons are added on their original layers/datatypes, overlaid at the origin (or shifted by `offset`).

        If no working cell exists yet, the imported design is flattened first (using `import_cells` from init).

        Inputs:
            - merge_paths - Path (or list of paths) to GDS files to merge in
            - offset - (Defaults to (0, 0)) (x, y) shift applied to the merged geometry. Either numbers in um,
                       or strings with units (e.g. ("1mm", "-200um")). A list of offsets (one per file) can also be passed.
            - merge_cells - (Optional) Cell name(s) to take from the merged file(s); defaults to all top-level cells
            - export - (Defaults to True) Export the merged design
            - export_path - (Optional) Export path (if different from self)

        Output:
            - cell_flattened - gdstk cell containing the merged design
        """
        if isinstance(merge_paths, str):
            merge_paths = [merge_paths]
        if isinstance(merge_cells, str):
            merge_cells = [merge_cells]

        # one offset for all files, or one per file
        if isinstance(offset, list):
            assert len(offset) == len(merge_paths), "If passing a list of offsets, provide one per merged file."
            offsets = offset
        else:
            offsets = [offset] * len(merge_paths)

        if self.cell_flattened is None:
            self.flatten_cells(cell_keys=self.import_cells)

        for path, off in zip(merge_paths, offsets):
            off_um = tuple(
                QUtilities.parse_value_length(v) / self.gds_units if isinstance(v, str) else v
                for v in off
            )
            other = gdstk.read_gds(path, unit=self.gds_units)
            other_cells = {c.name: c for c in other.cells}
            keys = (
                [c.name for c in other.top_level()] if merge_cells is None else merge_cells
            )
            self._check_cells(keys, other_cells)

            n_added = 0
            for key in keys:
                polys = other_cells[key].get_polygons()
                labels = other_cells[key].get_labels()
                for p in polys:
                    p.translate(off_um)
                for lbl in labels:
                    lbl.origin = (lbl.origin[0] + off_um[0], lbl.origin[1] + off_um[1])
                self.cell_flattened.add(*polys, *labels)
                n_added += len(polys)
            print(f"Merged {path} (cells {keys}, {n_added} polygons, offset {off_um} um)")

        if export:
            out_path = self.outfile if export_path is None else export_path
            self._write_cells([self.cell_flattened], out_path)
            print(
                f"\nOutputting cell '{self.cell_flattened.name}' ({len(self.cell_flattened.polygons)} polygons)\n--> {out_path}\n"
            )

        return self.cell_flattened

    def make_array_onChip(
        self,
        columns,
        rows,
        spacing=("50um", "50um"),
        chip_dimension=("20mm", "20mm"),
        export=True,
        export_path=None,
        use_cells=None,
        add_labels=False,
        label_string=None,
        label_offset=None,
        label_size=50,
        label_layer=100,
    ):
        """
        Makes an array of a design input and exports it as a new GDS. Optionally labels each arrayed structure.

        Inputs:
            - columns - number of columns in the array
            - rows - number of rows in the array
            - spacing - (Defaults to ("50um", "50um")) Tuple containing x- and y-pitch of the array (origin-to-origin) as strings with units
            - chip_dimension - (Defaults to ("20mm", "20mm")) Dimension of the chip on which to array the structure (WIP)
            - export - (Defaults to True) Choose whether to export the generated file (default), or just to operate on the design
            - export_path - (Optional) Export path (if different from self)
            - use_cells - (Optional) List of cells from the input GDS to include in the arrayed structure
            - add_labels - (Defaults to False) Add a label beneath each arrayed structure
            - label_string - (Optional) If None, defaults to R1C1, R1C2 etc. (R{row_index}C{column_index}). If a single string is passed, it is printed with an incrementing number. If a list (length rows * columns) is passed, labels go along row 1 from left to right, then row 2 etc.
            - label_offset - (Optional) Distance beneath each array origin for the label. Defaults to half of the y-spacing. Can be passed as a string with units (e.g. "300um").
            - label_size - (Defaults to 50) Text height in um
            - label_layer - (Defaults to 100) GDS layer for the labels

        Output:
            - cell_array - gdstk cell containing the arrayed structures
        """
        assert isinstance(chip_dimension, tuple) and isinstance(
            chip_dimension[0], str
        ), r"Input argument 'chip_dimension' should be a tuple contining strings of the chip's (x, y) dimensions with units - e.g. 'chip_dimension=('20mm', '20mm')', as in qiskit and SQDMetal."
        if isinstance(label_string, list):
            assert len(label_string) == (
                columns * rows
            ), "If you are passing a list of labels, please ensure it is the same length as the total number of arrayed structures (i.e. rows * columns)."

        sp_x = QUtilities.parse_value_length(spacing[0]) / self.gds_units
        sp_y = QUtilities.parse_value_length(spacing[1]) / self.gds_units
        lbl_offset = (
            sp_y / 2
            if label_offset is None
            else QUtilities.parse_value_length(label_offset) / self.gds_units
        )

        # use_cells re-flattens from the input file; otherwise reuse the working design
        # (so merges/remaps/etc. done beforehand are kept), flattening import_cells if needed
        if use_cells is not None:
            self.flatten_cells(cell_keys=use_cells)
        elif self.cell_flattened is None:
            self.flatten_cells(cell_keys=self.import_cells)

        # Array cell references the flattened cell (previously it referenced itself)
        if self.cell_array is not None and self.cell_array in self.lib.cells:
            self.lib.remove(self.cell_array)
        self.cell_array = self.lib.new_cell("Array")
        self.cell_array.add(
            gdstk.Reference(
                self.cell_flattened, columns=columns, rows=rows, spacing=(sp_x, sp_y)
            )
        )
        self.cell_array.flatten()
        self.cell_current = self.cell_array

        # Labels: row 1 is the top row, columns go left to right
        if add_labels:
            idx = 0
            for r in range(rows):
                y = (rows - 1 - r) * sp_y
                for c in range(columns):
                    x = c * sp_x
                    if label_string is None:
                        txt = f"R{r + 1}C{c + 1}"
                    elif isinstance(label_string, list):
                        txt = str(label_string[idx])
                    else:
                        txt = f"{label_string}{idx + 1}"
                    self.cell_array.add(
                        *gdstk.text(
                            txt,
                            label_size,
                            (x + self.origin[0], y + self.origin[1] - lbl_offset),
                            layer=label_layer,
                        )
                    )
                    idx += 1

        array_export = self.outfile if export_path is None else export_path
        if export:
            self._write_cells([self.cell_array], array_export)
            print(
                f"\nOutputting cell 'Array' ({len(self.cell_array.polygons)} polygons)\n --> {array_export}\n"
            )

        print(
            f"\nFrom {self.infile}:\n A {columns} x {rows} array of the design was made\n Spacing: {spacing}\n"
            + (f"\nThe array was exported to: {array_export}" if export else "")
        )

        return self.cell_array

    def remap_layers(self, layer_map, export=False, export_path=None):
        """
        Moves shapes from one layer/datatype to another in the current working cell.

        Inputs:
            - layer_map - dict of {source: target}. Each key/value is either a layer number (int) or a
                          (layer, datatype) tuple. An int source matches all datatypes on that layer; an int
                          target keeps each shape's datatype. E.g. {1: 10, (2, 0): (20, 5)}
            - export - (Defaults to False) Export the result
            - export_path - (Optional) Export path (if different from self)

        Output:
            - counts - dict of {source: number of shapes moved}
        """
        cell = self._working_cell()
        # resolve every shape against the ORIGINAL layers first, so swaps like {1: 2, 2: 1} work
        moves = []
        counts = {}
        for src, dst in layer_map.items():
            src_key = self._layer_key(src)
            hits = [s for s in list(cell.polygons) + list(cell.labels) if self._matches(s, [src_key])]
            moves.append((hits, dst))
            counts[src] = len(hits)
        for hits, dst in moves:
            for shape in hits:
                if isinstance(dst, int):
                    shape.layer = dst
                else:
                    shape.layer = int(dst[0])
                    if hasattr(shape, "datatype"):
                        shape.datatype = int(dst[1])
                    else:
                        shape.texttype = int(dst[1])
        for src, dst in layer_map.items():
            print(f"Remapped {src} -> {dst} ({counts[src]} shapes)")

        if export:
            self.export(export_path)
        return counts

    def remove_layers(self, layers, export=False, export_path=None):
        """
        Deletes all shapes on the given layers from the current working cell.

        Inputs:
            - layers - layer number (int), (layer, datatype) tuple, or a list of these.
                       An int removes every datatype on that layer.
            - export - (Defaults to False) Export the result
            - export_path - (Optional) Export path (if different from self)

        Output:
            - n_removed - number of shapes removed
        """
        cell = self._working_cell()
        if isinstance(layers, (int, tuple)):
            layers = [layers]
        specs = [self._layer_key(l) for l in layers]

        polys = [p for p in cell.polygons if self._matches(p, specs)]
        labels = [l for l in cell.labels if self._matches(l, specs)]
        if polys or labels:
            cell.remove(*polys, *labels)
        n_removed = len(polys) + len(labels)
        print(f"Removed {n_removed} shapes on layers {layers}")

        if export:
            self.export(export_path)
        return n_removed

    def merge_overlaps(self, layers=None, precision=1e-3, export=False, export_path=None):
        """
        Fuses overlapping/touching polygons into single polygons, separately for each layer/datatype
        (shapes on different layers are never merged together).

        Inputs:
            - layers - (Optional) layer number, (layer, datatype) tuple, or list of these to merge. Defaults to all layers.
            - precision - (Defaults to 1e-3) Boolean operation precision in um
            - export - (Defaults to False) Export the result
            - export_path - (Optional) Export path (if different from self)

        Output:
            - counts - dict of {(layer, datatype): (polygons before, polygons after)}
        """
        cell = self._working_cell()
        if layers is not None and isinstance(layers, (int, tuple)):
            layers = [layers]
        specs = None if layers is None else [self._layer_key(l) for l in layers]

        groups = {}
        for p in cell.polygons:
            if specs is None or self._matches(p, specs):
                groups.setdefault((p.layer, p.datatype), []).append(p)

        counts = {}
        for (layer, datatype), polys in sorted(groups.items()):
            merged = gdstk.boolean(
                polys, [], "or", precision=precision, layer=layer, datatype=datatype
            )
            cell.remove(*polys)
            cell.add(*merged)
            counts[(layer, datatype)] = (len(polys), len(merged))
            print(f"Layer {layer}/{datatype}: {len(polys)} -> {len(merged)} polygons")

        if export:
            self.export(export_path)
        return counts

    def report(self, verbose=True):
        """
        Summarises the current working cell: overall bounding box and size, plus polygon count,
        vertex count, area and bounding box for every layer/datatype.

        Inputs:
            - verbose - (Defaults to True) Print the report

        Output:
            - info - dict with keys 'cell', 'bounding_box', 'size', 'n_polygons', 'n_labels', 'layers'
                     ('layers' maps (layer, datatype) to a dict of polygons, vertices, area_um2, bounding_box)
        """
        cell = self._working_cell()
        bbox = cell.bounding_box()

        layers = {}
        for p in cell.polygons:
            key = (p.layer, p.datatype)
            d = layers.setdefault(key, {"polygons": 0, "vertices": 0, "area_um2": 0.0, "_boxes": []})
            d["polygons"] += 1
            d["vertices"] += p.size
            d["area_um2"] += p.area()  # note: overlapping polygons are double counted
            d["_boxes"].append(p.bounding_box())
        for d in layers.values():
            boxes = d.pop("_boxes")
            d["bounding_box"] = (
                (min(b[0][0] for b in boxes), min(b[0][1] for b in boxes)),
                (max(b[1][0] for b in boxes), max(b[1][1] for b in boxes)),
            )

        size = None if bbox is None else (bbox[1][0] - bbox[0][0], bbox[1][1] - bbox[0][1])
        info = {
            "cell": cell.name,
            "bounding_box": bbox,
            "size": size,
            "n_polygons": len(cell.polygons),
            "n_labels": len(cell.labels),
            "layers": dict(sorted(layers.items())),
        }

        if verbose:
            print(f"Cell '{cell.name}'")
            if bbox is None:
                print("  (empty)")
            else:
                print(f"  Bounding box: ({bbox[0][0]:.3f}, {bbox[0][1]:.3f}) to ({bbox[1][0]:.3f}, {bbox[1][1]:.3f}) um")
                print(f"  Size: {size[0]:.3f} x {size[1]:.3f} um")
            print(f"  {info['n_polygons']} polygons, {info['n_labels']} labels\n")
            print(f"  {'Layer/DT':<10}{'Polygons':>10}{'Vertices':>11}{'Area (um^2)':>16}")
            for (layer, dt), d in info["layers"].items():
                print(f"  {f'{layer}/{dt}':<10}{d['polygons']:>10}{d['vertices']:>11}{d['area_um2']:>16.3f}")
            print("  (area sums polygon areas; run merge_overlaps first for exact coverage)\n")
        return info

    def preview(self, path=None, scaling=1, background="#222222"):
        """
        Saves an SVG preview of the current working cell (open it in any browser).

        Inputs:
            - path - (Optional) Output .svg path. Defaults to the export path with a .svg extension.
            - scaling - (Defaults to 1) SVG pixels per um. Lower this for large chips.
            - background - (Defaults to "#222222") Background colour, or None for transparent

        Output:
            - path - the path the SVG was written to
        """
        cell = self._working_cell()
        if path is None:
            base = str(self.outfile)
            path = (base[:-4] if base.lower().endswith(".gds") else base) + ".svg"
        cell.write_svg(path, scaling=scaling, background=background if background else "none")
        print(f"Preview of '{cell.name}' --> {path}")
        return path

    def export(self, export_path=None):
        """
        Writes the current working cell to GDS.

        Inputs:
            - export_path - (Optional) Export path (if different from self)
        """
        cell = self._working_cell()
        path = self.outfile if export_path is None else export_path
        self._write_cells([cell], path)
        print(f"\nOutputting cell '{cell.name}' ({len(cell.polygons)} polygons)\n--> {path}\n")
        return path