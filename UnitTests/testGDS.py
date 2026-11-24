from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
from qiskit_metal.qlibrary.qubits.transmon_pocket import TransmonPocket

from SQDMetal.Utilities.MakeGDS import MakeGDS
from SQDMetal.Utilities.ManipulateGDS import ManipulateGDS

import contextlib
import importlib.util
import io
import os
import shutil
import unittest

import gdstk
import shapely

HAS_KLAYOUT = importlib.util.find_spec("klayout") is not None

# Test outputs are kept for visual inspection, next to this file (independent of the working directory)
OUTPUT_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gds_test_outputs")

# Chip used by the Qiskit Metal design below (um)
CHIP_X_UM, CHIP_Y_UM = 2800.0, 2000.0


# ----------------------------------------------------------------- helpers

def build_test_design():
    """Transmon + meandered readout terminated in an open, on a 2.8 mm x 2 mm chip."""
    design = designs.DesignPlanar({}, True)
    design.chips.main.size.center_x = '0.5mm'
    design.chips.main.size.center_y = '0.1mm'
    design.chips.main.size['size_x'] = '2.8mm'
    design.chips.main.size['size_y'] = '2mm'

    TransmonPocket(design, 'Q1', options=dict(
        pad_width='425 um',
        pocket_height='650um',
        connection_pads=dict(
            readout=dict(loc_W=+1, loc_H=+1, pad_width='200um')
        )))

    OpenToGround(design, 'open_to_ground', options=dict(pos_x='1.75mm', pos_y='0um', orientation='0'))
    RouteMeander(design, 'readout', Dict(
        total_length='6 mm',
        hfss_wire_bonds=True,
        fillet='90 um',
        lead=dict(start_straight='100um'),
        pin_inputs=Dict(
            start_pin=Dict(component='Q1', pin='readout'),
            end_pin=Dict(component='open_to_ground', pin='open')), ))
    design.rebuild()
    return design


def read_top_cell(path):
    lib = gdstk.read_gds(path, unit=1e-6)
    tops = lib.top_level()
    assert len(tops) == 1, f"Expected one top cell in {path}, found {[c.name for c in tops]}"
    return tops[0]


def layers_in(cell):
    return sorted({p.layer for p in cell.get_polygons()})


def layer_area(cell, layer):
    return sum(p.area() for p in cell.get_polygons(layer=layer, datatype=0))


def quiet():
    """Silence the (verbose) print output of the classes under test."""
    return contextlib.redirect_stdout(io.StringIO())


class _OutputDirTest(unittest.TestCase):
    """
    Writes each test's files to UnitTests/gds_test_outputs/<TestClass>/<test_name>/ so they can be opened
    and checked by eye afterwards (add gds_test_outputs/ to .gitignore). Each test's folder is cleared when
    that test starts, so its contents always come from the latest run.
    """

    def setUp(self):
        self._folder_path = os.path.join(OUTPUT_ROOT, type(self).__name__, self._testMethodName)
        shutil.rmtree(self._folder_path, ignore_errors=True)
        os.makedirs(self._folder_path)

    def path(self, name):
        return os.path.join(self._folder_path, name)


# ------------------------------------------------------------------ MakeGDS

class TestGDS(_OutputDirTest):
    ERR_TOL = 5e-13

    @classmethod
    def setUpClass(cls):
        # building the design is the slow part, so do it once for all MakeGDS tests
        with quiet():
            cls.design = build_test_design()

    def make_gds(self, **kwargs):
        with quiet():
            return MakeGDS(self.design, **kwargs)

    def save_result(self, gds_out, name="result.gds"):
        """Write the MakeGDS object's current state (including in-memory edits) for visual inspection."""
        path = self.path(name)
        with quiet():
            gds_out.lib.write_gds(path)
        return path

    def export(self, gds_out, name, **kwargs):
        path = self.path(name)
        with quiet():
            gds_out.export(path, **kwargs)
        self.assertTrue(os.path.isfile(path), f"{path} was not written")
        return read_top_cell(path)

    def test_GDSexport(self):
        gds_out = self.make_gds()
        for kwargs in [dict(export_type="positive"), dict(export_type="negative"),
                       dict(export_type="all"), dict(export_layers=[0, 1]), dict()]:
            with self.subTest(**kwargs):
                self.export(gds_out, "unit_test.gds", **kwargs)

    def test_export_type_layers(self):
        gds_out = self.make_gds()
        self.assertEqual(layers_in(self.export(gds_out, "pos.gds", export_type="positive")), [0])
        self.assertEqual(layers_in(self.export(gds_out, "neg.gds", export_type="negative")), [0])
        # ground (0), metal (1) and their negatives (offset by highest layer + 10)
        self.assertEqual(layers_in(self.export(gds_out, "all.gds", export_type="all")), [0, 1, 11, 12])

    def test_positive_and_negative_tile_the_chip(self):
        gds_out = self.make_gds()
        pos = layer_area(self.export(gds_out, "pos.gds", export_type="positive"), 0)
        neg = layer_area(self.export(gds_out, "neg.gds", export_type="negative"), 0)
        self.assertGreater(pos, 0)
        self.assertGreater(neg, 0)
        # metal + etched regions should cover the whole chip exactly once
        self.assertAlmostEqual(pos + neg, CHIP_X_UM * CHIP_Y_UM, delta=1.0)

    def test_export_layers_subset(self):
        gds_out = self.make_gds()
        self.assertEqual(layers_in(self.export(gds_out, "sub.gds", export_layers=[0, 1])), [0, 1])
        self.assertEqual(layers_in(self.export(gds_out, "sub1.gds", export_layers=[1])), [1])

    def test_get_cell_layers(self):
        gds_out = self.make_gds()
        self.assertEqual(gds_out.get_cell_layers(), [0, 1, 11, 12])
        self.save_result(gds_out)

    def test_delete_layers(self):
        gds_out = self.make_gds()
        gds_out.delete_layers([11, 12])
        self.assertEqual(gds_out.get_cell_layers(), [0, 1])
        self.save_result(gds_out, "after_delete_11_12.gds")
        gds_out.delete_layers(1)  # single int is accepted too
        self.assertEqual(gds_out.get_cell_layers(), [0])
        self.save_result(gds_out, "after_delete_1.gds")

    def test_add_boolean_layer(self):
        gds_out = self.make_gds()
        area_1 = sum(p.area() for p in gds_out.cell.polygons if p.layer == 1)

        # default output layer is one above the current highest layer
        out = gds_out.add_boolean_layer(1, 1, "and")
        self.assertEqual(out, 13)
        area_and = sum(p.area() for p in gds_out.cell.polygons if p.layer == out)
        self.assertAlmostEqual(area_and, area_1, delta=1.0)

        # metal and its own negative should not overlap (allow sub-um^2 slivers from curved-edge rounding)
        with quiet():
            out = gds_out.add_boolean_layer(1, 12, "and", output_layer=50)
        self.assertLess(sum(p.area() for p in gds_out.cell.polygons if p.layer == out), 1.0)

        # user layers survive an "all" export
        cell = self.export(gds_out, "bool.gds", export_type="all")
        self.assertIn(13, layers_in(cell))

    def test_add_boolean_layer_rejects_bad_operation(self):
        gds_out = self.make_gds()
        with self.assertRaises(AssertionError):
            gds_out.add_boolean_layer(0, 1, "nand")

    def test_add_text_on_new_layer(self):
        gds_out = self.make_gds()
        gds_out.add_text("SQDMetal", layer=5, size=100, position=(0.1, 0.1))
        self.assertIn(5, gds_out.get_cell_layers())
        cell = self.export(gds_out, "text.gds", export_type="all")
        self.assertGreater(layer_area(cell, 5), 0)

    def test_add_text_cut_from_ground_plane(self):
        gds_out = self.make_gds(export_type="positive")
        before = sum(p.area() for p in gds_out.cell.polygons if p.layer == 0)
        gds_out.add_text("SQDMetal", layer=0, size=100, position=(0.1, 0.1))
        after = sum(p.area() for p in gds_out.cell.polygons if p.layer == 0)
        self.assertLess(after, before)
        self.save_result(gds_out)

    def test_add_text_rejects_bad_position(self):
        gds_out = self.make_gds()
        with self.assertRaises(AssertionError):
            gds_out.add_text("x", layer=5, position=(1.5, 0.5))

    def test_perforate_layer(self):
        gds_out = self.make_gds()
        before = sum(p.area() for p in gds_out.cell.polygons if p.layer == 0)
        result = gds_out.perforate_layer(0, x_space=100e-6, y_space=100e-6, hole_size=5e-6)
        self.assertIsNone(result)  # hole stencil is removed by default
        after = sum(p.area() for p in gds_out.cell.polygons if p.layer == 0)
        self.assertLess(after, before)
        # each hole removes at most 5x5 um^2, and the grid has at most 28 x 20 holes
        self.assertLessEqual(before - after, 28 * 20 * 25 + 1)
        self.save_result(gds_out)

    def test_perforate_layer_updates_negative_layer(self):
        # holes cut from the ground plane (0) must appear as filled squares in its negative layer (11)
        gds_out = self.make_gds(export_type="all")
        area = lambda layer: sum(p.area() for p in gds_out.cell.polygons if p.layer == layer)
        gnd_before, neg_before = area(0), area(11)
        gds_out.perforate_layer(0, x_space=100e-6, y_space=100e-6, hole_size=5e-6)
        removed = gnd_before - area(0)
        self.assertGreater(removed, 0)
        self.assertAlmostEqual(area(11) - neg_before, removed, delta=1.0)
        self.assertAlmostEqual(area(0) + area(11), CHIP_X_UM * CHIP_Y_UM, delta=1.0)
        self.save_result(gds_out)

    def test_perforate_layer_negative_export(self):
        # in a negative export, layer 0 is the etched pattern, so perforating adds the holes to it
        gds_out = self.make_gds(export_type="negative")
        before = sum(p.area() for p in gds_out.cell.polygons if p.layer == 0)
        gds_out.perforate_layer(0, x_space=100e-6, y_space=100e-6, hole_size=5e-6)
        self.assertGreater(sum(p.area() for p in gds_out.cell.polygons if p.layer == 0), before)
        self.save_result(gds_out)

    @staticmethod
    def _hole_polygons(gds_out, hole_layer):
        holder = gdstk.Cell("holes")
        holder.add(gds_out._layer_metals[hole_layer])
        return [shapely.Polygon(p.points) for p in holder.get_polygons()]

    @staticmethod
    def _to_shapely(polys):
        return shapely.union_all([shapely.Polygon(p.points).buffer(0) for p in polys])

    def test_perforate_respects_clearance(self):
        clearance_um = 10
        gds_out = self.make_gds(export_type="all")
        ground = self._to_shapely([p for p in gds_out.cell.polygons if p.layer == 0])
        metal_1 = self._to_shapely([p for p in gds_out.cell.polygons if p.layer == 1])
        hole_layer = gds_out.perforate_layer(0, x_space=50e-6, y_space=50e-6, hole_size=5e-6,
                                             clearance=clearance_um * 1e-6, keep_hole_layer=True)
        holes = self._hole_polygons(gds_out, hole_layer)
        self.assertGreater(len(holes), 100)
        all_holes = shapely.union_all(holes)
        self.assertAlmostEqual(all_holes.difference(ground).area, 0, delta=1e-6)  # every hole is inside the metal
        self.assertGreaterEqual(min(ground.boundary.distance(h) for h in holes), clearance_um - 1e-6)
        self.assertGreaterEqual(metal_1.distance(all_holes), clearance_um - 1e-6)
        self.save_result(gds_out)

    def test_perforate_avoids_other_layers(self):
        # text on its own layer sits on top of the ground plane; holes must keep clear of it
        gds_out = self.make_gds(export_type="all")
        gds_out.add_text("KEEP OUT", layer=5, size=150, position=(0.05, 0.05))
        text = self._to_shapely([p for p in gds_out.cell.polygons if p.layer == 5])
        hole_layer = gds_out.perforate_layer(0, x_space=20e-6, y_space=20e-6, hole_size=4e-6,
                                             clearance=5e-6, keep_hole_layer=True)
        holes = shapely.union_all(self._hole_polygons(gds_out, hole_layer))
        self.assertGreaterEqual(text.distance(holes), 5 - 1e-6)

        # avoid_layers=[] only keeps clear of the ground plane's own edges, so holes land under the text
        gds_out = self.make_gds(export_type="all")
        gds_out.add_text("KEEP OUT", layer=5, size=150, position=(0.05, 0.05))
        hole_layer = gds_out.perforate_layer(0, x_space=20e-6, y_space=20e-6, hole_size=4e-6,
                                             clearance=5e-6, keep_hole_layer=True, avoid_layers=[])
        holes = shapely.union_all(self._hole_polygons(gds_out, hole_layer))
        self.assertGreater(text.intersection(holes).area, 0)
        self.save_result(gds_out)

    def test_perforate_avoid_geometry_off_uses_full_grid(self):
        gds_out = self.make_gds(export_type="all")
        hole_layer = gds_out.perforate_layer(0, x_space=100e-6, y_space=100e-6, hole_size=5e-6,
                                             avoid_geometry=False, keep_hole_layer=True)
        self.assertEqual(len(self._hole_polygons(gds_out, hole_layer)), 28 * 20)

    def test_perforate_clearance_too_large(self):
        gds_out = self.make_gds(export_type="all")
        before = sum(p.area() for p in gds_out.cell.polygons if p.layer == 0)
        with quiet():
            result = gds_out.perforate_layer(0, x_space=100e-6, y_space=100e-6, hole_size=5e-6, clearance=5e-3)
        self.assertIsNone(result)
        self.assertEqual(sum(p.area() for p in gds_out.cell.polygons if p.layer == 0), before)  # unchanged

    def test_perforate_layer_keep_hole_layer(self):
        gds_out = self.make_gds()
        hole_layer = gds_out.perforate_layer(0, x_space=100e-6, y_space=100e-6, hole_size=5e-6, keep_hole_layer=True)
        self.assertIsInstance(hole_layer, int)
        self.assertGreater(hole_layer, 12)
        self.save_result(gds_out)

    @unittest.skipUnless(HAS_KLAYOUT, "KLayout (pip install klayout) is required")
    def test_merge_gds_per_layer_inplace(self):
        # two overlapping rectangles on layer 1 should become one polygon
        path = self.path("overlap.gds")
        lib = gdstk.Library()
        cell = lib.new_cell("Top")
        cell.add(gdstk.rectangle((0, 0), (10, 10), layer=1), gdstk.rectangle((5, 0), (15, 10), layer=1))
        lib.write_gds(path)

        with quiet():
            MakeGDS.merge_gds_per_layer_inplace(path, inplace=False)
        merged_path = self.path("overlap_merged.gds")
        self.assertTrue(os.path.isfile(merged_path))
        self.assertEqual(len(read_top_cell(merged_path).polygons), 1)
        self.assertEqual(len(read_top_cell(path).polygons), 2)  # original untouched

        # extension is optional, and inplace overwrites the input
        with quiet():
            MakeGDS.merge_gds_per_layer_inplace(path[:-4])
        merged = read_top_cell(path)
        self.assertEqual(len(merged.polygons), 1)
        self.assertAlmostEqual(merged.polygons[0].area(), 150.0, places=3)


# ------------------------------------------------------------- ManipulateGDS

class TestManipulateGDS(_OutputDirTest):

    def setUp(self):
        super().setUp()
        # hierarchical design: Top references Sub; Top also has its own shapes and a label
        lib = gdstk.Library()
        sub = lib.new_cell("Sub")
        sub.add(gdstk.rectangle((0, 0), (10, 10), layer=1))
        top = lib.new_cell("Top")
        top.add(gdstk.Reference(sub, (20, 0)))
        top.add(gdstk.rectangle((0, 0), (10, 10), layer=1),
                gdstk.rectangle((5, 0), (15, 10), layer=1),       # overlaps the one above
                gdstk.rectangle((0, 20), (5, 25), layer=2, datatype=3),
                gdstk.rectangle((0, 30), (5, 35), layer=7))
        top.add(gdstk.Label("ID", (0, 40), layer=7))
        self.input_gds = self.path("input.gds")
        lib.write_gds(self.input_gds)

        # second file saved in nm, to check unit handling on merge
        lib_nm = gdstk.Library(unit=1e-9)
        lib_nm.new_cell("Markers").add(gdstk.rectangle((0, 0), (3000, 3000), layer=5))
        self.nm_gds = self.path("markers_nm.gds")
        lib_nm.write_gds(self.nm_gds)

        self.output_gds = self.path("output.gds")

    def manip(self, **kwargs):
        with quiet():
            return ManipulateGDS(self.input_gds, self.output_gds, **kwargs)

    def run_quiet(self, func, *args, **kwargs):
        with quiet():
            return func(*args, **kwargs)

    def save_result(self, m, name="result.gds"):
        """Write the current working cell for visual inspection."""
        return self.run_quiet(m.export, self.path(name))

    def test_flatten_includes_subcells_and_labels(self):
        m = self.manip()
        cell = self.run_quiet(m.flatten_cells, export=True)
        # 4 shapes in Top + 1 from the referenced Sub cell (not double counted)
        self.assertEqual(len(cell.polygons), 5)
        self.assertEqual(len(cell.labels), 1)
        out = read_top_cell(self.output_gds)
        self.assertEqual(len(out.polygons), 5)
        bboxes = sorted(tuple(map(tuple, p.bounding_box())) for p in out.polygons if p.layer == 1)
        self.assertIn(((20.0, 0.0), (30.0, 10.0)), bboxes)  # Sub placed at its reference offset

    def test_flatten_single_cell_name_and_origin(self):
        m = self.manip(origin=(100, 50))
        cell = self.run_quiet(m.flatten_cells, "Sub")  # a plain string, not split into letters
        self.assertEqual(len(cell.polygons), 1)
        (x0, y0), _ = cell.polygons[0].bounding_box()
        self.assertAlmostEqual(x0, 100)
        self.assertAlmostEqual(y0, 50)
        self.save_result(m)

    def test_missing_cell_raises(self):
        with self.assertRaises(KeyError):
            self.manip(import_cells=["NotACell"])
        m = self.manip()
        with self.assertRaises(KeyError):
            self.run_quiet(m.flatten_cells, ["NotACell"])

    def test_merge_gds_units_and_offset(self):
        m = self.manip()
        self.run_quiet(m.merge_gds, self.nm_gds, offset=("100um", "-0.05mm"))
        out = read_top_cell(self.output_gds)
        self.assertEqual(len(out.polygons), 6)
        marker = [p for p in out.polygons if p.layer == 5]
        self.assertEqual(len(marker), 1)
        (x0, y0), (x1, y1) = marker[0].bounding_box()
        # 3000 nm square, shifted by (100, -50) um
        self.assertAlmostEqual(x0, 100, places=6)
        self.assertAlmostEqual(y0, -50, places=6)
        self.assertAlmostEqual(x1 - x0, 3, places=6)

    def test_merge_gds_multiple_offsets(self):
        m = self.manip()
        self.run_quiet(m.merge_gds, [self.nm_gds, self.nm_gds], offset=[(0, 0), (500, 0)], export=False)
        xs = sorted(p.bounding_box()[0][0] for p in m.cell_current.polygons if p.layer == 5)
        self.assertEqual(len(xs), 2)
        self.assertAlmostEqual(xs[0], 0)
        self.assertAlmostEqual(xs[1], 500)
        self.save_result(m)
        with self.assertRaises(AssertionError):
            m.merge_gds([self.nm_gds, self.nm_gds], offset=[(0, 0)])

    def test_remap_layers_with_swap(self):
        m = self.manip()
        counts = self.run_quiet(m.remap_layers, {1: 2, (2, 3): (1, 0)})
        # layer 1 holds 3 shapes: two in Top plus the flattened Sub square
        self.assertEqual(counts, {1: 3, (2, 3): 1})
        cell = m.cell_current
        self.assertEqual(sum(1 for p in cell.polygons if p.layer == 2), 3)
        moved = [p for p in cell.polygons if p.layer == 1]
        self.assertEqual(len(moved), 1)
        self.assertEqual(moved[0].datatype, 0)
        self.save_result(m)

    def test_remove_layers(self):
        m = self.manip()
        n = self.run_quiet(m.remove_layers, 7)  # polygon and label on layer 7
        self.assertEqual(n, 2)
        self.assertNotIn(7, {p.layer for p in m.cell_current.polygons})
        self.assertEqual(len(m.cell_current.labels), 0)
        # tuple spec only matches that datatype
        self.assertEqual(self.run_quiet(m.remove_layers, (2, 0)), 0)
        self.assertEqual(self.run_quiet(m.remove_layers, [(2, 3)]), 1)
        self.save_result(m)

    def test_merge_overlaps(self):
        m = self.manip()
        counts = self.run_quiet(m.merge_overlaps, 1)
        # two overlapping squares fuse; the separate Sub square stays on its own
        self.assertEqual(counts[(1, 0)], (3, 2))
        info = m.report(verbose=False)
        self.assertAlmostEqual(info["layers"][(1, 0)]["area_um2"], 150 + 100, places=3)
        self.assertEqual(info["layers"][(2, 3)]["polygons"], 1)  # other layers untouched
        self.save_result(m)

    def test_report(self):
        m = self.manip()
        info = self.run_quiet(m.report, verbose=False)
        self.assertEqual(info["n_polygons"], 5)
        self.assertEqual(info["n_labels"], 1)
        self.assertEqual(set(info["layers"]), {(1, 0), (2, 3), (7, 0)})
        self.assertEqual(info["layers"][(1, 0)]["polygons"], 3)
        (x0, y0), (x1, y1) = info["bounding_box"]
        self.assertAlmostEqual(x1 - x0, info["size"][0])
        self.assertAlmostEqual(x1, 30)
        self.save_result(m)

    def test_preview(self):
        m = self.manip()
        path = self.run_quiet(m.preview)
        self.assertEqual(path, self.output_gds[:-4] + ".svg")
        with open(path) as f:
            self.assertIn("<svg", f.read(200))

    def test_export(self):
        m = self.manip()
        custom = self.path("custom.gds")
        self.assertEqual(self.run_quiet(m.export, custom), custom)
        self.assertEqual(len(read_top_cell(custom).polygons), 5)

    def test_make_array(self):
        m = self.manip()
        arr = self.run_quiet(m.make_array_onChip, 3, 2, spacing=("100um", "200um"))
        self.assertEqual(len(arr.polygons), 5 * 3 * 2)
        out = read_top_cell(self.output_gds)
        self.assertEqual(out.name, "Array")
        (x0, y0), (x1, y1) = out.bounding_box()
        self.assertAlmostEqual(x1 - x0, 2 * 100 + 30)
        self.assertAlmostEqual(y1 - y0, 1 * 200 + 40)

    def test_make_array_labels(self):
        m = self.manip()
        arr = self.run_quiet(m.make_array_onChip, 2, 2, add_labels=True, label_layer=100, export_path=self.path("labelled_array.gds"))
        self.assertIn(100, {p.layer for p in arr.polygons})
        with self.assertRaises(AssertionError):
            m.make_array_onChip(2, 2, label_string=["only", "three", "labels"], export=False)

    def test_make_array_uses_edited_design(self):
        m = self.manip()
        self.run_quiet(m.remove_layers, 7)
        arr = self.run_quiet(m.make_array_onChip, 2, 1, export=False)
        self.assertNotIn(7, {p.layer for p in arr.polygons})
        self.assertEqual(len(arr.polygons), 4 * 2)
        self.save_result(m)


class TestMakeGDSToManipulateGDS(_OutputDirTest):
    """End-to-end: export a Qiskit Metal design, then post-process it with ManipulateGDS."""

    def test_pipeline(self):
        with quiet():
            design = build_test_design()
            gds_path = self.path("design.gds")
            MakeGDS(design).export(gds_path, export_type="all")

            m = ManipulateGDS(gds_path, self.path("processed.gds"))
            m.remove_layers([11, 12])
            m.remap_layers({1: 2})
            m.merge_gds(gds_path, offset=("3mm", "0mm"), export=False)
            m.merge_overlaps()
            info = m.report(verbose=False)
            arr = m.make_array_onChip(2, 2, spacing=("7mm", "3mm"))

        layers = {layer for layer, _ in info["layers"]}
        self.assertEqual(layers, {0, 1, 2, 11, 12})  # second (unedited) copy brings 1, 11, 12 back
        self.assertAlmostEqual(info["size"][0], CHIP_X_UM + 3000, delta=1.0)
        self.assertEqual(len(arr.polygons), 4 * info["n_polygons"])
        self.assertTrue(os.path.isfile(self.path("processed.gds")))


if __name__ == '__main__':
    unittest.main()