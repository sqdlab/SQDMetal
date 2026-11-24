import shutil
import tempfile
import unittest
from pathlib import Path

import nbformat
import numpy as np
import shapely
import qiskit_metal as metal
from qiskit_metal import Dict
from qiskit_metal.qlibrary.couplers.coupled_line_tee import CoupledLineTee
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
from qiskit_metal.qlibrary.tlines.straight_path import RouteStraight

from SQDMetal.Comps.Qubits import TransmonTaperedInsets
from SQDMetal.Utilities.SingleQubitCells import SingleQubitCells

CPW = dict(trace_width="10um", trace_gap="6um")


def make_chip():
    """Two qubits on a shared feedline: Q1 above it (orientation 180), Q2 below it (orientation 0)."""
    design = metal.designs.DesignPlanar(metadata=Dict(design_name="test_chip"))
    design.overwrite_enabled = True
    design.chips.main.size.size_x = "6mm"
    design.chips.main.size.size_y = "6mm"
    qopts = Dict(pad_gap="20um", taper_height="0um", pocket_width="940um", pocket_height="480um",
                 pocket_lower_tighten="100um",
                 connection_pads=Dict(readout=Dict(pad_gap="10um", cpw_width="10um", cpw_gap="6um")))
    tee = dict(coupling_space="5um", coupling_length="220um", down_length="200um", prime_width="10um",
               prime_gap="6um", second_width="10um", second_gap="6um", open_termination=True)
    for name, x, y, orient in (("Q1", -1.0, 1.5, 180), ("Q2", 1.0, -1.5, 0)):
        TransmonTaperedInsets(design, name, options=Dict(**qopts, pos_x=f"{x}mm", pos_y=f"{y}mm",
                                                         orientation=str(orient)))
        CoupledLineTee(design, f"{name}_Read_T", options=dict(**tee, pos_x=f"{x - 0.3}mm", pos_y="0mm",
                                                              orientation=str(orient)))
        RouteMeander(design, f"Read_{name}", options=Dict(
            pin_inputs=Dict(start_pin=Dict(component=name, pin="readout"),
                            end_pin=Dict(component=f"{name}_Read_T", pin="second_end")),
            lead=Dict(start_straight="50um", end_straight="50um"), meander=Dict(spacing="130um", asymmetry="0um"),
            fillet="40um", total_length="3mm", **CPW))
    OpenToGround(design, "otg_left", options=dict(pos_x="-2.5mm", pos_y="0mm", orientation="180", width="10um", gap="6um"))
    OpenToGround(design, "otg_right", options=dict(pos_x="2.5mm", pos_y="0mm", orientation="0", width="10um", gap="6um"))
    RouteStraight(design, "feedline_main", options=dict(pin_inputs=dict(
        start_pin=dict(component="otg_left", pin="open"), end_pin=dict(component="otg_right", pin="open")), **CPW))
    return design


class TestSingleQubitCells(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp())
        cls.design = make_chip()
        cls.cells = SingleQubitCells(cls.design, cls.tmp / "sims", user_options=dict(fillet_resolution=12),
                                     mesh=dict(min_size=20e-6, max_size=400e-6), name="test")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_cell_components(self):
        self.assertEqual(self.cells.qubit_names(), ["Q1", "Q2"])
        self.assertEqual(self.cells.components("Q1"), ["Q1", "Q1_Read_T", "Read_Q1"])
        self.assertEqual(self.cells.components("Q2"), ["Q2", "Q2_Read_T", "Read_Q2"])

    def test_shared_component_must_be_excluded(self):
        cells = SingleQubitCells(self.design, self.tmp / "unused", user_options={}, exclude=r"^$")
        RouteStraight(self.design, "link", options=dict(pin_inputs=dict(
            start_pin=dict(component="Q1_Read_T", pin="prime_end"),
            end_pin=dict(component="Q2_Read_T", pin="prime_end")), **CPW))
        try:
            with self.assertRaises(ValueError):
                cells.components("Q1")
        finally:
            self.design.delete_component("link")

    def test_copied_geometry_matches_chip(self):
        for q in ("Q1", "Q2"):
            cell = self.cells.cell(q)
            self.assertEqual(len(cell.stubs), 2)
            self.assertLess(max(cell.geometry_difference_mm2.values()), 1e-12)
            for n in cell.components:
                self.assertEqual(type(cell.design.components[n]), type(self.design.components[n]))

    def test_conductors_identified_by_geometry(self):
        for q in ("Q1", "Q2"):
            rec = self.cells.simulate(q, run=False)
            self.assertEqual(rec["status"], "prepared")
            idx = rec["conductor_indices"]
            self.assertEqual(sorted(idx), ["feed", "ground", "pad1", "pad2", "res"])
            self.assertEqual(sorted(idx.values()), [0, 1, 2, 3, 4])
            self.assertTrue(Path(rec["plot"]).exists())
            # pad1 is the readout-side pad: nearer the readout pin than pad2
            sim = self.cells._prepare(f"test_{q}_check", self.cells.cell(q).design)
            t = sim._cur_cap_terminals
            pin = shapely.Point(np.array(self.cells.cell(q).design.components[q].pins["readout"]["middle"]))
            self.assertLess(t[idx["pad1"]].distance(pin), t[idx["pad2"]].distance(pin))
            self.assertEqual(idx, self.cells._conductors(q, sim))

    def test_params_needs_results(self):
        self.cells.simulate("Q1", run=False)
        with self.assertRaises(ValueError):
            self.cells.params("Q1", f_q_GHz=5.0, f_r_GHz=7.0)

    def test_to_catalogue(self):
        from SQDMetal.Utilities.QubitCatalogue import QubitCatalogue
        from UnitTests.testQubitCatalogue import expected_ec, write_terminal_csv
        cells = SingleQubitCells(self.design, self.tmp / "sims_cat", user_options=dict(fillet_resolution=12),
                                 mesh=dict(min_size=20e-6, max_size=400e-6), name="test")
        rec = cells.simulate("Q2", label="qubit2", run=False)
        caps = {("pad1", "ground"): 78.4, ("pad2", "ground"): 130.6, ("pad1", "pad2"): 79.7, ("pad1", "res"): 25.6,
                ("pad2", "res"): 1.6, ("res", "feed"): 3.4, ("res", "ground"): 900.0, ("feed", "ground"): 300.0}
        write_terminal_csv(Path(rec["output_dir"]).parent, rec["conductor_indices"], caps)   # stands in for Palace
        self.assertEqual(cells.simulate("Q2", label="qubit2")["status"], "cached")
        cat = QubitCatalogue(self.tmp / "catalogue")
        key = cells.to_catalogue(cat, "Q2", f_q_GHz=5.0, f_r_GHz=7.0)
        self.assertEqual(cells.to_catalogue(cat, "Q2", f_q_GHz=5.0, f_r_GHz=7.0), key)
        row = cat.collect().set_index("key").loc[key]
        self.assertAlmostEqual(row["E_C_GHz"], expected_ec(caps), places=9)
        self.assertEqual(row["source.qubit"], "Q2")
        self.assertEqual(row["source.label"], "qubit2")
        self.assertEqual(row["cell.builder"], "chip")
        self.assertEqual(row["cell.components"], "CoupledLineTee,RouteMeander")
        self.assertEqual(row["pad_gap"], "20um")
        self.assertTrue(np.isfinite(row["result.g_MHz"]) and row["result.g_MHz"] != 0)
        self.assertTrue((cat._folder(key) / "conductors.png").exists())

    def test_junction_for_round_trip(self):
        for f01, EC in ((2.0, 0.14), (5.0, 0.22)):
            j = SingleQubitCells.junction_for(f01, EC)
            self.assertAlmostEqual(j["f01_GHz"], f01, places=9)
            self.assertAlmostEqual(SingleQubitCells.transmon_levels(j["E_J_GHz"], EC)[1], f01, places=9)
            # L_J = (Phi0 / 2pi)^2 / E_J
            h, e = 6.62607015e-34, 1.602176634e-19
            self.assertAlmostEqual(j["L_J_nH"], (h / (2 * e) / (2 * np.pi)) ** 2 / (h * j["E_J_GHz"] * 1e9) * 1e9, places=9)
            self.assertLess(j["alpha_MHz"], 0)
        with self.assertRaises(ValueError):
            SingleQubitCells.junction_for(1e-4, 0.2)       # below what E_J/E_C >= 0.5 can give

    def test_fingerprint_tracks_cell_geometry(self):
        fp = self.cells.fingerprint("Q1")
        self.assertEqual(fp, self.cells.fingerprint("Q1"))
        self.assertNotEqual(fp, self.cells.fingerprint("Q2"))
        old = self.design.components["Read_Q1"].options.total_length
        self.design.components["Read_Q1"].options.total_length = "3.2mm"
        try:
            self.assertNotEqual(fp, self.cells.fingerprint("Q1"))
        finally:
            self.design.components["Read_Q1"].options.total_length = old

    def test_transmon_levels(self):
        E = SingleQubitCells.transmon_levels(EJ=20.0, EC=0.25)
        f01, alpha = E[1], E[2] - 2 * E[1]
        self.assertAlmostEqual(f01, np.sqrt(8 * 20.0 * 0.25) - 0.25, delta=0.02)
        self.assertAlmostEqual(alpha, -0.25, delta=0.03)

    def test_resonator_length_and_capacitance(self):
        from SQDMetal.Utilities.QUtilities import QUtilities
        L = self.cells.resonator_length_mm("Q1")
        route = QUtilities.calc_length_of_path(self.design, "Read_Q1")
        tee = QUtilities.calc_length_of_path(self.design, "Q1_Read_T", trace_name="second_cpw")
        pad = self.design.components["Q1"].get_resonator_length_mm()
        self.assertAlmostEqual(L, route + tee + pad, places=9)          # design units are mm
        self.assertAlmostEqual(route, 3.0, places=6)
        # capacitance: Maxwell self-capacitance of the resonator conductor
        from UnitTests.testQubitCatalogue import write_terminal_csv
        cells = SingleQubitCells(self.design, self.tmp / "sims_res", user_options=dict(fillet_resolution=12),
                                 mesh=dict(min_size=20e-6, max_size=400e-6), name="test")
        rec = cells.simulate("Q1", run=False)
        caps = {("pad1", "res"): 25.0, ("pad2", "res"): 2.0, ("res", "feed"): 3.0, ("res", "ground"): 600.0,
                ("pad1", "ground"): 80.0, ("pad2", "ground"): 120.0, ("pad1", "pad2"): 70.0, ("feed", "ground"): 200.0}
        write_terminal_csv(Path(rec["output_dir"]).parent, rec["conductor_indices"], caps)
        cells.simulate("Q1")
        self.assertAlmostEqual(cells.resonator_capacitance_fF("Q1"), 630.0, places=6)

    def test_resonator_calibration(self):
        # resonators with the same CPW: f = K / sqrt(l C) + pull
        K, lengths, caps = 870.0, (12.0, 9.5), (2500.0, 1990.0)
        pulls = (0.002, 0.004)
        f_true = [K / np.sqrt(l * c) + p for l, c, p in zip(lengths, caps, pulls)]
        cal = SingleQubitCells.resonator_calibration(f_true[0], lengths[0], caps[0], pulls[0], qubit="Q1")
        self.assertAlmostEqual(cal["K"], K, places=9)
        self.assertEqual(cal["qubit"], "Q1")
        self.assertAlmostEqual(SingleQubitCells.predict_resonator(cal, lengths[1], caps[1], pulls[1]), f_true[1], places=12)
        # design correction: built for 5.0, predicted 4.8 -> build for 5.0 * 4.9 / 4.8 to land on 4.9
        self.assertAlmostEqual(SingleQubitCells.corrected_design_frequency(4.9, 5.0, 4.8), 5.0 * 4.9 / 4.8)
        self.assertAlmostEqual(SingleQubitCells.qubit_pull_GHz(5.0, 3.0, -60.0), 0.06 ** 2 / 2.0)

    def test_resonator_capacitance_fit(self):
        # C = 215 l - 70 with +-1 % alternating mesh scatter: the fit recovers the line
        lengths = np.array([9.6, 11.5, 12.0, 12.3, 12.6, 12.9])
        scatter = np.array([1, -1, 1, -1, 1, -1]) * 0.01
        caps = (215.0 * lengths - 70.0) * (1 + scatter)
        fit = SingleQubitCells.fit_resonator_capacitance(lengths, caps)
        self.assertEqual(fit["n"], 6)
        self.assertAlmostEqual(fit["C_per_mm_fF"], 215.0, delta=10)
        self.assertAlmostEqual(fit["rms_pct"], 1.0, delta=0.2)
        with self.assertRaises(ValueError):
            SingleQubitCells.fit_resonator_capacitance([12.0, 12.0], [2500.0, 2510.0])
        # exact line: the length solved for a target frequency gives back that frequency
        fit = SingleQubitCells.fit_resonator_capacitance(lengths, 215.0 * lengths - 70.0)
        C_cal = SingleQubitCells.fitted_capacitance_fF(fit, 12.685)
        cal = SingleQubitCells.resonator_calibration(4.18386, 12.685, C_cal, 0.0032)
        l = SingleQubitCells.resonator_length_for(cal, fit, 4.3, 0.0025)
        f = SingleQubitCells.predict_resonator(cal, l, SingleQubitCells.fitted_capacitance_fF(fit, l), 0.0025)
        self.assertAlmostEqual(f, 4.3, places=10)
        self.assertLess(l, 12.685)

    def test_eig_rows(self):
        f = self.tmp / "eig.csv"
        self.assertEqual(SingleQubitCells._eig_rows(f), 0)
        f.write_text("        m,  Re{f} (GHz),  Im{f} (GHz),  Q\n")         # header only: interrupted run
        self.assertEqual(SingleQubitCells._eig_rows(f), 0)
        f.write_text(f.read_text() + "1, 2.8, 1e-6, 1e6\n2, 4.6, 1e-4, 2e4\n")
        self.assertEqual(SingleQubitCells._eig_rows(f), 2)

    def test_design_from_notebook(self):
        nb = nbformat.v4.new_notebook()
        nb.cells = [nbformat.v4.new_code_cell("%matplotlib inline\nx = 1"),
                    nbformat.v4.new_code_cell("y = x + 1  # BUILD"),
                    nbformat.v4.new_code_cell("z = 3  # after the build cell")]
        path = self.tmp / "design.ipynb"
        nbformat.write(nb, path)
        ns = SingleQubitCells.design_from_notebook(path, "# BUILD")
        self.assertEqual(ns["y"], 2)
        self.assertNotIn("z", ns)       # stops at the build cell
        with self.assertRaises(ValueError):
            SingleQubitCells.design_from_notebook(path, "missing marker")


if __name__ == "__main__":
    unittest.main()
