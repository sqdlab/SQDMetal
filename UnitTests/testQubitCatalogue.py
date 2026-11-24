import json
import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from qiskit_metal import Dict

from SQDMetal.Utilities.QubitCatalogue import QubitCatalogue

E, H = 1.602176634e-19, 6.62607015e-34


def write_terminal_csv(folder, idx, caps_fF):
    """Writes a Palace-style terminal-C.csv (Maxwell matrix, F) for mutual capacitances caps_fF
    {(role_a, role_b): C}, with conductors numbered by idx {role: matrix index}."""
    n = len(idx)
    M = np.zeros((n, n))
    for (a, b), c in caps_fF.items():
        i, j = idx[a], idx[b]
        M[i, j] = M[j, i] = -c * 1e-15
    for i in range(n):
        M[i, i] = -M[i].sum()
    out = folder / "outputFiles"
    out.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(M, columns=[f"C[i][{k + 1}] (F)" for k in range(n)])
    df.insert(0, "i", np.arange(1, n + 1))
    df.to_csv(out / "terminal-C.csv", index=False)
    (out / "out.log").write_text("palace log")
    (out / "paraview").mkdir(exist_ok=True)
    (out / "paraview" / "big.vtu").write_text("x" * 100)


def expected_ec(c):
    C1, C2 = c[("pad1", "ground")] + c[("pad1", "res")], c[("pad2", "ground")] + c[("pad2", "res")]
    Cs = (C1 * C2 / (C1 + C2) + c[("pad1", "pad2")]) * 1e-15
    return E ** 2 / (2 * Cs) / H / 1e9


class TestQubitCatalogueImport(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.cat = QubitCatalogue(self.tmp / "cat", base_options=Dict(pad_gap="20um"))
        self.idx = dict(ground=0, res=1, pad2=2, pad1=3, feed=4)     # a non-default order, as Palace gives
        self.caps = {("pad1", "ground"): 78.4, ("pad2", "ground"): 130.6, ("pad1", "pad2"): 79.7,
                     ("pad1", "res"): 25.6, ("pad2", "res"): 1.6, ("res", "feed"): 3.4, ("res", "ground"): 900.0,
                     ("feed", "ground"): 300.0, ("pad1", "feed"): 1e-5, ("pad2", "feed"): 1e-5}

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def add(self, pad_gap="20um", caps=None, **kw):
        folder = self.tmp / f"sim_{pad_gap}_{len(list(self.tmp.iterdir()))}"
        write_terminal_csv(folder, self.idx, caps or self.caps)
        opts = Dict(pad_gap=pad_gap, pos_x="1mm", pos_y="2mm", orientation="180",
                    connection_pads=Dict(readout=Dict(pad_gap="10um")))
        return self.cat.add_simulation(folder, opts, self.idx, cell=dict(builder="chip"),
                                       mesh=dict(min_size=1e-5), sim_options=dict(solver_order=2, palace_dir="/x"),
                                       source=dict(type="chip", design="D", qubit="Q1", cell_key="abc"),
                                       results=dict(g_MHz=70.0, kappa_MHz=0.05), **kw)

    def test_entry_files_and_record(self):
        key = self.add()
        folder = self.cat._folder(key)
        self.assertTrue((folder / "outputFiles" / "terminal-C.csv").exists())
        self.assertTrue((folder / "outputFiles" / "out.log").exists())
        self.assertFalse((folder / "outputFiles" / "paraview").exists())
        rec = json.loads((folder / "variant.json").read_text())
        self.assertNotIn("pos_x", rec["spec"]["qubit"])
        self.assertNotIn("palace_dir", rec["spec"]["sim"])
        self.assertEqual(rec["spec"]["cell"]["builder"], "chip")
        # the variant lists what differs from the class defaults + base options
        self.assertNotIn("pad_gap", rec["spec"]["variant"])
        self.assertEqual(rec["spec"]["variant"]["connection_pads"]["readout"]["pad_gap"], "10um")

    def test_collect_reads_imported_entry(self):
        self.add()
        row = self.cat.collect().iloc[0]
        self.assertEqual(row["status"], "done")
        self.assertAlmostEqual(row["E_C_GHz"], expected_ec(self.caps), places=9)
        self.assertAlmostEqual(row["C1_readout_fF"], 25.6, places=6)
        self.assertAlmostEqual(row["C1_feed_fF"], 1e-5, places=9)
        self.assertEqual(row["source.design"], "D")
        self.assertEqual(row["result.g_MHz"], 70.0)
        self.assertEqual(row["cell.builder"], "chip")

    def test_same_simulation_is_not_duplicated(self):
        self.assertEqual(self.add(), self.add())
        self.assertEqual(len(self.cat.keys()), 1)

    def test_run_refuses_imported_entry_without_results(self):
        key = self.add()
        (self.cat._folder(key) / "outputFiles" / "terminal-C.csv").unlink()
        self.cat.sim_options.palace_dir = shutil.which("true")
        with self.assertRaises(ValueError):
            self.cat.run([key])

    def test_geometry_for_keeps_chip_entries_apart(self):
        for gap, c12 in (("20um", 79.7), ("40um", 60.0)):
            self.add(gap, caps={**self.caps, ("pad1", "pad2"): c12})
        self.cat.collect()
        # by default only the catalogue's own unit cells are used -> none here
        with self.assertRaises(ValueError):
            self.cat.geometry_for(0.15, "pad_gap")
        # chip entries on their own
        table = self.cat.table.sort_values("pad_gap_um")
        target = float(table["E_C_GHz"].mean())
        out = self.cat.geometry_for(target, "pad_gap", use_current=False, fixed={"cell.builder": "chip"})
        self.assertTrue(20 < float(out.pad_gap[:-2]) < 40)

    def test_missing_results_raise(self):
        with self.assertRaises(FileNotFoundError):
            self.cat.add_simulation(self.tmp / "nothing", Dict(), self.idx, cell=dict(builder="chip"))


if __name__ == "__main__":
    unittest.main()
