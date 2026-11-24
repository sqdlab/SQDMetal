import contextlib
import io
import json
import re
import unittest

import numpy as np

from SQDMetal.Utilities.QubitDesigner import (
    FloatingTransmonDesigner,
    ResonatorQuarterWave,
    XmonDesigner,
)


class TestQubitOptimisation(unittest.TestCase):
    def test_fixed_xmon_reports_original_constraint_residuals(self):
        resonator = ResonatorQuarterWave(7e9)
        designer = XmonDesigner(resonator)
        params = {"fQubit": 4e9, "C_g": 5e-15, "C_J": 60e-15, "chi": -2e6, "Ej/Ec": 50}
        cr = resonator.get_res_capacitance()
        chi = designer._chi_hertz(
            params["C_g"],
            params["C_J"],
            params["fQubit"],
            resonator.get_res_frequency(),
            cr,
            resonator.get_res_inductance(),
        )
        ratio = designer.EJonEC(
            designer._eff_CJ(params["C_J"], params["C_g"], cr), params["fQubit"]
        )
        expected = 50 * (
            abs(chi / params["chi"] - 1) + abs(ratio / params["Ej/Ec"] - 1)
        )

        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            designer.optimise(params)
        # Check the numerical objective, not the report's wording.
        reported = float(
            re.search(
                r"[-+]?\d+(?:\.\d*)?(?:e[-+]?\d+)?",
                output.getvalue().splitlines()[0],
                re.IGNORECASE,
            ).group()
        )
        self.assertAlmostEqual(reported, expected, places=10)
        self.assertGreater(expected, 1)

    def test_floating_result_preserves_target_cost_and_serializes(self):
        designer = FloatingTransmonDesigner(ResonatorQuarterWave(7e9))
        params = {
            "fQubit": 4e9,
            "C_q1": 30e-15,
            "C_q2": 35e-15,
            "C_g1": 5e-15,
            "C_g2": 1e-15,
            "C_J": 40e-15,
            "chi": -2e6,
            "Ej/Ec": 50,
            "C_sigma": 100e-15,
            "beta": 0.1,
        }
        result = designer.optimise(params, print_results=False)
        residuals = [
            ("chi_Hz", "chi"),
            ("Ej_Ec_ratio", "Ej/Ec"),
            ("C_sigma_F", "C_sigma"),
            ("beta", "beta"),
        ]
        expected = 50 * sum(
            abs(result[field] / params[target] - 1) for field, target in residuals
        )
        self.assertTrue(np.isfinite(expected))
        self.assertAlmostEqual(result["cost_error_%"], expected, places=10)
        decoded = json.loads(json.dumps(result))
        self.assertAlmostEqual(decoded["cost_error_%"], expected, places=10)
        for field, target in [
            ("fQubit_Hz", "fQubit"),
            ("C_q1_F", "C_q1"),
            ("C_q2_F", "C_q2"),
            ("C_g1_F", "C_g1"),
            ("C_g2_F", "C_g2"),
            ("C_J_F", "C_J"),
        ]:
            self.assertEqual(result[field], params[target])


if __name__ == "__main__":
    unittest.main()
