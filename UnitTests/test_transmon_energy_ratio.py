import unittest

import numpy as np
from scipy.constants import e, h

from SQDMetal.Utilities.QubitDesigner import TransmonBase


class TestTransmonEnergyRatio(unittest.TestCase):
    def test_recovers_josephson_to_charging_energy_ratio(self):
        transmon = TransmonBase()
        for capacitance in (30e-15, 60e-15, 120e-15):
            charging_energy = e**2 / (2 * capacitance)
            for ratio in (20, 50, 100):
                with self.subTest(capacitance=capacitance, ratio=ratio):
                    frequency = (np.sqrt(8 * ratio) - 1) * charging_energy / h
                    # Allow for the rounded electron charge used by the designer.
                    np.testing.assert_allclose(
                        transmon.EJonEC(capacitance, frequency), ratio, rtol=1e-7
                    )


if __name__ == "__main__":
    unittest.main()
