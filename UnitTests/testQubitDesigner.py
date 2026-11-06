import unittest

import numpy as np

from SQDMetal.Utilities.QubitDesigner import ResonatorLC


class TestResonatorLC(unittest.TestCase):
    def assertRelAlmostEqual(self, actual, expected, rel=1e-12):
        self.assertAlmostEqual(actual/expected, 1.0, delta=rel)

    def test_from_inductance_and_capacitance(self):
        L, C = 1e-9, 1e-12
        resonator = ResonatorLC(L=L, C=C)
        self.assertRelAlmostEqual(resonator.get_res_inductance(), L)
        self.assertRelAlmostEqual(resonator.get_res_capacitance(), C)
        self.assertRelAlmostEqual(resonator.get_res_frequency(), 1/(2*np.pi*np.sqrt(L*C)))
        self.assertRelAlmostEqual(resonator.get_res_impedance(), np.sqrt(L/C))
        self.assertTrue(resonator.is_res_parallelLC())

    def test_from_frequency_and_capacitance(self):
        f0, C = 7.25e9, 1e-12
        resonator = ResonatorLC(f0=f0, C=C)
        self.assertEqual(resonator.get_res_frequency(), f0)
        self.assertRelAlmostEqual(resonator.get_res_capacitance(), C)
        self.assertRelAlmostEqual(resonator.get_res_inductance(), 1/((2*np.pi*f0)**2*C))
        self.assertRelAlmostEqual(1/(2*np.pi*np.sqrt(resonator.get_res_inductance()
                                                      * resonator.get_res_capacitance())), f0)

    def test_from_frequency_and_inductance(self):
        f0, L = 7.25e9, 1e-9
        resonator = ResonatorLC(f0=f0, L=L)
        self.assertEqual(resonator.get_res_frequency(), f0)
        self.assertRelAlmostEqual(resonator.get_res_inductance(), L)
        self.assertRelAlmostEqual(resonator.get_res_capacitance(), 1/((2*np.pi*f0)**2*L))
        self.assertRelAlmostEqual(1/(2*np.pi*np.sqrt(resonator.get_res_inductance()
                                                      * resonator.get_res_capacitance())), f0)

    def test_all_three_parameters(self):
        L, C = 1e-9, 1e-12
        f0 = 1/(2*np.pi*np.sqrt(L*C))
        resonator = ResonatorLC(f0=f0, L=L, C=C)
        self.assertEqual(resonator.get_res_frequency(), f0)
        self.assertRelAlmostEqual(resonator.get_res_inductance(), L)
        self.assertRelAlmostEqual(resonator.get_res_capacitance(), C)
        rounded = ResonatorLC(f0=f0*(1+1e-14), L=L, C=C)
        self.assertRelAlmostEqual(rounded.get_res_inductance(), L)
        self.assertRelAlmostEqual(rounded.get_res_capacitance(), C)
        with self.assertRaises(AssertionError):
            ResonatorLC(f0=f0*(1+1e-8), L=L, C=C)
        with self.assertRaises(AssertionError):
            ResonatorLC(f0=10e9, L=L, C=C)

    def test_missing_parameters(self):
        with self.assertRaises(AssertionError):
            ResonatorLC(L=1e-9)
        with self.assertRaises(AssertionError):
            ResonatorLC(f0=7.25e9)


if __name__ == '__main__':
    unittest.main()
