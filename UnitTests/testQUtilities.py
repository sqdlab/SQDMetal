import unittest

from SQDMetal.Utilities.QUtilities import QUtilities


class TestParseValueLength(unittest.TestCase):
    def test_prefixed_units(self):
        cases = {
            '5mm': 5e-3,
            '12 um': 12e-6,
            '0.5mm': 0.5e-3,
            '5um': 5e-6,
            '5nm': 5e-9,
            '5pm': 5e-12,
            '5fm': 5e-15,
        }
        for text, expected in cases.items():
            with self.subTest(text=text):
                self.assertAlmostEqual(QUtilities.parse_value_length(text), expected,
                                       delta=abs(expected)*1e-12)

    def test_metre_suffix(self):
        # Regression: the single-'m' branch used to strip two characters, so
        # '5m' and '1m' raised ValueError and '10m'/'0.5m' lost their value.
        cases = {
            '1m': 1.0,
            '5m': 5.0,
            '10m': 10.0,
            '0.5m': 0.5,
            '2.5 m': 2.5,
        }
        for text, expected in cases.items():
            with self.subTest(text=text):
                self.assertAlmostEqual(QUtilities.parse_value_length(text), expected,
                                       delta=abs(expected)*1e-12)

    def test_numeric_passthrough(self):
        self.assertEqual(QUtilities.parse_value_length(5), 5)
        self.assertEqual(QUtilities.parse_value_length(2.5), 2.5)

    def test_invalid_input(self):
        # An unsupported unit is returned unchanged, as documented.
        self.assertEqual(QUtilities.parse_value_length('5cm'), '5cm')
        with self.assertRaises(AssertionError):
            QUtilities.parse_value_length('m')


if __name__ == '__main__':
    unittest.main()
