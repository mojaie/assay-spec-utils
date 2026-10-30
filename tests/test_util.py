import unittest

from assay_spec_utils.util import *


class TestConvertConcUnits(unittest.TestCase):
    def test_molar(self):
        for unit, expected in [("M", 2), ("mM", 2e-3), ("uM", 2e-6),
                               ("nM", 2e-9), ("pM", 2e-12)]:
            with self.subTest(unit=unit):
                value, base = convert_conc_units(2, unit)
                self.assertAlmostEqual(value, expected)
                self.assertEqual(base, "M")

    def test_invalid_molar(self):
        with self.assertRaises(ValueError):
            convert_conc_units(-1, "M")

    def test_mass_per_volume(self):
        self.assertEqual(convert_conc_units(5, "mg/mL"), (5, "g/L"))
        self.assertAlmostEqual(convert_conc_units(5, "ng/uL")[0], 5e-3)

    def test_ratio(self):
        self.assertAlmostEqual(convert_conc_units(0.5, "v/v%")[0], 0.005)
        self.assertEqual(convert_conc_units(0.5, "v/v%")[1], "ratio")
        self.assertEqual(convert_conc_units(200, "dilution"), (0.005, "ratio"))
        self.assertEqual(convert_conc_units(0.2, "ratio"), (0.2, "ratio"))

    def test_invalid_values(self):
        for value, unit in [(101, "v/v%"), (0.5, "dilution"), (1.5, "ratio"),
                            (1, "mol")]:
            with self.subTest(value=value, unit=unit):
                with self.assertRaises(ValueError):
                    convert_conc_units(value, unit)


class TestConvertTimeUnits(unittest.TestCase):
    def test_units(self):
        for unit, expected in [("sec", 2), ("min", 120), ("hour", 7200),
                               ("day", 172800)]:
            with self.subTest(unit=unit):
                self.assertEqual(convert_time_units(2, unit), (expected, "sec"))

    def test_invalid_values(self):
        with self.assertRaises(ValueError):
            convert_time_units(-1, "sec")
        with self.assertRaises(ValueError):
            convert_time_units(1, "week")


class TestIsConvertible(unittest.TestCase):
    def test_int(self):
        self.assertTrue(is_convertible_to_int("10"))
        self.assertTrue(is_convertible_to_int(1.5))
        self.assertFalse(is_convertible_to_int("1.5"))
        self.assertFalse(is_convertible_to_int("abc"))

    def test_float(self):
        self.assertTrue(is_convertible_to_float("1.5"))
        self.assertTrue(is_convertible_to_float("1e-3"))
        self.assertFalse(is_convertible_to_float("abc"))


if __name__ == "__main__":
    unittest.main()
