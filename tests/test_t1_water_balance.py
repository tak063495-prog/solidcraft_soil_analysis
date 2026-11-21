import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from t1_water_balance import Account, LIMITS


class WaterBalanceTest(unittest.TestCase):
    def test_equal_throughflow_and_zero_inventory_change(self):
        for profile in LIMITS:
            row = Account(1.).record(1., 1., .01, .01, profile)
            self.assertEqual((row["e"], row["D"], row["E"], row["U"]), (0., .01, 0., 0.))
            self.assertTrue(row["passed"])
            self.assertIsNone(Account(0).record(0, 0, 0, 0, profile)["ratio"])

    def test_absolute_limit_inclusive_in_both_profiles(self):
        for profile, (_, _, floor) in LIMITS.items():
            self.assertTrue(Account(0).record(0, floor, 0, 0, profile)["passed"])
            above = math.nextafter(floor, math.inf)
            self.assertFalse(Account(0).record(0, above, 0, 0, profile)["passed"])

    def test_local_and_temporal_cancellation(self):
        first = Account(0.)
        first.record(0., 8e-13, 0., 0., "verification")
        row = first.record(8e-13, 0., 0., 0., "verification")
        self.assertEqual(row["E"], 0.)
        self.assertEqual(row["U"], 1.6e-12)
        self.assertFalse(row["passed"])
        self.assertTrue(Account(1.).record(1., 1., 0., 0., "verification")["passed"])
        for sign in (-1., 1.):
            self.assertFalse(Account(.5).record(.5, .5+sign*1e-8, 0., 0., "verification")["passed"])

    def test_relative_limits_and_event_units(self):
        for profile, (step, cumulative, _) in LIMITS.items():
            row = Account(0).record(0, .002, .002, 0, profile)
            self.assertEqual(row["D"], .002)
            self.assertEqual(row["step_limit"], step*.002)
            self.assertEqual(row["cumulative_limit"], cumulative*.002)
            self.assertTrue(row["passed"])
