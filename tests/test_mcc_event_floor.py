"""The diagnostic floor changes a numerical stop, never MCC equations/gates."""
import sys
from pathlib import Path
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'tools'), str(ROOT/'interim/2026-11-09-nc-pseudo-arclength')]
from independent_branches import StrainRate
from mcc_rate_rk4 import advance_rk4 as original
from mcc_rate_event_floor import advance_rk4 as diagnostic


class EventFloor(unittest.TestCase):
    def test_nc_isotropic_compression_keeps_existing_lambda_closed_form(self):
        state = np.zeros((1, 14)); state[0, :3] = 100.; state[0, 6:8] = np.log([100., 2.])
        strain = np.array([[.001, .001, .001, 0., 0., 0.]])
        parameters = np.array([[1.44, .2, .04, .3]])
        expected_v = 2*np.exp(-.003)
        expected_p = 100*np.exp((2-expected_v)/.2)
        result, _, audit = diagnostic(StrainRate(parameters), state, np.zeros(1), strain, 32)
        np.testing.assert_allclose(result[0, :3], expected_p, rtol=0, atol=1e-9)
        self.assertAlmostEqual(np.exp(result[0, 6]), expected_p, places=9)
        self.assertAlmostEqual(np.exp(result[0, 7]), expected_v, places=12)
        self.assertLess(audit['maximum_normalized_yield_drift'], 2e-6)
        baseline, _, _ = original(StrainRate(parameters), state, np.zeros(1), strain, 32)
        np.testing.assert_array_equal(result, baseline)

    def test_elastic_unloading_and_zero_increment_keep_hardening_and_history(self):
        state = np.zeros((1, 14)); state[0, :3] = 100.; state[0, 6:8] = np.log([200., 2.])
        parameters = np.array([[1.44, .2, .04, .3]])
        for strain in (np.zeros((1, 6)), np.array([[-1e-4, -1e-4, -1e-4, 2e-5, 0., 0.]])):
            result, _, _ = diagnostic(StrainRate(parameters), state, np.zeros(1), strain, 8)
            baseline, _, _ = original(StrainRate(parameters), state, np.zeros(1), strain, 8)
            np.testing.assert_array_equal(result, baseline)
            np.testing.assert_array_equal(result[:, 6], state[:, 6])
            np.testing.assert_array_equal(result[:, 8:], state[:, 8:])

    def test_invalid_floor_is_not_an_acceptance_override(self):
        for floor in (0., -1., 1e-10, float('nan')):
            with self.assertRaisesRegex(ValueError, 'event floor'):
                diagnostic(None, None, None, None, 4, minimum_event_step=floor)


if __name__ == '__main__':
    unittest.main()
