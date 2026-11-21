"""Check step recovery against the existing integrator and an elastic solution."""
import importlib.util
from pathlib import Path
import sys
import time
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'tools'), str(ROOT/'interim/2026-11-09-nc-pseudo-arclength')]
from independent_branches import StrainRate
from mcc_rate_rk4 import advance_rk4 as fixed
SPEC = importlib.util.spec_from_file_location('nc_regrowth', ROOT/'interim/2026-11-10-nc-mean-stress-release/rk4_regrowth.py')
recovery = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(recovery)


class Rk4Regrowth(unittest.TestCase):
    def state_and_system(self):
        # Existing continuum MCC state layout; a material-only analytical check.
        state = np.zeros((1, 14)); state[0, :3] = 100.
        state[0, 6:8] = np.log([100., 2.])
        return state, StrainRate(np.array([[1.2**2, .2, .05, .3]]))

    def test_elastic_isotropic_expansion_converges_to_analytical_pressure(self):
        state, system = self.state_and_system(); strain = np.array([[-.01/3]*3+[0.]*3])
        exact = 100.*np.exp(2./.05*(1-np.exp(.01)))
        errors = []
        for steps in (4, 8):
            actual, _, _ = recovery.advance_rk4(system, state, np.zeros(1), strain, steps)
            errors.append(abs(actual[0, 0]-exact))
            self.assertAlmostEqual(float(np.exp(actual[0, 7])), 2.*np.exp(.01), places=13)
            np.testing.assert_array_equal(actual[:, 8:], 0.)
            self.assertEqual(actual[0, 6], state[0, 6])
        self.assertLess(errors[1], errors[0]/10)

    def test_unshortened_path_is_exactly_the_existing_rk4(self):
        state, system = self.state_and_system(); strain = np.array([[-.001/3]*3+[0.]*3])
        a, ua, _ = fixed(system, state, np.zeros(1), strain, 4)
        b, ub, _ = recovery.advance_rk4(system, state, np.zeros(1), strain, 4)
        np.testing.assert_array_equal(a, b); np.testing.assert_array_equal(ua, ub)

    def test_deadline_inside_material_integration_preserves_input(self):
        state, system = self.state_and_system(); before = state.copy(); system.deadline = time.monotonic()-1
        with self.assertRaisesRegex(ValueError, 'wall time budget inside'):
            recovery.advance_rk4(system, state, np.zeros(1), np.zeros((1, 6)), 4)
        np.testing.assert_array_equal(state, before)

    def test_over_drift_state_is_rejected_without_projection(self):
        state, system = self.state_and_system(); state[0, 6] = np.log(99.); before = state.copy()
        with self.assertRaisesRegex(ValueError, 'drift cannot be controlled without a projection'):
            recovery.advance_rk4(system, state, np.zeros(1), np.zeros((1, 6)), 4)
        np.testing.assert_array_equal(state, before)


if __name__ == '__main__':
    unittest.main()
