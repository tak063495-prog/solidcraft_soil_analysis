"""Analytical checks of read-only MCC stress-path diagnostics."""
import copy
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
import diagnose_nc_material_path as diag

PARAM = dict(critical_state_m=1.4, **{'lambda': .1, 'kappa': .02, 'poisson': .25})


def state(stress, pc=100., plastic=None):
    return dict(effective_stress_compression_positive_kpa=stress,
                plastic_strain_compression_positive=plastic or [0.]*6,
                preconsolidation_pressure_kpa=pc, void_ratio=1.2)


class MaterialPath(unittest.TestCase):
    def test_hydrostatic_and_pure_shear_invariants(self):
        a = diag.invariants([state([100., 100., 100., 0., 0., 0.]),
                             state([50., 50., 50., 2., 0., 0.])], 'native', 1.4)
        np.testing.assert_array_equal(a['p'], [100., 50.])
        np.testing.assert_allclose(a['q'], [0., np.sqrt(12.)], rtol=1e-15)
        np.testing.assert_array_equal(a['flow_volume_normal'], [100., 0.])

    def test_tensor_rotation_preserves_p_q_and_flow_side(self):
        tensor = np.array([[60., 4., -2.], [4., 30., 3.], [-2., 3., 20.]])
        angle = .731
        c, s = np.cos(angle), np.sin(angle)
        rotation = np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])
        def packed(a):
            return [a[0, 0], a[1, 1], a[2, 2], a[0, 1], a[1, 2], a[2, 0]]
        a = diag.invariants([state(packed(tensor)), state(packed(rotation@tensor@rotation.T))], 'native', 1.4)
        for key in ('p', 'q', 'flow_volume_normal', 'yield_over_pc_squared'):
            self.assertAlmostEqual(a[key][0], a[key][1], places=12)

    def test_rate_normal_matches_independent_yield_finite_difference(self):
        stress = np.array([50., 30., 20., 4., 3., -2.])
        result = diag.continuum_terms([state(stress.tolist())], 'native', PARAM)
        def yield_value(x):
            tensor = np.array([[x[0], x[3], x[5]], [x[3], x[1], x[4]], [x[5], x[4], x[2]]])
            p = np.trace(tensor)/3
            dev = tensor-p*np.eye(3)
            return 1.5*np.trace(dev@dev)/1.4**2+p*(p-100.)
        h = 1e-4
        numerical = [(yield_value(stress+h*np.eye(6)[i])-yield_value(stress-h*np.eye(6)[i]))/(2*h) for i in range(6)]
        np.testing.assert_allclose(result['normal'][0], numerical, rtol=1e-9, atol=1e-8)

    def test_dry_softening_and_wet_hardening_have_opposite_plastic_eigenvalue_signs(self):
        samples = []
        for p in (25., 75.):
            q = 1.4*np.sqrt(p*(100.-p))
            samples.append(state([p+2*q/3, p-q/3, p-q/3, 0., 0., 0.]))
        r = diag.continuum_terms(samples, 'native', PARAM)
        self.assertLess(r['hardening'][0], 0.)
        self.assertGreater(r['hardening'][1], 0.)
        self.assertTrue(np.all(r['denominator'] > 0.))
        self.assertTrue(np.all(np.linalg.eigvalsh(r['elastic']) > 0.))
        self.assertLess(np.linalg.eigvalsh(r['plastic_loading'])[0, 0], 0.)
        self.assertGreater(np.linalg.eigvalsh(r['plastic_loading'])[1, 0], 0.)

    def test_native_and_independent_representations_report_same_physics(self):
        native = [state([60., 40., 30., 2., 3., 4.], plastic=[.01, -.02, .003, .004, 0., 0.])]
        x = native[0]
        independent = [x['effective_stress_compression_positive_kpa']+[np.log(100.), np.log1p(1.2)]+x['plastic_strain_compression_positive']]
        a, b = (diag.invariants(samples, kind, 1.4) for samples, kind in ((native, 'native'), (independent, 'independent')))
        for key in a:
            np.testing.assert_allclose(a[key], b[key], rtol=1e-14, atol=1e-14)

    def test_measured_dilation_and_softening_do_not_grant_acceptance(self):
        before = [state([30., 30., 30., 0., 0., 0.])]
        after = [state([29., 29., 29., 0., 0., 0.], pc=99., plastic=[-.001]*3+[0.]*3)]
        original = copy.deepcopy(after)
        r = diag.describe(after, 'native', 1.4, [0], [1.], before)
        self.assertEqual(r['increment']['dilating_and_softening_points'], 1)
        self.assertEqual(r['dry_side_volume_fraction'], 1.)
        self.assertFalse(r['numerical_acceptance'])
        self.assertFalse(r['physical_failure_confirmed'])
        self.assertEqual(after, original)

    def test_nonfinite_invalid_inventory_and_unknown_map_are_refused(self):
        good = [state([30., 30., 30., 0., 0., 0.])]
        for indices, weights in (([0, 0], [1., 1.]), ([1], [1.]), ([0], [-1.]), ([0], [float('nan')])):
            with self.assertRaises(ValueError):
                diag.describe(good, 'native', 1.4, indices, weights)
        bad = copy.deepcopy(good); bad[0]['void_ratio'] = float('nan')
        with self.assertRaises(ValueError):
            diag.invariants(bad, 'native', 1.4)
        with self.assertRaises(ValueError):
            diag.invariants(good, 'unknown', 1.4)
        with self.assertRaises(ValueError):
            diag.continuum_terms(good, 'native', {**PARAM, 'lambda': .02})


if __name__ == '__main__':
    unittest.main()
