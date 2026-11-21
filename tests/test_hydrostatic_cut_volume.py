"""Cut-volume moments, force work and causal ablation on the existing column."""
import copy
import json
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import fem_reference as fem
from hydrostatic_cut_volume_reference import force
from diagnose_wet_joint_cut_volume import components, linear_initial_correction, recorded_inputs
from wet_joint_cycle_reference import cut_volume_reference, reference


class HydrostaticCutVolume(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mesh, cls.ground, cls.case = [json.loads((ROOT/'examples'/(name+'.json')).read_text()) for name in (
            'joint-water-column-tet10-mesh', 'ground-joint-water-mc', 'stages-joint-water-cycle')]
        cls.ground = fem.normalize_ground(cls.ground); cls.case = fem.normalize_stage_case(cls.case)

    def test_cube_analytic_moments_and_force_virtual_work(self):
        pressure = 0.; volume = 0.; moment = np.zeros(3)
        xyz = np.asarray(self.mesh['nodes'])
        for ei in range(12):
            value = force(self.mesh, ei, .75, 10)
            f = value['force']; p = value['pressure_integral']
            pressure += p; volume += value['wet_volume']; moment += value['pressure_first_moment']
            np.testing.assert_allclose(f.sum(axis=0), 0, atol=2e-14)
            np.testing.assert_allclose(np.cross(xyz[value['nodes']], f).sum(axis=0), 0, atol=2e-14)
            # Exact work against every affine displacement gradient.
            np.testing.assert_allclose(xyz[value['nodes']].T@f, -p*np.eye(3), atol=2e-14)
            for order in (3, 5):
                np.testing.assert_allclose(f, force(self.mesh, ei, .75, 10, order)['force'], rtol=1e-12, atol=2e-14)
        self.assertAlmostEqual(volume, .75, places=13)
        self.assertAlmostEqual(pressure, 10*.75**2/2, places=13)
        np.testing.assert_allclose(moment, [10*.75**2/4, 10*.75**2/4, 10*.75**3/6], atol=2e-14)

    def test_dry_vertex_face_partial_and_full_wetting_tet4_and_tet10(self):
        for quadratic in (False, True):
            mesh = copy.deepcopy(self.mesh)
            if not quadratic:
                del mesh['tetrahedra10']
            for level in (0., .5, .75, 1., 1.25):
                rows = [force(mesh, ei, level, 10) for ei in range(12)]
                expected_volume = min(level, 1.)
                expected_pressure = 10*(level**2/2 if level <= 1 else level-.5)
                self.assertAlmostEqual(sum(r['wet_volume'] for r in rows), expected_volume, places=13)
                self.assertAlmostEqual(sum(r['pressure_integral'] for r in rows), expected_pressure, places=13)
                for ei, r in enumerate(rows):
                    np.testing.assert_allclose(r['force'], force(mesh, ei, level, 10, 3)['force'], atol=2e-14)

    def test_curved_or_nonfinite_geometry_is_outside_reference(self):
        bad = copy.deepcopy(self.mesh); bad['nodes'][bad['tetrahedra10'][0][4]][2] += .001
        with self.assertRaisesRegex(AssertionError, 'curved Tet10'):
            force(bad, 0, .75, 10)
        with self.assertRaisesRegex(AssertionError, 'finite water'):
            force(self.mesh, 0, float('nan'), 10)
        with self.assertRaisesRegex(AssertionError, 'unit weight'):
            force(self.mesh, 0, .75, -1)

    def test_only_partial_pore_term_accounts_for_initial_displacement_change(self):
        original = copy.deepcopy((self.mesh, self.ground, self.case))
        sampled = components(self.mesh, self.ground, 'gl', 7)
        cut = components(self.mesh, self.ground, 'cut', 7)
        np.testing.assert_array_equal(sampled['stiffness'], cut['stiffness'])
        np.testing.assert_array_equal(sampled['effective_seed'], cut['effective_seed'])
        self.assertGreater(np.linalg.norm(sampled['pore_force']-cut['pore_force']), 1e-5)
        for row in sampled['elements']:
            if not row['partially_wet']:
                self.assertLess(row['force_error']['maximum_absolute'], 2e-14)
        a = reference(self.mesh, self.ground, self.case, 4, 8, 7)
        b = cut_volume_reference(self.mesh, self.ground, self.case, 4, 8, 7)
        correction = linear_initial_correction(self.mesh, self.ground, self.case, sampled, cut)
        actual = b['frames'][0]['u']-a['frames'][0]['u']
        np.testing.assert_allclose(actual, correction, rtol=1e-6, atol=1e-12)
        self.assertEqual((self.mesh, self.ground, self.case), original)
        self.assertEqual(b['volume_rule'], 'affine_cut_hydrostatic')
        self.assertFalse(b['reference_uses_native_results'])
        self.assertEqual(len(b['path']), 28)
        self.assertLess(b['maximum_relative_force_residual'], 1e-10)
        for x, y in zip(a['frames'], b['frames']):
            np.testing.assert_array_equal(x['physical_pressure'], y['physical_pressure'])
            np.testing.assert_array_equal(x['water_load'], y['water_load'])
            np.testing.assert_array_equal(x['pore'], y['pore'])

    def test_recorded_256_histories_are_read_with_provenance_not_recomputed(self):
        values, provenance = recorded_inputs()
        self.assertFalse(provenance['rerun'])
        self.assertEqual(set(values), {'discrete', 'gl7'})
        self.assertEqual(len(provenance['entries']), 5)
        for r in values.values():
            self.assertEqual(r['steps'], 256)
            self.assertEqual(len(r['frames']), 7)


if __name__ == '__main__':
    unittest.main()
