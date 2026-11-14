"""Analytic counterexamples and tensor/sign checks for read-only NC diagnostics."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from nc_physical_limit_audit import stress_values, cone_lambda, resultant, face_label, run, sign


def state(stress, pc=150., plastic=None):
    return dict(effective_stress_compression_positive_kpa=stress, preconsolidation_pressure_kpa=pc,
                plastic_strain_compression_positive=plastic or [0.]*6, void_ratio=1.2)


class PhysicalLimitAudit(unittest.TestCase):
    def test_mcc_ellipse_allows_dry_state_above_critical_line(self):
        # This contradicts treating q<=Mp as a universal MCC admissibility rule.
        r = stress_values(state([100., 100., 0., 0., 0., 0.]), 1.4)
        self.assertAlmostEqual(r['p_kpa'], 200/3)
        self.assertAlmostEqual(r['q_kpa'], 100.)
        self.assertAlmostEqual(r['q_over_Mp'], 15/14)
        self.assertEqual(r['stress_side'], 'dry')
        self.assertLess(r['normalized_ellipse_yield'], 0.)
        self.assertGreater(r['q_over_Mp'], 1.)

    def test_tensor_q_keeps_shear_and_rotation_invariance(self):
        s = np.array([[80., 13., -7.], [13., 100., 11.], [-7., 11., 30.]])
        axis = np.array([1., 2., 3.]); axis /= np.linalg.norm(axis)
        cross = np.array([[0., -axis[2], axis[1]], [axis[2], 0., -axis[0]], [-axis[1], axis[0], 0.]])
        rot = np.eye(3)+np.sin(.7)*cross+(1-np.cos(.7))*(cross@cross)
        def voigt(a):
            return [a[0, 0], a[1, 1], a[2, 2], a[0, 1], a[1, 2], a[2, 0]]
        a, b = [stress_values(state(voigt(x)), 1.4) for x in (s, rot@s@rot.T)]
        self.assertAlmostEqual(a['p_kpa'], b['p_kpa'])
        self.assertAlmostEqual(a['q_kpa'], b['q_kpa'])
        self.assertAlmostEqual(stress_values(state([100., 100., 100., 10., 0., 0.]), 1.4)['q_kpa'], np.sqrt(300))

    def test_plastic_volume_uses_three_normal_components_only(self):
        original = state([80., 80., 80., 0., 0., 0.], plastic=[.01, -.02, .003, 7., 8., 9.])
        before = copy.deepcopy(original); r = stress_values(original, 1.4)
        self.assertAlmostEqual(r['plastic_volume_strain_compression_positive'], -.007)
        self.assertEqual(sign(r['plastic_volume_strain_compression_positive']), 'negative')
        self.assertEqual(original, before)

    def test_uniform_formula_and_self_weight_are_diagnostic_only(self):
        a = cone_lambda(1.4, 100., 100.)
        self.assertAlmostEqual(a['required_vertical_kpa'], 50/11)
        self.assertAlmostEqual(a['lambda_uniform'], 21/22)
        b = cone_lambda(1.4, 100., 110., 10.)
        self.assertAlmostEqual(b['lambda_self_weight_corrected'], 1-(50/11-10)/100)
        self.assertFalse(a['exact_bound_for_mcc'])
        self.assertEqual(cone_lambda(1.5, 100., 100.)['lambda_uniform'], 1.)

    def test_force_and_moment_are_additive_without_face_node_recounting(self):
        xyz = np.array([[1., 0., 0.], [0., 2., 0.]])
        f = np.array([[0., 0., -3.], [0., 0., -4.]])
        r = resultant(xyz, f)
        self.assertEqual(r['force_kn'], [0., 0., -7.])
        self.assertEqual(r['moment_about_origin_kn_m'], [-8., 3., 0.])
        self.assertEqual(face_label(np.array([[0., 0., .5], [1., 0., .5], [0., 1., .5]])), 'z=0.5')

    def test_nonfinite_stress_and_invalid_pressure_are_rejected(self):
        for value in ([float('nan')]*6, [0.]*6, [1., 2.]):
            with self.assertRaises(ValueError):
                stress_values(state(value), 1.4)
        with self.assertRaises(ValueError):
            stress_values(state([100., 100., 100., 0., 0., 0.]), 0.)

    def test_weak_mode_lift_preserves_recorded_arc_orientation(self):
        from make_free_contact import make
        from nc_physical_limit_audit import lift_recorded_mode
        mesh, _, _ = make(mcc=True, pc=100.)
        transform = np.zeros((len(mesh['nodes'])*3, 1)); transform[4*3+2, 0] = 1.
        geometry = dict(mesh=mesh, transform=transform.tolist(), active_elements=[True]*12)
        result = lift_recorded_mode(geometry, dict(progress=.9, tangent=[-1., 0.]),
                                    dict(right_weak_mode=[1.], stiffness={}, augmented={}))
        self.assertEqual(result['nodes'][4]['displacement_direction'], [0., 0., -1.])
        self.assertFalse(result['exact_null_vector'])
        self.assertFalse(result['augmented_right_null_vector_exists'])
        self.assertFalse(result['dynamic_stability_classified'])

    def test_changed_evidence_is_rejected_before_creating_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp); source = path/'source'; source.mkdir()
            (source/'evidence-sha256.json').write_text(json.dumps({'geometry.json': '0'*64}))
            (source/'geometry.json').write_text('{}')
            with self.assertRaisesRegex(ValueError, 'checksum differs'):
                run(source, path/'unread-origin.json', path/'output')
            self.assertFalse((path/'output').exists())


if __name__ == '__main__':
    unittest.main()
