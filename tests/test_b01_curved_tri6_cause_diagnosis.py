"""Small B-01 contracts. Sealed matrix/global solves are not repeated in CI."""
import hashlib
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT/'interim/2026-12-25-b01-curved-tri6-cause-diagnosis'
_path = sys.path[:]
try:
    sys.path.insert(0, str(HERE))
    from b01cause_build import generate
    from b01cause_cli_controls import clean_run
    from b01cause_evidence import verify
    from b01cause_geometry import PARAM, document_text, dot, frame, geometry, mesh, metric, norm, shape, surface
    from b01cause_report import prior_assessment
finally:
    sys.path[:] = _path


class B01CurvedTri6CauseDiagnosis(unittest.TestCase):
    def test_saved_run_time_exclusion_is_recursive_and_preserves_physics(self):
        a = dict(attempts=[dict(linear_solves=[dict(elapsed_seconds=1., residual=1e-11)])],
                 provenance=dict(solver_binary_sha256='old'), value=-0.)
        b = dict(attempts=[dict(linear_solves=[dict(elapsed_seconds=2., residual=1e-11)])],
                 provenance=dict(solver_binary_sha256='new'), value=-0.)
        self.assertEqual(json.dumps(clean_run(a), sort_keys=True), json.dumps(clean_run(b), sort_keys=True))
        b['attempts'][0]['linear_solves'][0]['residual'] = 2e-11
        self.assertNotEqual(clean_run(a), clean_run(b))
        self.assertEqual(a['attempts'][0]['linear_solves'][0]['elapsed_seconds'], 1.)

    def test_fixed_order_norm_and_original_normalization(self):
        self.assertEqual(dot([1e16, 1., -1e16], [1., 1., 1.]), 1.)
        self.assertEqual(norm([[3., 4.], [0., 12.]]), 13.)
        self.assertEqual(metric([1e-18], [0.])['relative_l2'], 1e-18/1e-12)
        with self.assertRaises(ValueError):
            metric([math.nan], [0.])
        with self.assertRaises(ValueError):
            metric([1.], [1., 2.])

    def test_quadratic_partition_and_interpolation_use_actual_geometry(self):
        for i, p in enumerate(PARAM):
            n, _, _ = shape(*p)
            np.testing.assert_array_equal(n, np.eye(6)[i])
        n, d, dd = shape(.21, .32)
        self.assertAlmostEqual(float(sum(n)), 1.)
        np.testing.assert_allclose(d.sum(axis=0), [0., 0.], atol=1e-15)
        np.testing.assert_allclose(dd.sum(axis=0), np.zeros((2, 2)), atol=1e-15)

    def test_normal_derivative_agrees_with_independent_central_difference(self):
        x = surface('sphere', PARAM); r, s, h = .24, .31, 1e-5
        f = frame(x, r, s)
        for axis, step in enumerate(((h, 0.), (0., h))):
            plus = frame(x, r+step[0], s+step[1])['normal']
            minus = frame(x, r-step[0], s-step[1])['normal']
            np.testing.assert_allclose(f['normal_derivative'][axis], (plus-minus)/(2*h), atol=2e-10, rtol=1e-8)
        self.assertGreater(norm(f['normal_derivative']), .1)

    def test_coplanarity_is_not_the_existing_straight_midpoint_contract(self):
        x = surface('plane', PARAM); x[3, 1] += .02
        a = geometry(x)
        self.assertEqual(max(map(abs, a['out_of_corner_plane_m'])), 0.)
        self.assertGreater(max(a['midpoint_distances_m']), a['guard_cutoff_m'])
        self.assertEqual(a['maximum_normal_derivative'], 0.)

    def test_mesh_boundary_and_resource_bound_are_fixed_before_solve(self):
        x, elements, boundary = mesh('sphere', 2)
        self.assertEqual((len(x), len(elements), len(boundary)), (25, 8, 16))
        self.assertEqual(len({v for e in elements for v in e}), 25)
        with self.assertRaises(ValueError):
            mesh('sphere', 8)
        self.assertLess(geometry(surface('sphere', PARAM))['exact_quadratic_rigid_covariant_membrane_max'], 1e-14)

    def test_transplant_is_work_only_and_has_only_reviewed_guard_omission(self):
        sources = [ROOT/'src/plates.rs', ROOT/'src/plates/mindlin.rs', ROOT/'src/plates/offset.rs']
        before = [p.read_bytes() for p in sources]
        with tempfile.TemporaryDirectory() as tmp:
            dest = Path(tmp); binding = generate(dest)
            diff = (dest/'guard-only.diff').read_text()
            self.assertIn('-        for (e, (i, j))', diff)
            self.assertIn('requires straight edges', diff)
            self.assertNotIn('-            bm[', diff)
            self.assertNotIn('-            bb[', diff)
            self.assertTrue(binding['production_guard_unchanged'])
        self.assertEqual(before, [p.read_bytes() for p in sources])

    def test_document_suffix_accepts_real_crlf(self):
        marker = 'B-01は未達、passed=false。\n'
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'README.ja.md'
            path.write_bytes(('記録\r\n'+marker.replace('\n', '\r\n')).encode())
            self.assertTrue(document_text(path).endswith(marker))
        self.assertTrue(document_text(HERE/'README.ja.md').endswith(marker))

    def test_prior_values_and_short_series_do_not_create_a_pass(self):
        bindings = json.loads((HERE/'prior-record-binding.json').read_text())
        for name, expected in bindings.items():
            self.assertEqual(hashlib.sha256((ROOT/name).read_bytes()).hexdigest(), expected, name)
        a = prior_assessment()
        self.assertIsNone(a['roof_observed_order'])
        self.assertTrue(a['native_minus_q4_sign_reversal'])
        self.assertEqual(a['legacy_references_unchanged']['open_top_hemisphere'], .094)
        self.assertFalse(a['reference_definitions_interchangeable'])
        self.assertIsNone(a['original_numerical_acceptance'])

    def test_sealed_evidence_and_false_status_are_separate(self):
        a = json.loads((HERE/'acceptance.json').read_text())
        self.assertFalse(a['passed'])
        self.assertFalse(a['b01_complete'])
        self.assertTrue(a['diagnostic_only'])
        self.assertIsNone(a['curved_element_acceptance'])
        self.assertFalse(a['production_scope_expanded'])
        result = verify(HERE)
        self.assertTrue(result['all_sha256_verified'])
        self.assertEqual(result['numerical_reruns'], 0)
        # Dedicated C-01 jobs exclude archives and never select this module;
        # ordinary core jobs retain every new B-01 evidence part.
        workflow = (ROOT/'.github/workflows/fem-ci.yml').read_text()
        self.assertIn("matrix.group == 'python-c01-cause'", workflow)
        self.assertNotIn('!/interim/2026-12-25-b01-curved-tri6-cause-diagnosis/', workflow)

    def test_algebraic_mode_count_does_not_override_physical_rigid_residual(self):
        a = json.loads((HERE/'measurements.json').read_text())
        data = {r['id']: r for r in a['main']}
        self.assertEqual(data['sphere']['inherited_zero_count'], 6)
        self.assertGreater(max(data['sphere']['physical_rigid_residual']), 1e-3)
        self.assertGreater(max(data['sphere']['physical_rigid_fraction_outside_algebraic_nullspace']), .01)
        self.assertTrue(data['outplane-bulge']['internal_K_exactly_equal_to_plane_same_corners'])
        self.assertEqual(a['global_equilibrium_solves'], 9)
        self.assertFalse(a['passed'])
        self.assertFalse(a['b01_complete'])


if __name__ == '__main__':
    unittest.main()
