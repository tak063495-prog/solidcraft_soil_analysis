"""Explicit new tests; deliberately outside existing test_*.py collection."""
import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from numpy.testing import assert_allclose
from scipy.spatial.transform import Rotation

import geometry_reference as g


class ClosedGeometry(unittest.TestCase):
    def test_kink_general_offsets(self):
        for theta in np.deg2rad([0, 15, 45, 90, 135, 179]):
            for e1, e2 in [(.05, .05), (-.05, .05), (0, .03), (.02, -.04)]:
                n2 = [np.sin(theta), 0, np.cos(theta)]
                got = g.mismatch([0, 0, 0], [0, 0, 1], e1, [0, 0, 0], n2, e2)
                assert_allclose(got['delta_m'], g.kink_closed(theta, e1, e2), rtol=0, atol=0)
                self.assertAlmostEqual(got['norm_m']**2, e1**2+e2**2-2*e1*e2*np.cos(theta), places=15)

    def test_equal_offset_chord(self):
        for theta in np.deg2rad([0, 45, 90, 179]):
            for e in [-.05, 0, .05]:
                self.assertAlmostEqual(np.linalg.norm(g.kink_closed(theta, e, e)), 2*abs(e)*np.sin(theta/2), places=15)

    def test_step_signed_vector(self):
        for h in [-.1, 0, .1]:
            for e1, e2 in [(.05, -.05), (-.05, .05), (0, 0)]:
                out = g.mismatch([0, 0, 0], [0, 0, 1], e1, [.03, 0, h], [0, 0, 1], e2)
                assert_allclose(out['delta_m'], g.step_closed(h, .03, e1, e2), atol=0, rtol=0)

    def test_signed_projections(self):
        theta, e1, e2 = .7, -.03, .05
        r = g.mismatch([0, 0, 0], [0, 0, 1], e1, [0, 0, 0], [np.sin(theta), 0, np.cos(theta)], e2)
        self.assertAlmostEqual(r['delta_dot_n1_m'], e2*np.cos(theta)-e1)
        self.assertAlmostEqual(r['delta_dot_n2_m'], e2-e1*np.cos(theta))

    def test_parallel_gap_overlap_touch(self):
        self.assertAlmostEqual(g.parallel_layer_gap(.2, .1, .1), .1)
        self.assertAlmostEqual(g.parallel_layer_gap(.05, .1, .1), -.05)
        self.assertEqual(g.parallel_layer_gap(.1, .1, .1), 0)
        self.assertAlmostEqual(g.parallel_layer_gap(0, .1, .2), -.1)

    def test_normal_and_offset_joint_flip_invariant(self):
        n = np.array([.6, 0, .8])
        a = g.mismatch([0, 0, 0], n, .04, [0, 0, .1], n, -.03)
        b = g.mismatch([0, 0, 0], -n, -.04, [0, 0, .1], -n, .03)
        self.assertEqual(a['delta_m'], b['delta_m'])
        self.assertEqual(a['norm_m'], b['norm_m'])

    def test_attachment_to_neutral_sign(self):
        n1, n2 = np.array([0, 0, 1]), np.array([1, 0, 0])
        e1, e2 = .02, .05
        xs = np.array([1, 2, 3])
        xp1, xp2 = xs-e1*n1, xs-e2*n2
        assert_allclose(xp2-xp1, -(e2*n2-e1*n1), atol=5e-16, rtol=0)

    def test_rotation_translation_covariance(self):
        r = Rotation.from_rotvec([.2, .3, -.4]).as_matrix()
        x = np.array([[0, 0, 0], [.03, .02, .1]])
        n = np.array([[0, 0, 1], [1, 0, 0]])
        e = np.array([.05, -.03])
        out = g.mismatch(x[0], n[0], e[0], x[1], n[1], e[1])
        y = x@r.T + [3, -2, 1]
        nn = n@r.T
        rot = g.mismatch(y[0], nn[0], e[0], y[1], nn[1], e[1])
        assert_allclose(rot['delta_m'], r@out['delta_m'], rtol=0, atol=5e-16)


class CandidateDiagnostics(unittest.TestCase):
    def setUp(self):
        self.x = np.zeros((2, 3))
        self.n = np.array([[0, 0, 1.], [1., 0, 0]])
        self.e = np.array([.05, .05])
        self.p, self.q = g.candidate_points(self.x, self.n, self.e)

    def test_independent_rigid_arms(self):
        d = g.diagnose_candidate(self.p, self.q['independent_rigid_arms'])
        self.assertEqual(d['max_geometry_error_m'], 0)
        self.assertEqual(d['original_rigid_subspace_dimension'], 6)
        self.assertEqual(d['constraint_nullity'], 6)
        self.assertEqual(d['original_work_defect_max_coefficient'], 0)
        self.assertIsNone(d['accepted'])

    def test_average_normal_geometry(self):
        error = np.linalg.norm(self.q['average_normal_common_point']-self.p, axis=1)
        assert_allclose(error, 2*.05*np.sin(np.pi/8), rtol=1e-14, atol=0)

    def test_bisector_projection_geometry(self):
        q = self.q['orthogonal_bisector_projection']
        assert_allclose(q, [[.025, 0, .025]]*2, atol=2e-17, rtol=0)
        assert_allclose(np.linalg.norm(q-self.p, axis=1), .05*np.sin(np.pi/4), rtol=1e-14)

    def test_moved_candidates_lose_original_rotations(self):
        for name in ['average_normal_common_point', 'orthogonal_bisector_projection']:
            d = g.diagnose_candidate(self.p, self.q[name])
            self.assertEqual(d['represented_rigid_map_rank'], 6)
            self.assertEqual(d['constraint_nullity'], 6)
            self.assertEqual(d['original_rigid_subspace_dimension'], 4)
            self.assertGreater(d['original_work_defect_max_coefficient'], .01)
            self.assertLess(d['represented_six_mode_max_residual'], 1e-15)

    def test_parallel_bisector_is_undefined(self):
        _, q = g.candidate_points(self.x, [[0, 0, 1]]*2, self.e)
        self.assertIsNone(q['orthogonal_bisector_projection'])
        self.assertFalse(g.diagnose_candidate(self.p, None)['defined'])

    def test_opposed_average_is_undefined(self):
        _, q = g.candidate_points(self.x, [[0, 0, 1], [0, 0, -1]], self.e)
        self.assertIsNone(q['average_normal_common_point'])

    def test_rigid_arm_work_all_basis_pairs(self):
        r = np.array([.2, -.3, .4])
        b = g.arm(r)
        for j in range(6):
            f = np.eye(6)[j]
            gen = np.r_[f[:3], np.cross(r, f[:3])+f[3:]]
            assert_allclose(b.T@f, gen, rtol=0, atol=0)
        assert_allclose(g.arm(-r)@b, np.eye(6), rtol=0, atol=0)

    def test_naive_equal_dof_negative_control(self):
        c = np.hstack([-np.eye(6), np.eye(6)])
        bp = np.vstack([g.arm(p) for p in self.p])
        self.assertEqual(g.numerical_rank(c@bp)[0], 2)
        self.assertGreater(np.max(np.abs(c@bp)), .01)

    def test_projection_unequal_offsets_does_not_promise_merge(self):
        _, q = g.candidate_points(self.x, self.n, [.02, .05])
        self.assertGreater(np.linalg.norm(np.diff(q['orthogonal_bisector_projection'], axis=0)), .01)

    def test_flat_identical_geometry_is_control(self):
        p, q = g.candidate_points(self.x, [[0, 0, 1]]*2, [.05, .05])
        d = g.diagnose_candidate(p, q['average_normal_common_point'])
        self.assertEqual(d['max_geometry_error_m'], 0)
        self.assertEqual(d['original_rigid_subspace_dimension'], 6)


class ContractChecks(unittest.TestCase):
    def setUp(self):
        self.doc = g.synthetic_document()

    def check_bad(self, mutator):
        bad = copy.deepcopy(self.doc)
        mutator(bad)
        with self.assertRaises((ValueError, TypeError)):
            g.compare_documents(self.doc, bad)

    def test_synthetic_json_csv_roundtrips(self):
        with tempfile.TemporaryDirectory() as td:
            for suffix in ['json', 'csv']:
                p = Path(td)/('fixture.'+suffix)
                g.write_document(self.doc, p)
                actual = g.read_document(p)
                self.assertEqual(actual, self.doc)
                result = g.compare_documents(self.doc, actual)
                self.assertTrue(result['contract_consistent'])
                self.assertIsNone(result['accepted'])

    def test_units_rejected(self):
        self.check_bad(lambda d: d['metadata'].update(length_unit='mm'))

    def test_basis_rejected(self):
        self.check_bad(lambda d: d['metadata'].update(coordinate_basis='local'))

    def test_normal_sign_metadata_rejected(self):
        self.check_bad(lambda d: d['metadata'].update(normal_sign='averaged'))

    def test_offset_sign_metadata_rejected(self):
        self.check_bad(lambda d: d['metadata'].update(offset_sign='neutral=attachment+e*normal'))

    def test_missing_id_rejected(self):
        self.check_bad(lambda d: d['rows'].pop())

    def test_duplicate_id_rejected(self):
        self.check_bad(lambda d: d['rows'].append(d['rows'][0].copy()))

    def test_missing_field_rejected(self):
        self.check_bad(lambda d: d['rows'][0].pop('neutral_m'))

    def test_nonfinite_rejected(self):
        self.check_bad(lambda d: d['rows'][0].update(e_m=float('nan')))

    def test_outside_e_rejected_without_slack(self):
        self.check_bad(lambda d: d['rows'][0].update(e_m=float(np.nextafter(.05, np.inf))))
        self.assertEqual(g.orientation_diagnostics([0, 0, 1], .05, .1, [0, 0, 1]), [])
        self.assertEqual(g.orientation_diagnostics([0, 0, 1], -.05, .1, [0, 0, 1]), [])
        self.assertIn('offset_outside_thickness', g.orientation_diagnostics([0, 0, 1], -.051, .1, [0, 0, 1]))

    def test_nonpositive_t_rejected(self):
        self.check_bad(lambda d: d['rows'][0].update(t_m=0))

    def test_normal_scaling_rejected(self):
        self.check_bad(lambda d: d['rows'][0].update(normal=[0, 0, 2]))

    def test_compensated_sign_flip_rejected(self):
        bad = copy.deepcopy(self.doc)
        bad['rows'][0].update(normal=[0, 0, -1], e_m=-.05)
        r = g.compare_documents(self.doc, bad)
        self.assertFalse(r['contract_consistent'])
        self.assertIn('normal_orientation_reversed', r['rows'][0]['errors'])
        self.assertIn('offset_sign_reversed', r['rows'][0]['errors'])

    def test_silent_average_rejected(self):
        bad = copy.deepcopy(self.doc)
        for row in bad['rows']:
            row['attachment_m'] = [.025, 0, .025]
        self.assertFalse(g.compare_documents(self.doc, bad)['contract_consistent'])

    def test_silent_clip_rejected(self):
        bad = copy.deepcopy(self.doc)
        bad['rows'][0]['e_m'] = .04
        bad['rows'][0]['attachment_m'] = [0, 0, .04]
        self.assertFalse(g.compare_documents(self.doc, bad)['contract_consistent'])

    def test_even_tiny_discrepancy_not_accepted(self):
        bad = copy.deepcopy(self.doc)
        bad['rows'][0]['attachment_m'][2] = float(np.nextafter(.05, np.inf))
        self.assertFalse(g.compare_documents(self.doc, bad)['contract_consistent'])

    def test_orientation_without_reference_unconfirmed(self):
        self.assertIn('normal_orientation_unconfirmed', g.orientation_diagnostics([0, 0, 1], 0, .1))

    def test_duplicate_json_key_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)/'bad.json'
            p.write_text('{"rows":[],"rows":[],"metadata":{}}')
            with self.assertRaises(ValueError):
                g.read_document(p)

    def test_csv_missing_coordinate_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td)/'bad.csv'
            g.write_document(self.doc, p)
            text = p.read_text().replace('normal_x,', '')
            p.write_text(text)
            with self.assertRaises(ValueError):
                g.read_document(p)


    def test_piecewise_face_constants_cannot_vary(self):
        bad = copy.deepcopy(self.doc)
        second = copy.deepcopy(bad['rows'][0])
        second.update(shared_node_id='corner-2', e_m=.04, attachment_m=[0, 0, .04])
        bad['rows'].append(second)
        with self.assertRaisesRegex(ValueError, 'constant'):
            g.validate_document(bad)

    def test_series_all_closed_forms_and_no_adoption(self):
        r = g.series()
        self.assertEqual(r['case_count'], 294)
        self.assertIsNone(r['accepted'])
        for case in r['cases']:
            self.assertLess(case['closed_vector_error_m'], 2e-15)
            self.assertLess(case.get('closed_norm_error_m', 0), 2e-15)
            self.assertIsNone(case['accepted'])

    def test_compare_cli_exit_codes(self):
        with tempfile.TemporaryDirectory() as td:
            ref, actual, out = [Path(td)/x for x in ['ref.json', 'actual.json', 'out.json']]
            g.write_document(self.doc, ref)
            g.write_document(self.doc, actual)
            cmd = [sys.executable, str(Path(g.__file__)), 'compare', str(ref), str(actual), str(out)]
            self.assertEqual(subprocess.run(cmd, capture_output=True).returncode, 2)
            bad = copy.deepcopy(self.doc)
            bad['rows'][0]['normal'] = [0, 0, -1]
            g.write_document(bad, actual)
            self.assertEqual(subprocess.run(cmd, capture_output=True).returncode, 1)
            self.assertIsNone(json.loads(out.read_text(encoding='utf-8'))['accepted'])



class ArchivedRefusals(unittest.TestCase):
    def test_archived_b03_geometry_and_original_refusals(self):
        from legacy_refusals import audit_repository
        root = Path(__file__).resolve().parents[2]
        report = audit_repository(root)
        self.assertTrue(all(r['same_rejection_confirmed'] for r in report['cases']))
        self.assertIsNone(report['accepted'])
        byname = {r['case']: r for r in report['cases']}
        self.assertIn('offset_outside_thickness', byname['outside_half_thickness']['independent_guard']['errors'])
        folded = byname['folded_offset_normals']['independent_guard']
        self.assertIn('shared_signed_normal_offset_mismatch', folded['errors'])
        self.assertTrue(all(abs(n['normal_dot']+1) < 1e-14 for n in folded['shared_nodes']))
        self.assertTrue(all(x['status'] == '未確認' for x in report['unconfirmed']))

    def test_archived_input_tamper_cannot_be_confirmed(self):
        from legacy_refusals import inspect_direct_geometry
        root = Path(__file__).resolve().parents[2]
        fixture = json.loads((root/'examples/plate6-offset-stages.json').read_text())
        original = copy.deepcopy(fixture)
        self.assertFalse(inspect_direct_geometry(fixture)['rejected'])
        fixture['beam_system']['plates'][0]['attachment_offset_m'] = .051
        self.assertTrue(inspect_direct_geometry(fixture)['rejected'])
        self.assertEqual(original['beam_system']['plates'][0]['attachment_offset_m'], .05)

if __name__ == '__main__':
    unittest.main(verbosity=2)
