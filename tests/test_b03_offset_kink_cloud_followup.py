"""Light saved-evidence contracts; no patch solves or online retrieval in CI."""
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

from diagnostic_source import diagnostic_bytes
from workflow_pin import prior_record_sha256

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT/'interim/2026-12-27-b03-offset-kink-cloud-followup'
_path = sys.path[:]
try:
    sys.path.insert(0, str(HERE))
    from b03follow_evidence import verify
    from b03follow_math import arm_map, document_text, dot, guard, mv, norm
finally:
    sys.path[:] = _path


class B03OffsetKinkCloudFollowup(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((HERE/'measurements.json').read_text(encoding='utf-8'))

    def test_compensated_c_order_reduction_preserves_cancelled_term(self):
        self.assertEqual(dot([[1e16, 1.], [-1e16, 3.]], [[1., 1.], [1., 1.]]), 4.)
        self.assertEqual(norm([[3., 4.], [0., 12.]]), 13.)
        with self.assertRaises(ValueError):
            dot([math.inf], [1.])

    def test_surface_force_couple_and_displacement_are_work_conjugate(self):
        import numpy as np
        a = arm_map([0., 0., .125], 1)
        d = np.array([.125, .25, -.125, .03125, -.0625, .015625])
        fp = np.array([2., -3., 1., .25, .125, -.5])
        fs = mv(a.T, fp)
        self.assertEqual(fs[3:].tolist(), [-.125, -.125, -.5])
        self.assertEqual(dot(fp, mv(a, d)), dot(fs, d))
        self.assertNotEqual(dot(fp, d), dot(fs, d))

    def test_current_thickness_and_nontransitive_chain_are_distinct(self):
        tau = 64*sys.float_info.epsilon*.1
        a = [[0., 0., 0.], [0., 0., .75*tau], [0., 0., 1.5*tau]]
        self.assertTrue(guard(a, [.1]*3)['accepted'])
        self.assertFalse(guard([a[i] for i in (0, 2, 1)], [.1]*3)['accepted'])
        self.assertTrue(guard([a[0], a[2]], [.1, .2])['accepted'])
        self.assertFalse(guard([a[2], a[0]], [.2, .1])['accepted'])

    def test_boundary_endpoints_and_all_face_orders_are_present(self):
        b = self.data['boundary']
        self.assertEqual(b['actual_model_constructions'], 5859)
        self.assertEqual(b['prediction_mismatches'], [])
        self.assertEqual(b['common_edge_node_order_sensitive_groups'], 0)
        self.assertGreater(b['face_order_sensitive_sampled_groups'], 0)
        expected = {'same': [409, 204, 204, 204, 204, 409],
                    'mixed': [204, 102, 102, 204, 204, 409],
                    'binary-exact': [512, 256, 256, 256, 256, 512]}
        for profile in b['profiles']:
            rows = profile['face_order_intervals']
            self.assertEqual([r['maximum_accepted_gap_ulp'] for r in rows], expected[profile['profile']])
            self.assertTrue(all(r['last_accepted_in_actual_sweep'] and r['first_rejected_in_actual_sweep'] for r in rows))

    def test_permitted_rotated_patch_inventory_and_rigid_modes(self):
        r = self.data['rotation']
        self.assertEqual((r['native_solves'], r['independent_solves'], r['max_tri6'], r['max_dofs']), (45, 45, 18, 294))
        self.assertEqual(len(r['full_patch_results']['cases']), 45)
        for case in r['full_patch_results']['cases']:
            for side in ('native', 'independent'):
                self.assertEqual(case[side]['algebraic_zero_modes'], 6)
                self.assertEqual(case[side]['additional_zero_modes'], 0)
        self.assertFalse(r['full_patch_results']['physical_reference_solution_certified'])

    def test_original_float64_exceedances_are_retained(self):
        r = self.data['rotation']
        self.assertEqual(r['original_rotation_absolute_comparison_magnitude_rad'], 1e-9)
        self.assertEqual(len(r['exceedances_against_original_absolute_magnitude']), 13)
        self.assertTrue(all(v['maximum_absolute_rotation_difference_rad'] > 1e-9 for v in r['exceedances_against_original_absolute_magnitude']))
        self.assertTrue(r['exceedances_not_erased_by_decimal_diagnostic'])
        self.assertFalse(self.data['passed'])

    def test_decimal_control_changes_arithmetic_without_new_physical_cases(self):
        p = self.data['rotation']['selected_decimal_results']
        self.assertEqual(p['new_physical_equilibrium_cases'], 0)
        self.assertEqual(p['decimal_saddle_solves'], 8)
        self.assertTrue(p['all_50_80_digit_roundtrips_identical'])
        for row in p['cases']:
            self.assertLess(row['decimal_rotation_covariance_max_abs_rad'], row['ordinary_covariance']['rotation_max_abs_rad'])
        self.assertFalse(p['passed'])

    def test_work_account_does_not_claim_missing_soil_or_external_allocation(self):
        w = self.data['work']
        self.assertEqual((w['free_saved_native'], w['free_saved_independent'], w['coupled_saved_runs'],
                          w['coupled_unique_stage_end_states'], w['coupled_operator_accounts']), (23, 23, 6, 30, 60))
        self.assertEqual(w['new_equilibrium_solves'], 0)
        self.assertFalse(w['cumulative_s11_work_recomputed'])
        self.assertFalse(w['free_external_face_allocation_identifiable'])
        self.assertIsNone(w['rejected_junction_physical_work_primary_cause'])
        self.assertLess(w['free_maximum_full_account_residual_kn_m'], w['free_maximum_omitted_arm_account_residual_kn_m'])
        self.assertLess(w['coupled_maximum_full_account_residual_kn_m'], w['coupled_maximum_naive_force_dot_d_minus_twice_u_kn_m'])

    def test_all_symptoms_and_undone_rows_remain_explicit(self):
        self.assertEqual(set(self.data['contributions']), {'S'+str(i) for i in range(1, 13)})
        self.assertTrue(all(v is False for v in self.data['other_rows_complete'].values()))
        a = json.loads((HERE/'acceptance.json').read_text(encoding='utf-8'))
        for key in ('passed', 'b03_complete', 'production_scope_expanded', 'default_enabled'):
            self.assertIs(a[key], False)
        self.assertIsNone(a['general_joint_acceptance'])
        self.assertFalse(self.data['disposition']['pc_required_for_current_three'])

    def test_unacquired_literature_is_not_promoted_to_reference(self):
        r = json.loads((HERE/'primary-source-review.json').read_text(encoding='utf-8'))
        self.assertEqual(r['new_acceptance_reference_values'], [])
        self.assertEqual(r['new_primary_fulltext_confirmations'], 0)
        self.assertFalse(r['macneal_harder']['body_obtained'])
        self.assertTrue(all(not v.get('body_obtained', False) for v in r['supplemental_retrieval_attempts']))

    def test_old_records_and_production_sources_are_unchanged(self):
        # The follow-up records its production and reference code at this base commit.
        commit = '672824272e2baa1fe6ea776c8fc5889eacfc9a15'
        for name, expected in json.loads((HERE/'prior-record-binding.json').read_text()).items():
            data = diagnostic_bytes(ROOT, commit, name)
            self.assertEqual(prior_record_sha256(name, data), expected, name)

    def test_document_suffix_accepts_actual_crlf(self):
        suffix = 'B-03は未達、passed=false。\n'
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'README.ja.md'
            path.write_bytes(('観測\r\n'+suffix.replace('\n', '\r\n')).encode())
            self.assertTrue(document_text(path).endswith(suffix))
        self.assertTrue(document_text(HERE/'README.ja.md').endswith(suffix))

    def test_seal_and_sparse_checkout_keep_evidence_without_rerunning_solves(self):
        result = verify(HERE)
        self.assertTrue(result['all_sha256_verified'])
        self.assertEqual(result['numerical_reruns'], 0)
        workflow = (ROOT/'.github/workflows/fem-ci.yml').read_text()
        self.assertNotIn('!/interim/2026-12-27-b03-offset-kink-cloud-followup/', workflow)
        self.assertIn('/*', workflow)
        self.assertIn("matrix.group == 'python-c01-cause'", workflow)


if __name__ == '__main__':
    unittest.main()
