"""Checks on actual saved 3D diagnostics; no acceptance/stability certification."""
import copy
import json
import os
import tempfile
import unittest
from pathlib import Path

import numpy as np
from branch_path import restore


class FurtherEvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m = json.loads(Path(os.environ['NC_FURTHER_MEASUREMENTS']).read_text())
        cls.e = Path(os.environ['NC_FURTHER_EVIDENCE'])

    def test_production_history_and_inconclusive_flags(self):
        m = self.m
        for key in ['passed', 'history_converged', 'static_limit_confirmed', 'production_adopted',
                    'full_removal_reached', 'full_release_freeze_verified', 'restart_to_completion_verified']:
            self.assertIs(m[key], False)
        self.assertEqual(m['native_highest_accepted_progress'], .9514000415802002)
        self.assertEqual(m['native_accepted_increments'], 266)
        self.assertEqual(m['native_new_accepted_increments'], 0)
        self.assertEqual(m['windows_ci_for_new_head'], 'unconfirmed')

    def test_both_extended_segments_have_real_saved_points_and_no_return(self):
        for kind in ['native', 'independent']:
            for label in ['long', 'further']:
                r = self.m['paths'][f'extended-{kind}-{label}']
                self.assertEqual(r['equilibrium_samples'], 4096)
                self.assertEqual(r['negative_progress_increments'], 4096)
                self.assertEqual(r['positive_progress_increments'], 0)
                self.assertEqual(r['stop_reason'], 'bounded step count; inconclusive')
            self.assertAlmostEqual(self.m['paths'][f'extended-{kind}-further']['arc_length_m'], .9216, places=11)

    def test_force_and_actual_linear_systems_are_remeasured(self):
        for r in self.m['paths'].values():
            if r['recounted_samples']:
                self.assertLess(r['maximum_reassembled_relative_force_residual'], 1e-8)
            self.assertLessEqual(r['maximum_recomputed_true_linear_residual'], 1e-10)
            self.assertEqual(r['linear_systems_above_1e_10'], 0)
        for r in self.m['local_root_searches'].values():
            self.assertLess(r['maximum_reassembled_relative_force_residual'], 1e-8)
            self.assertLessEqual(r['maximum_recomputed_true_linear_residual'], 1e-10)
        for kind in ['native', 'independent']:
            self.assertLess(self.m[f'remaining_region_lambda_zero_{kind}']['relative_force_residual'], 1e-8)

    def test_removed_material_and_contacts_frozen_in_every_saved_field(self):
        for r in list(self.m['paths'].values())+list(self.m['local_root_searches'].values()):
            self.assertTrue(r['removed_material_frozen'])
            self.assertTrue(r['contact_history_frozen'])
        for r in self.m['paths'].values():
            self.assertEqual(r['recounted_samples'], r['equilibrium_samples'])
            self.assertEqual(r['active_contact_points'], 0)

    def test_reverse_probes_keep_real_histories_and_stop_inconclusive(self):
        for kind in ['native', 'independent']:
            parent = json.loads((self.e/f'parents/{kind}-whole.json').read_text())['payload']
            initial = json.loads((self.e/f'extended-{kind}-reverse/initial.json').read_text())['payload']
            for key in ['states', 'contacts', 'displacement_m', 'progress']:
                self.assertEqual(initial[key], parent[key])
            np.testing.assert_array_equal(initial['tangent'], -np.array(parent['tangent']))
            r = self.m['paths'][f'extended-{kind}-reverse']
            self.assertEqual(r['equilibrium_samples'], 0)
            self.assertEqual(r['stop_reason'], 'Newton limit; inconclusive')

    def test_same_anchor_half_step_is_only_partial_sensitivity(self):
        for kind in ['native', 'independent']:
            r = self.m['step_sensitivity'][kind+'-localhalf']
            self.assertTrue(r['same_full_field_start'])
            self.assertEqual(r['matched_points'], 512)
            self.assertEqual(r['fine_points'], 1024)
            self.assertFalse(r['full_extension_sensitivity_verified'])
            self.assertEqual(set(r['endpoint']['material']), {'stress', 'plastic_strain', 'pc', 'void_ratio'})
            self.assertFalse(self.m['step_sensitivity'][kind+'-fine1024']['same_full_field_start'])

    def test_local_root_inventory_does_not_certify_absence_or_full_release(self):
        self.assertEqual(len(self.m['local_root_searches']), 4)
        for r in self.m['local_root_searches'].values():
            self.assertEqual(r['trials'], 42)
            self.assertEqual(r['distinct_local_root_clusters'], 1)
            self.assertEqual(r['nonzero_increment_clusters'], 0)
            self.assertFalse(r['full_removal_root_found'])
            self.assertIn('inconclusive', r['stop_reason'])
            self.assertFalse(r['static_limit_confirmed'])

    def test_real_checkpoint_restarts_match_all_saved_fields(self):
        for kind, count in [('native', 128), ('independent', 64)]:
            r = self.m['restarts'][kind]
            self.assertEqual(r['compared_samples'], count)
            self.assertEqual(r['resumed_samples'], count)
            self.assertTrue(r['all_saved_fields_exact'])
            self.assertFalse(r['production_job_restart'])
            self.assertFalse(r['restart_to_full_removal_verified'])

    def test_zero_increment_predictors_trace_unloading_without_certifying_full_removal(self):
        for kind in ['native', 'independent']:
            parent = json.loads((self.e/f'parents/{kind}-whole.json').read_text())['payload']
            for label in ['unloading', 'reloading']:
                initial = json.loads((self.e/f'extended-{kind}-{label}/initial.json').read_text())['payload']
                for key in ['states', 'contacts', 'displacement_m', 'progress']:
                    self.assertEqual(initial[key], parent[key])
                r = self.m['paths'][f'extended-{kind}-{label}']
                self.assertLessEqual(r['zero_increment_predictor_system']['recomputed_true_relative_residual'], 1e-10)
                self.assertFalse(r['zero_increment_predictor_system']['stability_evidence'])
                self.assertEqual(r['diagnostic_verdict'], 'inconclusive')
            self.assertEqual(self.m['paths'][f'extended-{kind}-unloading']['equilibrium_samples'], 1024)
            r = self.m['paths'][f'extended-{kind}-reloading']
            self.assertEqual(r['equilibrium_samples'], 0)
            self.assertIn('strain bound', r['stop_reason'])
            self.assertFalse(r['full_removal_reached'])

    def test_checkpoint_binding_and_payload_tampering_are_rejected(self):
        p = self.e/'parents/independent-whole.json'
        original = json.loads(p.read_text())
        self.assertEqual(restore(p, original['binding']), original['payload'])
        wrong = copy.deepcopy(original['binding'])
        wrong['material_steps'] = 8
        with self.assertRaises(ValueError):
            restore(p, wrong)
        with tempfile.TemporaryDirectory() as d:
            q = Path(d)/'changed.json'
            changed = copy.deepcopy(original)
            changed['payload']['progress'] += .001
            q.write_text(json.dumps(changed))
            with self.assertRaises(ValueError):
                restore(q, original['binding'])


if __name__ == '__main__':
    unittest.main()
