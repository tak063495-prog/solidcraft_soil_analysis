"""Checks on real retained histories; no numerical acceptance values invented."""
import copy
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'2026-11-09-nc-beyond-fold'))
from continue_path import restore
from measure import stream, fields


class Evidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m = json.loads(Path(os.environ['NC_CONTINUATION_MEASUREMENTS']).read_text())
        cls.e = Path(os.environ['NC_CONTINUATION_EVIDENCE'])

    def test_original_history_and_unreached_acceptance_remain_explicit(self):
        for key in ['passed', 'history_converged', 'static_limit_confirmed', 'production_adopted',
                    'full_removal_reached', 'full_release_freeze_verified', 'restart_to_completion_verified']:
            self.assertIs(self.m[key], False)
        self.assertEqual(self.m['native_highest_accepted_progress'], .9514000415802002)
        self.assertEqual(self.m['native_accepted_increments'], 266)
        self.assertEqual(self.m['native_new_accepted_increments'], 0)
        self.assertEqual(self.m['condition']['preload_kpa'], 5)
        self.assertEqual(self.m['condition']['pc_kpa'], 100)

    def test_additional_main_interval_has_4096_real_equilibria_per_map(self):
        for kind in ['native', 'independent']:
            r = self.m['paths'][f'extended-{kind}-next']
            self.assertEqual(r['recounted_samples'], 4096)
            self.assertEqual(r['negative_progress_increments'], 4096)
            self.assertEqual(r['positive_progress_increments'], 0)
            self.assertAlmostEqual(r['arc_length_m'], 1.3312, places=10)
            self.assertEqual(r['stop_reason'], 'bounded step count; inconclusive')
            self.assertFalse(r['full_removal_reached'])

    def test_force_and_actual_solved_linear_systems_keep_existing_gates(self):
        for r in self.m['paths'].values():
            self.assertLess(r['maximum_reassembled_relative_force_residual'], 1e-8)
            self.assertLessEqual(r['maximum_recomputed_true_linear_residual'], 1e-10)
            self.assertEqual(r['linear_systems_above_1e_10'], 0)
            self.assertEqual(r['recounted_samples'], r['equilibrium_samples'])

    def test_both_perturbed_zero_load_solves_are_real_and_still_not_full_release(self):
        for r in self.m['remaining_region_lambda_zero_perturbed_solves'].values():
            self.assertTrue(r['diagnostic_equilibrium_solved'])
            self.assertGreater(r['initial_increment_norm_m'], 0.)
            self.assertGreater(r['initial_relative_force_residual'], 1e-8)
            self.assertGreaterEqual(r['newton_evaluations'], 2)
            self.assertEqual(r['actual_linear_systems'], r['newton_evaluations'])
            self.assertLess(r['force_reassembly'], 1e-8)
            self.assertLessEqual(r['maximum_recomputed_true_linear_residual'], 1e-10)
            self.assertTrue(r['removed_material_frozen'])
            self.assertTrue(r['contact_history_frozen'])
            self.assertFalse(r['passed'])
            self.assertFalse(r['full_release_freeze_verified'])
        self.assertIn('identity', self.m['remaining_region_lambda_zero_independent']['note'])

    def test_every_saved_state_retains_all_3d_fields_and_frozen_removed_histories(self):
        for name, r in self.m['paths'].items():
            self.assertTrue(r['removed_material_frozen'])
            self.assertTrue(r['contact_history_frozen'])
            self.assertEqual(r['active_contact_points'], 0)
            native = '-native-' in name
            rows = stream(self.e/name/'samples.jsonl')
            first = next(rows)
            for row in [first]:
                self.assertEqual(np.asarray(row['displacement_m']).shape, (54, 3))
                self.assertEqual(np.asarray(fields(row, native)['stress']).shape, (324, 6))
                self.assertEqual(len(row['contacts']['points' if native else 'jump_m']), 32)

    def test_half_step_covers_entire_additional_interval_without_certifying_pre_history(self):
        for r in self.m['step_sensitivity'].values():
            self.assertEqual(r['matched_coarse_samples'], 4096)
            self.assertEqual(r['fine_samples'], 8192)
            self.assertTrue(r['same_parent_full_fields'])
            self.assertFalse(r['full_history_converged'])
            self.assertAlmostEqual(r['parent_arc_length_m'], .9216, places=10)
            self.assertEqual(set(r['endpoint']['material']), {'stress', 'plastic_strain', 'pc', 'void_ratio'})
            self.assertTrue(all(np.isfinite(v) for v in r['path_maximum_differences'].values()))

    def test_actual_runtime_checkpoint_restarts_match_every_saved_field_and_audit(self):
        for r in self.m['restarts'].values():
            self.assertEqual(r['whole_rows'], 128)
            self.assertEqual(r['restarted_rows'], 128)
            self.assertTrue(r['all_saved_fields_and_audits_exact'])
            self.assertTrue(r['capture']['actual_runtime_checkpoint'])
            self.assertTrue(r['capture']['source_report_absent_before_and_after'])
            self.assertFalse(r['production_job_restart'])
            self.assertFalse(r['restart_to_completion_verified'])

    def test_rebinding_or_payload_changes_cannot_resume(self):
        p = self.e/'parents/native.json'
        c = json.loads(p.read_text())
        self.assertEqual(restore(p, c['binding']), c['payload'])
        altered = copy.deepcopy(c['binding']); altered['material_steps'] = 8
        with self.assertRaises(ValueError): restore(p, altered)
        with tempfile.TemporaryDirectory() as d:
            q = Path(d)/'changed.json'; changed = copy.deepcopy(c)
            changed['payload']['progress'] += .001
            q.write_text(json.dumps(changed))
            with self.assertRaises(ValueError): restore(q, c['binding'])


if __name__ == '__main__': unittest.main()
