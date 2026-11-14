"""Saved 3D diagnostic evidence checks; these are not acceptance tests."""
import json
import os
import unittest
from pathlib import Path


class EvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m = json.loads(Path(os.environ['NC_EXTENDED_MEASUREMENTS']).read_text())

    def test_inconclusive_flags_and_native_production_progress(self):
        m = self.m
        for key in ['passed', 'static_limit_confirmed', 'history_converged', 'production_adopted',
                    'full_removal_reached', 'full_release_freeze_verified', 'restart_to_completion_verified']:
            self.assertIs(m[key], False)
        self.assertEqual(m['native_accepted_increments'], 266)
        self.assertEqual(m['native_highest_accepted_progress'], .9514000415802002)
        self.assertEqual(m['native_new_accepted_increments'], 0)

    def test_both_full_3d_maps_extend_beyond_fold(self):
        for kind in ['native', 'independent']:
            r = self.m['paths'][f'extended-{kind}-whole']
            self.assertEqual(r['equilibrium_samples'], 1024)
            self.assertGreater(r['positive_progress_increments'], 0)
            self.assertGreater(r['negative_progress_increments'], 0)
            self.assertLess(r['last_progress'], r['maximum_sampled_progress'])
            self.assertGreater(r['positive_to_negative_eigen_crossings'], 0)
            self.assertLess(r['minimum_symmetric_eigenvalue_kn_m'], 0)
            self.assertEqual(r['stop_reason'], 'bounded step count; inconclusive')

    def test_existing_force_and_linear_gates_reassembled(self):
        for r in self.m['paths'].values():
            self.assertLess(r['maximum_reassembled_relative_force_residual'], 1e-8)
            self.assertLessEqual(r['maximum_recomputed_true_linear_residual'], 1e-10)
            self.assertEqual(r['linear_systems_above_1e_10'], 0)
        for kind in ['native', 'independent']:
            self.assertLess(self.m[f'remaining_region_lambda_zero_{kind}']['relative_force_residual'], 1e-8)

    def test_material_and_contact_history_frozen_at_every_saved_point(self):
        for r in self.m['paths'].values():
            self.assertTrue(r['removed_material_frozen'])
            self.assertTrue(r['contact_history_frozen'])
            self.assertEqual(r['active_contact_points'], 0)
            self.assertEqual(r['recounted_samples'], r['equilibrium_samples'])

    def test_same_arc_extent_and_full_field_step_comparisons_exist(self):
        for r in self.m['step_sensitivity'].values():
            self.assertTrue(r['same_requested_total_arc_length'])
            self.assertEqual(r['coarse_samples'], 1024)
            self.assertEqual(r['fine_samples'], 2048)
            self.assertEqual(r['matched_coarse_samples'], 1024)
            self.assertEqual(set(r['endpoint']['material']), {'stress', 'plastic_strain', 'pc', 'void_ratio'})
        self.assertFalse(self.m['history_converged'])

    def test_actual_both_map_checkpoints_reproduce_every_saved_field(self):
        for r in self.m['restarts'].values():
            self.assertTrue(r['all_saved_fields_exact'])
            self.assertEqual(r['continuous_samples'], 1024)
            self.assertEqual(r['split_samples'], 1024)
            self.assertFalse(r['production_job_restart'])
            self.assertFalse(r['restart_to_full_removal_verified'])


if __name__ == '__main__':
    unittest.main()
