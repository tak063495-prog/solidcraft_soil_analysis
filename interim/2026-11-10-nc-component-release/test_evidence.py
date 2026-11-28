"""Published diagnostic evidence contracts; none confers physical acceptance."""
import json
import unittest
from pathlib import Path

HERE=Path(__file__).resolve().parent


class Evidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.m=json.loads((HERE/'measurements.json').read_text())

    def test_diagnostic_never_certifies_failure_or_completion(self):
        for k in ['passed','history_converged','static_limit_confirmed','full_removal_accepted','production_adopted']:
            self.assertFalse(self.m[k])
        self.assertEqual(self.m['native_committed_increments'],0)

    def test_two_component_orders_and_both_maps_are_present(self):
        self.assertEqual(len(self.m['paths']),10)
        for kind in ['native','independent']:
            for order in ['horizontal','vertical']:
                for step in [128,256]:self.assertIn(f'{kind}-{order}{step}',self.m['paths'])

    def test_all_saved_stress_fields_meet_original_force_gate(self):
        for p in self.m['paths'].values():
            self.assertGreater(p['saved_fields'],1)
            self.assertLess(p['maximum_reassembled_force_residual'],1e-8)

    def test_removed_material_and_contacts_remain_frozen(self):
        for p in self.m['paths'].values():self.assertTrue(p['removed_material_and_contacts_frozen'])

    def test_real_plastic_fields_and_actual_linear_systems_are_saved(self):
        for p in self.m['paths'].values():
            self.assertGreater(p['active_plastic_points_changed'],0)
            self.assertGreater(p['linear_systems'],p['saved_fields'])

    def test_half_steps_have_fullfield_comparisons_without_new_tolerance(self):
        self.assertEqual(len(self.m['step_sensitivity']),4)
        for p in self.m['step_sensitivity'].values():
            self.assertGreater(p['common_coordinates'],1)
            self.assertEqual(set(p['maximum_relative_L2']),{'displacement','stress','plastic','pc','void'})
            self.assertFalse(p['complete_history_converged'])

    def test_separate_process_restarts_preserve_every_common_field(self):
        self.assertEqual(len(self.m['restart']),2)
        for p in self.m['restart'].values():
            self.assertEqual(p['exact_common_full_fields'],65)
            self.assertFalse(p['completion_verified'])

    def test_no_missing_full_release_is_promoted_to_pass(self):
        for p in self.m['paths'].values():
            if not p['full_release']:
                self.assertIsNone(p['full_release_freeze'])
                self.assertGreater(p['remaining_holding_norm_kn'],0)
            self.assertFalse(p['report']['passed'])
            self.assertFalse(p['report']['full_release_acceptance'])


if __name__=='__main__':unittest.main()
