"""Remeasure complete NC removal fields, including lossless actual systems."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT/'interim/2026-11-10-nc-mean-stress-release'
SPEC = importlib.util.spec_from_file_location('nc_mean_remeasure', HERE/'remeasure.py')
replay = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(replay)


class MeanStressEvidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(); cls.addClassCleanup(cls.temp.cleanup)
        cls.folder = Path(cls.temp.name)/'replay'
        cls.integrity = replay.run(cls.folder, allow_platform_float_differences=sys.platform == 'win32')
        cls.measured = json.loads((cls.folder/'measurements.json').read_text())
        print(json.dumps(dict(mean_stress_evidence_replay=cls.integrity), ensure_ascii=True))

    def test_archive_and_exact_raw_json_reconstruction_are_verified(self):
        self.assertTrue(self.integrity['archive_verified']); self.assertTrue(self.integrity['discrete_evidence_exact'])
        if sys.platform != 'win32':
            self.assertTrue(self.integrity['measurements_byte_identical'])
        self.assertFalse(self.integrity['numerical_acceptance'])

    def test_exploratory_paths_reach_zero_holding_but_keep_extra_confinement_separate(self):
        for kind, suffix in [('native', ''), ('independent', '-rate')]:
            for step in (128, 256):
                name = f'{kind}-complement{step}{suffix}'; result = self.measured['paths'][name]
                self.assertTrue(result['full_release']); self.assertEqual(result['remaining_holding_norm_kn'], 0.)
                self.assertTrue(result['removed_material_and_contacts_frozen']); self.assertLess(result['maximum_reassembled_force_residual'], 1e-8)
                self.assertFalse(result['report']['passed'])
                self.assertFalse(result['force_decomposition_bounded'])
                self.assertGreater(result['maximum_extra_holding_component_kn'],0.)
                post = self.measured['post_full_release'][name]
                self.assertEqual(post['actual_zero_increment_checks'], 2)
                self.assertTrue(post['all_material_contact_displacement_fields_exact'])
        for name,result in self.measured['paths'].items():
            if name.startswith('bounded-'):
                self.assertTrue(result['force_decomposition_bounded'])
                self.assertTrue(result['removed_material_and_contacts_frozen'])
                self.assertFalse(result['report']['passed'])

    def test_actual_checkpoints_restart_to_completion_with_exact_fields_and_systems(self):
        for name in ('native-complement128', 'native-complement256', 'independent-complement128-rate'):
            result = self.measured['restart'][name]
            self.assertTrue(result['completion_verified']); self.assertTrue(result['entire_accepted_rows_exact'])
            self.assertTrue(result['actual_linear_systems_exact']); self.assertGreater(result['compared_actual_linear_systems'], 0)
        self.assertTrue(self.measured['bounded_withdrawal']['unbounded_auxiliary_force_paths_not_transferred_to_goal'])
        self.assertTrue(self.measured['fresh_native_prefix']['fresh_unload_entire_state_exact'])
        self.assertEqual(set(self.measured['independent_prehistories']), {'400', '800'})

    def test_forged_acceptance_is_rejected(self):
        folder = self.folder/'evidence/native-complement128'; path = folder/'report.json'; raw = path.read_bytes()
        try:
            report = json.loads(raw); report['passed'] = True; path.write_text(json.dumps(report))
            with self.assertRaisesRegex(ValueError, 'falsely certifies'):
                replay.measure.audit(folder)
        finally:
            path.write_bytes(raw)


if __name__ == '__main__':
    unittest.main()
