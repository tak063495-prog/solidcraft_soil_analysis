"""Replay saved 3D evidence, without treating diagnostic completion as acceptance."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT/'interim/2026-11-10-nc-spatial-release'
SPEC = importlib.util.spec_from_file_location('nc_spatial_remeasure', HERE/'remeasure.py')
replay = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(replay)


class SpatialEvidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temp.cleanup)
        cls.folder = Path(cls.temp.name)/'replay'
        cls.integrity = replay.run(cls.folder,allow_platform_float_differences=sys.platform=='win32')
        cls.measured = json.loads((cls.folder/'measurements.json').read_text())
        print(json.dumps(dict(spatial_evidence_replay=cls.integrity),ensure_ascii=True))

    def test_all_archive_bytes_and_recomputed_measurements_match(self):
        self.assertTrue(self.integrity['archive_verified'])
        self.assertTrue(self.integrity['discrete_evidence_exact'])
        if sys.platform!='win32':self.assertTrue(self.integrity['measurements_byte_identical'])
        self.assertFalse(self.integrity['numerical_acceptance'])

    def test_both_maps_and_orders_have_full_saved_fields_and_frozen_histories(self):
        for kind in ('native', 'independent'):
            for order in ('left', 'center'):
                for step in (128, 256):
                    result = self.measured['paths'][f'{kind}-{order}{step}']
                    self.assertGreater(result['saved_fields'], 1)
                    self.assertGreater(result['linear_systems'], 0)
                    self.assertTrue(result['removed_material_and_contacts_frozen'])
                    self.assertLess(result['maximum_reassembled_force_residual'], 1e-8)
                    self.assertFalse(result['full_release'])
                    self.assertFalse(result['report']['passed'])
        for key in ('native-left', 'native-center', 'independent-left', 'independent-center'):
            self.assertGreater(self.measured['step_sensitivity'][key]['common_coordinates'], 1)
            self.assertFalse(self.measured['step_sensitivity'][key]['complete_history_converged'])

    def test_separate_process_restarts_match_actual_corresponding_fields(self):
        for key in ('native-left', 'independent-center'):
            result = self.measured['restart'][key]
            self.assertEqual(result['exact_common_full_fields'], 13)
            self.assertFalse(result['completion_verified'])
        phases=self.measured['production_phase_restarts']
        self.assertEqual(len(phases['verification']['comparisons']),4)
        self.assertTrue(phases['remeasured_entire_run_except_elapsed_seconds'])
        self.assertFalse(phases['completion_verified'])

    def test_forged_acceptance_flag_is_rejected_by_field_audit(self):
        folder = self.folder/'evidence/independent-left128'
        path = folder/'report.json'; raw = path.read_bytes()
        try:
            report = json.loads(raw); report['passed'] = True
            path.write_text(json.dumps(report))
            with self.assertRaisesRegex(ValueError, 'falsely certifies'):
                replay.measure.audit(folder)
        finally:
            path.write_bytes(raw)


if __name__ == '__main__':
    unittest.main()
