"""Check lossless accepted-field evidence without promoting numerical acceptance."""
import copy
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT/'interim/2026-11-08-wet-joint-finer-full-history'
sys.path.insert(0, str(HERE))
spec = importlib.util.spec_from_file_location('wet_joint_finer_study', HERE/'study.py')
study = importlib.util.module_from_spec(spec)
spec.loader.exec_module(study)
from wet_joint_cycle_reference import cut_volume_reference, serializable


class FinerFullHistory(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.values = study.inputs()
        cls.temporary = tempfile.TemporaryDirectory()
        cls.folder = Path(cls.temporary.name)/'accepted'
        sink = study.HistoryWriter(cls.folder, chunk_size=3)
        cls.observed = study.cut_volume_reference(*cls.values, 4, 8, 7, history_sink=sink)
        cls.manifest = sink.close()
        cls.original = cut_volume_reference(*cls.values, 4, 8, 7)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_observer_leaves_all_original_results_bitwise_unchanged(self):
        self.assertEqual(serializable(self.observed), serializable(self.original))

    def test_every_accepted_state_is_saved_and_endpoint_bits_match(self):
        report = study.check_history(self.folder, self.observed, self.values)
        self.assertEqual(report['accepted_states'], 28)
        self.assertEqual(report['fields_per_state'], 18)
        self.assertTrue(report['endpoint_fields_bitwise_equal'])
        self.assertFalse(report['full_history_acceptance'])

    def test_codec_preserves_random_bits_signed_zero_and_extremes(self):
        bits = np.random.default_rng(20261108).integers(0, 2**64, (7, 4, 6), dtype=np.uint64)
        bits.ravel()[:5] = [0, 1 << 63, 1, 0x7fefffffffffffff, 0x7ff8000000000001]
        actual = study.decode(study.encode(bits.view('<f8')))
        np.testing.assert_array_equal(actual.view('<u8'), bits)

    def test_common_accepted_comparison_is_a_measurement_without_acceptance(self):
        report = study.compare_accepted(self.folder, self.folder, ratio=1)
        self.assertEqual(report['common_accepted_states'], 28)
        self.assertFalse(report['full_history_acceptance'])
        for row in report['fields'].values():
            self.assertEqual(row['maximum_absolute'], 0.)
            self.assertEqual(row['pooled_relative_l2'], 0.)

    def test_doubled_path_compares_only_aligned_accepted_progress(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)/'fine'
            sink = study.HistoryWriter(folder, chunk_size=5)
            study.cut_volume_reference(*self.values, 8, 8, 7, history_sink=sink)
            sink.close()
            report = study.compare_accepted(self.folder, folder, ratio=2)
            self.assertEqual(report['common_accepted_states'], 28)
            self.assertFalse(report['full_history_acceptance'])
            self.assertFalse(report['finer_unmatched_states_compared'])
            with self.assertRaisesRegex(AssertionError, 'misaligned history'):
                study.compare_accepted(self.folder, folder, ratio=1)

    def test_modified_or_missing_chunk_is_rejected(self):
        import shutil
        for change in ('modify', 'missing'):
            with tempfile.TemporaryDirectory() as temp:
                folder = Path(temp)/'accepted'; shutil.copytree(self.folder, folder)
                target = folder/self.manifest['chunks'][0]['file']
                if change == 'modify':
                    target.write_bytes(target.read_bytes()+b'changed')
                else:
                    target.unlink()
                with self.assertRaises((AssertionError, FileNotFoundError)):
                    study.check_history(folder, self.observed, self.values)

    def test_missing_accepted_row_and_changed_schedule_are_rejected(self):
        for change in ('missing', 'schedule'):
            bad = copy.deepcopy(self.observed)
            if change == 'missing':
                bad['path'].pop()
            else:
                bad['path'][0]['factor'] += .1
            with self.assertRaises(AssertionError):
                study.check_history(self.folder, bad, self.values)

    def test_interior_nonfinite_or_coulomb_violation_is_rejected(self):
        system = study.Water(*self.values).contact
        fields = {key: np.asarray([self.observed['frames'][0][key]]) for key in study.FIELDS}
        bad = copy.deepcopy(fields); bad['stress'][0, 0, 0] = np.nan
        with self.assertRaises(AssertionError):
            study.check_fields(bad, system)
        bad = copy.deepcopy(fields); bad['traction'][0, 0, 0] += 1e6
        with self.assertRaises(AssertionError):
            study.check_fields(bad, system)


if __name__ == '__main__':
    unittest.main()
