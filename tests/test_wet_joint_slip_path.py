"""Diagnostics must preserve the solver and expose cancellation, not pass it."""
import copy
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT/'interim/2026-11-08-wet-joint-slip-path-diagnosis'
sys.path.insert(0, str(HERE))
spec = importlib.util.spec_from_file_location('slip_diagnostics', HERE/'slip_study.py')
study = importlib.util.module_from_spec(spec)
spec.loader.exec_module(study)
from wet_joint_cycle_reference import cut_volume_reference, serializable
import history_analysis


class WetJointSlipPath(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.values = study.inputs()
        cls.original = cut_volume_reference(*cls.values, 4, 2, 7)
        cls.observed = []
        cls.linear = []
        cls.instrumented = study.cut_volume_reference(*cls.values, 4, 2, 7,
            lambda row, fields: cls.observed.append((row, fields)), cls.linear.append)

    def test_observation_keeps_every_numeric_result_exact(self):
        self.assertEqual(serializable(self.original), serializable(self.instrumented))
        self.assertEqual(len(self.observed), 28)
        self.assertEqual(len(self.original['frames']), 7)
        for index, frame in enumerate(self.instrumented['frames']):
            row, fields = self.observed[(index+1)*4-1]
            self.assertEqual(row['progress'], 1.)
            for key in study.OBSERVED:
                np.testing.assert_array_equal(np.asarray(frame[key]).view('u8'), fields[key].view('u8'))

    def test_true_linear_and_force_residuals_keep_original_limits(self):
        self.assertGreater(len(self.linear), 0)
        self.assertTrue(all(np.isfinite(r['true_relative_residual']) and
                            r['true_relative_residual'] <= 1e-10 for r in self.linear))
        self.assertLess(self.original['maximum_relative_force_residual'], 1e-10)
        self.assertEqual(self.values[2]['settings']['relative_force_tolerance'], 1e-8)
        self.assertEqual(self.values[2]['settings']['relative_linear_tolerance'], 1e-10)

    def test_zero_difference_still_does_not_certify_history(self):
        result = study.pair(self.original, self.original)
        self.assertFalse(result['passed'] or result['history_converged'])
        self.assertEqual(len(result['plastic_slip_detail']), 7)
        for row in result['plastic_slip_detail']:
            self.assertEqual(row['numerator_l2_m'], 0.)
            self.assertEqual(len(row['all_points']), 32)

    def test_changed_slip_is_exposed_even_with_identical_displacement(self):
        other = copy.deepcopy(self.original)
        other['frames'][4]['plastic_slip'][0, 0] += 1e-6
        result = study.pair(other, self.original)
        row = result['plastic_slip_detail'][4]
        self.assertGreater(row['numerator_l2_m'], 0.)
        self.assertEqual(row['all_points'][0]['point'], 0)
        self.assertEqual(result['fields'][4]['u']['maximum_absolute'], 0.)
        self.assertFalse(result['history_converged'])

    def test_nonfinite_and_missing_frames_are_rejected(self):
        bad = copy.deepcopy(self.original)
        bad['frames'][4]['plastic_slip'][0, 0] = np.nan
        with self.assertRaises(AssertionError):
            study.pair(bad, self.original)
        bad = copy.deepcopy(self.original)
        bad['frames'].pop()
        with self.assertRaises(AssertionError):
            study.pair(bad, self.original)

    def test_contact_observation_codec_is_bitwise_lossless(self):
        with tempfile.TemporaryDirectory() as directory:
            sink = study.Sink(Path(directory)/'accepted')
            for row, fields in self.observed:
                sink(row, fields)
            sink.flush()
            self.assertEqual(sink.count, 28)
            self.assertEqual(len(sink.chunks), 7)
            first = sink.chunks[0]
            self.assertEqual(first['sha256'], study.sha(sink.folder/first['file']))
            with np.load(sink.folder/first['file'], allow_pickle=False) as saved:
                for key in study.OBSERVED:
                    actual = study.decode(saved[key])
                    expected = np.stack([f[key] for _, f in self.observed[:4]])
                    np.testing.assert_array_equal(actual.view('u8'), expected.view('u8'))

    def test_common_history_and_partial_prefix_are_distinguished(self):
        whole = history_analysis.compare(iter(self.observed), iter(self.observed))
        prefix = history_analysis.compare(iter(self.observed), iter(self.observed), 7)
        self.assertEqual(whole['common_accepted_states'], 28)
        self.assertFalse(whole['partial_prefix_only'] or whole['history_converged'])
        self.assertEqual(prefix['common_accepted_states'], 7)
        self.assertTrue(prefix['partial_prefix_only'])
        self.assertFalse(prefix['history_converged'])

    def test_changed_or_missing_common_history_is_rejected(self):
        bad = copy.deepcopy(self.observed)
        bad[0][0]['progress'] += 1e-6
        with self.assertRaisesRegex(AssertionError, 'changed path'):
            history_analysis.compare(iter(self.observed), iter(bad))
        with self.assertRaisesRegex(AssertionError, 'missing matching'):
            history_analysis.compare(iter(self.observed), iter(self.observed[:-1]))


if __name__ == '__main__':
    unittest.main()
