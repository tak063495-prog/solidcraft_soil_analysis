"""Check a real nonlinear rejection and the subsequent reassembled equilibria."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

HERE = Path(__file__).resolve().parents[1] / 'interim/2026-11-22-impl-c01-c03-nonlinear'
spec = importlib.util.spec_from_file_location('nonlinear_c03_probe', HERE/'probe_c03.py')
study = importlib.util.module_from_spec(spec)
spec.loader.exec_module(study)


class NonlinearJointProposals(unittest.TestCase):
    def test_failed_open_increment_reaches_same_anchor_with_unchanged_force_gate(self):
        old = study.measure(False, False)
        new = study.measure(False, True)
        self.assertIn('independent line search', old['stop'])
        self.assertEqual(old['rejected'][0], new['rejected'][0])
        self.assertEqual(old['prefix_sha256'], new['prefix_sha256'])
        self.assertEqual(new['stop'], 'observation point')
        self.assertEqual(new['last']['progress'], old['rejected'][0]['target'])
        self.assertGreater(len(new['observed_phase']), 1)
        self.assertLess(new['maximum_reassembled_force'], 1e-10)
        self.assertFalse(new['numerical_acceptance'])

    def test_documents_keep_unachieved_decision_under_windows_line_endings(self):
        for name in ('README.ja.md', 'REPRODUCE.ja.md'):
            original = (HERE/name).read_text(encoding='utf-8')
            for ending in ('\n', '\r\n'):
                with self.subTest(name=name, ending=repr(ending)), tempfile.TemporaryDirectory() as tmp:
                    path = Path(tmp)/name
                    path.write_bytes(original.replace('\n', ending).encode())
                    self.assertTrue(path.read_text(encoding='utf-8').endswith('C-01・C-03 は未達、passed=false。\n'))
