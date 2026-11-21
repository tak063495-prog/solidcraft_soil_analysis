"""Algorithm tests use a scalar polynomial, never a soil reference solution."""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from continue_path import checkpoint, follow, restore, validate_controls


class Polynomial:
    kind = 'independent'
    R = np.eye(1)
    L = np.ones(1)
    old = np.zeros(1)
    u = np.zeros((1, 1))
    lam = .95
    t = np.array([1., .01]) / np.linalg.norm([1., .01])
    contacts = []

    def __init__(self):
        self.m = self

    def evaluate(self, du, old, label):
        x = old[0] + du[0]
        return np.array([x**3-.3*x*x+.01*x+.95]), np.array([[3*x*x-.6*x+.01]]), np.array([x])

    def load(self, lam):
        return np.array([lam])

    def serial(self, old):
        return old.tolist()

    def state(self, old):
        return np.array(old)

    def audit(self, *args):
        return {}


def rows(folder):
    return [json.loads(line) for line in (folder/'samples.jsonl').read_text().splitlines()]


class ContinuationTests(unittest.TestCase):
    def test_crosses_two_turns_and_lands_on_exact_endpoint_without_acceptance(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            result = follow(Polynomial(), out, .01, 100, 30., {'test': 1})
            samples = rows(out)
            self.assertTrue(any(r['progress_increment'] < 0 for r in samples))
            self.assertTrue(any(r['progress_increment'] > 0 for r in samples[25:]))
            self.assertEqual(samples[-1]['progress'], 1.)
            self.assertTrue(result['full_removal_reached'])
            self.assertFalse(result['passed'])
            self.assertFalse(result['static_limit_confirmed'])
            self.assertTrue(json.loads((out/'full-release.json').read_text())['zero_increment_material_exact'])

    def test_checkpoint_resume_reproduces_all_accepted_fields(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for name in ['whole', 'a', 'b']:
                (root/name).mkdir()
            binding = {'test': 2}
            follow(Polynomial(), root/'whole', .01, 100, 30., binding)
            follow(Polynomial(), root/'a', .01, 20, 30., binding)
            follow(Polynomial(), root/'b', .01, 80, 30., binding, root/'a/restart.json')
            self.assertEqual(rows(root/'whole'), rows(root/'a')+rows(root/'b'))

    def test_wrong_map_source_or_step_binding_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'cp.json'
            checkpoint(p, {'map': 'native', 'step': 1e-4}, {'state': [0.]})
            with self.assertRaisesRegex(ValueError, 'binding'):
                restore(p, {'map': 'independent', 'step': 1e-4})
            with self.assertRaisesRegex(ValueError, 'binding'):
                restore(p, {'map': 'native', 'step': 5e-5})

    def test_modified_history_payload_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'cp.json'
            checkpoint(p, {}, {'states': [1.]})
            cp = json.loads(p.read_text())
            cp['payload']['states'][0] = 2.
            p.write_text(json.dumps(cp))
            with self.assertRaisesRegex(ValueError, 'payload'):
                restore(p, {})

    def test_budget_exhaustion_keeps_inconclusive_flags(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = follow(Polynomial(), Path(tmp), .01, 100, 0., {})
            self.assertEqual(result['equilibrium_samples'], 0)
            self.assertIn('budget', result['stop_reason'])
            self.assertFalse(result['passed'])
            self.assertFalse(result['full_removal_reached'])
            self.assertFalse(result['static_limit_confirmed'])

    def test_strict_bounds_reject_nonfinite_and_unbounded_controls(self):
        for values in [(float('nan'), 1, 30), (0., 1, 30), (2e-4, 1, 30),
                       (1e-4, 4097, 30), (1e-4, True, 30), (1e-4, 1, float('inf'))]:
            with self.assertRaises(ValueError):
                validate_controls(*values)
        validate_controls(1e-4, 4096, 1800)


if __name__ == '__main__':
    unittest.main()
