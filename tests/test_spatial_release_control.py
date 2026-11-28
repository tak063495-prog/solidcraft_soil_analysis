"""Checks for diagnostic ordering; success never certifies NC full removal."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import time
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from make_free_contact import make
from spatial_release_control import fractions, partition, load

SPEC = importlib.util.spec_from_file_location('nc_spatial_diagnostic', ROOT/'interim/2026-11-10-nc-spatial-release/study.py')
study = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(study)


class SpatialReleaseControl(unittest.TestCase):
    def test_fixture_is_exact_existing_5_kpa_case(self):
        expected = make(mcc=True, pc=100., preload_kpa=5., force=(1., .5, 0.))
        for name, value in zip(('mesh', 'ground', 'stages'), expected):
            self.assertEqual(json.loads((ROOT/'tests/fixtures/low-preload-nc'/(name+'.json')).read_text()), value)
        self.assertNotEqual(expected[2], make(mcc=True, pc=100., preload_kpa=100., force=(1., .5, 0.))[2])

    def test_spatial_order_preserves_all_components_and_endpoints(self):
        nodes = np.array(make(mcc=True, pc=100., force=(1., .5, 0.))[0]['nodes'])
        h = np.arange(nodes.size, dtype=float)-17.
        final = np.sin(h)
        for strategy in ('left-first', 'center-first'):
            a, b, chosen = partition(nodes, h, strategy)
            np.testing.assert_array_equal(a+b, h)
            self.assertTrue(chosen)
            np.testing.assert_array_equal(load(final, a, b, 1.), final)
            np.testing.assert_array_equal(load(final, a, b, .5), final+b)
            np.testing.assert_allclose(load(final, a, b, 0.), final+h, rtol=0, atol=3e-14)
            # Every component follows its original signed holding contribution.
            for s in np.linspace(0., 1., 17):
                rem = load(np.zeros_like(final), a, b, s)
                self.assertTrue(np.all(rem*h >= 0))
                self.assertTrue(np.all(abs(rem) <= abs(h)))

    def test_spatial_partition_cannot_silently_become_horizontal_vertical(self):
        nodes = np.array(make(mcc=True, pc=100.)[0]['nodes'])
        h = np.ones(nodes.size)
        for strategy in ('left-first', 'center-first'):
            a, b, _ = partition(nodes, h, strategy)
            self.assertTrue(np.all(a.reshape(-1, 3)[:, 0] == a.reshape(-1, 3)[:, 2]))
            self.assertTrue(np.all(b.reshape(-1, 3)[:, 0] == b.reshape(-1, 3)[:, 2]))
        self.assertFalse(np.array_equal(partition(nodes, h, 'left-first')[0], partition(nodes, h, 'center-first')[0]))

    def test_invalid_coordinate_or_degenerate_groups_are_rejected(self):
        for s in (-.01, 1.01, float('nan'), float('inf')):
            with self.assertRaises(ValueError):
                fractions(s)
        with self.assertRaises(ValueError):
            partition([[0, 0, 0]], [1., 1., 1.], 'left-first')
        with self.assertRaises(ValueError):
            partition([[0, 0, 0], [1, 1, 1]], np.ones(5), 'center-first')

    def test_real_newton_corrects_perturbed_zero_release(self):
        class Elastic:
            R = np.ones((1, 1))
            B = np.zeros((1, 6, 1))
            def load(self, s):
                return np.array([s])
            def evaluate(self, du, old, label):
                return du.copy(), np.ones((1, 1)), old.copy()
        with tempfile.TemporaryDirectory() as folder:
            new, du, events, reason = study.solve(Elastic(), [1.], 0., Path(folder), 'zero', time.monotonic()+5, np.array([1e-4]))
            self.assertIsNone(reason)
            self.assertGreaterEqual(len(events), 2)
            self.assertGreater(events[0]['relative_force_residual'], 1e-8)
            self.assertEqual(new, [1.])
            np.testing.assert_array_equal(du, [0.])

    def test_failed_trial_does_not_mutate_committed_material_history(self):
        class NoRoot:
            R = np.ones((1, 1))
            B = np.zeros((1, 6, 1))
            def load(self, s):
                return np.array([1.])
            def evaluate(self, du, old, label):
                return np.zeros(1), np.ones((1, 1)), [99.]
        old = [2.]
        with tempfile.TemporaryDirectory() as folder:
            new, du, _, reason = study.solve(NoRoot(), old, 1., Path(folder), 'failed', time.monotonic()+5)
            self.assertIsNone(new)
            self.assertIsNone(du)
            self.assertIn('line search exhausted', reason)
            self.assertEqual(old, [2.])

    def test_time_budget_is_inconclusive_without_equilibrium(self):
        class Uncalled:
            R = np.ones((1, 1))
            def evaluate(self, *args):
                raise AssertionError('expired solve must not evaluate')
        with tempfile.TemporaryDirectory() as folder:
            new, _, events, reason = study.solve(Uncalled(), [], 1., Path(folder), 'expired', time.monotonic()-1)
            self.assertIsNone(new)
            self.assertEqual(events, [])
            self.assertEqual(reason, 'wall time budget; inconclusive')


if __name__ == '__main__':
    unittest.main()
