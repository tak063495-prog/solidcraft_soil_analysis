"""Independent virtual-work and endpoint checks for the diagnostic force split."""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from make_free_contact import make
from mean_stress_release_control import partition, load


class MeanStressReleaseControl(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mesh, _, _ = make(mcc=True, pc=100., preload_kpa=5., force=(1., .5, 0.))
        upper = next(r['id'] for r in cls.mesh['regions'] if r['body'] == 'UpperClay')
        cls.active = np.array(cls.mesh['region_ids']) != upper
        cls.stress = np.tile([12., 18., 30., 7., -2., 3.], (27*int(cls.active.sum()), 1))
        cls.holding = np.sin(np.arange(3*len(cls.mesh['nodes'])))

    def test_mean_force_has_analytical_affine_virtual_work(self):
        first, _ = partition(self.mesh, self.active, self.stress, self.holding, 'mean-first')
        xyz = np.array(self.mesh['nodes'])
        # Existing lower cuboid is 1 x 1 x .5 m; trace(grad u)=.02.
        virtual = xyz@np.diag([.01, -.02, .03])
        self.assertAlmostEqual(float(first@virtual.ravel()), -20.*.02*.5, places=12)
        self.assertLess(np.linalg.norm(first.reshape(-1, 3).sum(axis=0)), 1e-12)
        self.assertLess(np.linalg.norm(np.cross(xyz, first.reshape(-1, 3)).sum(axis=0)), 1e-12)

    def test_trace_free_stress_and_shear_do_not_change_mean_force(self):
        a, _ = partition(self.mesh, self.active, self.stress, self.holding, 'mean-first')
        stress = self.stress+np.array([3., -5., 2., 19., 23., -29.])
        b, _ = partition(self.mesh, self.active, stress, self.holding, 'mean-first')
        np.testing.assert_array_equal(a, b)

    def test_complement_contains_exact_original_vector_with_roundoff_recorded(self):
        a, b = partition(self.mesh, self.active, self.stress, self.holding, 'mean-first')
        np.testing.assert_allclose(a+b, self.holding, atol=4e-15, rtol=0)
        reverse = partition(self.mesh, self.active, self.stress, self.holding, 'complement-first')
        np.testing.assert_array_equal(reverse[0], b)
        np.testing.assert_array_equal(reverse[1], a)

    def test_both_endpoint_loads_are_exact_without_cancellation(self):
        final = np.cos(self.holding)
        for order in ('mean-first', 'complement-first'):
            a, b = partition(self.mesh, self.active, self.stress, self.holding, order)
            np.testing.assert_array_equal(load(final, self.holding, a, b, 0.), final+self.holding)
            np.testing.assert_array_equal(load(final, self.holding, a, b, 1.), final)
            np.testing.assert_array_equal(load(final, self.holding, a, b, .5), final+b)

    def test_mean_nodal_force_is_not_a_component_or_node_mask(self):
        a, _ = partition(self.mesh, self.active, self.stress, self.holding, 'mean-first')
        self.assertTrue(np.any((abs(a) > 1e-10) & (abs(a-self.holding) > 1e-10)))
        self.assertGreater(np.linalg.norm(a.reshape(-1, 3)[:, :2]), 0.)
        self.assertGreater(np.linalg.norm(a.reshape(-1, 3)[:, 2]), 0.)

    def test_invalid_field_order_or_progress_is_rejected(self):
        for stress in (self.stress[:-1], self.stress*float('nan')):
            with self.assertRaises(ValueError):
                partition(self.mesh, self.active, stress, self.holding, 'mean-first')
        with self.assertRaises(ValueError):
            partition(self.mesh, self.active, self.stress, self.holding, 'left-first')
        for s in (-.1, 1.1, float('nan')):
            with self.assertRaises(ValueError):
                load(self.holding, self.holding, self.holding, self.holding, s)


if __name__ == '__main__':
    unittest.main()
