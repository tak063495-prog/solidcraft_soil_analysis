"""Ensure the alternative withdrawal never introduces extra nodal holding."""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from make_free_contact import make
from bounded_mean_release_control import partition, load


class BoundedMeanReleaseControl(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mesh, _, _ = make(mcc=True, pc=100., preload_kpa=5., force=(1., .5, 0.))
        upper = next(r['id'] for r in cls.mesh['regions'] if r['body'] == 'UpperClay')
        cls.active = np.array(cls.mesh['region_ids']) != upper
        cls.stress = np.tile([100., 100., 100., 0., 0., 0.], (27*int(cls.active.sum()), 1))
        cls.holding = np.sin(np.arange(3*len(cls.mesh['nodes']))); cls.holding[::7] = 0.

    def test_every_intermediate_component_is_bounded_by_original_signed_holding(self):
        for order in ('mean-first', 'complement-first'):
            a, b = partition(self.mesh, self.active, self.stress, self.holding, order)
            for s in np.linspace(0., 1., 129):
                f = load(np.zeros_like(a), self.holding, a, b, s)
                self.assertTrue(np.all(f*self.holding >= 0.))
                self.assertTrue(np.all(abs(f) <= abs(self.holding)+2e-16))
                np.testing.assert_array_equal(f[self.holding == 0.], 0.)

    def test_holding_parts_and_order_preserve_original_endpoints(self):
        a, b = partition(self.mesh, self.active, self.stress, self.holding, 'mean-first')
        np.testing.assert_allclose(a+b, self.holding, atol=1e-16, rtol=0.)
        c, d = partition(self.mesh, self.active, self.stress, self.holding, 'complement-first')
        np.testing.assert_array_equal(a, d); np.testing.assert_array_equal(b, c)
        final = np.cos(self.holding)
        np.testing.assert_array_equal(load(final,self.holding,a,b,0.),final+self.holding)
        np.testing.assert_array_equal(load(final,self.holding,a,b,1.),final)

    def test_original_zero_holding_cannot_generate_a_confining_load_loop(self):
        zero = np.zeros_like(self.holding)
        a, b = partition(self.mesh, self.active, self.stress, zero, 'mean-first')
        np.testing.assert_array_equal(a, 0.); np.testing.assert_array_equal(b, 0.)

    def test_invalid_order_is_rejected(self):
        with self.assertRaises(ValueError):
            partition(self.mesh, self.active, self.stress, self.holding, 'left-first')


if __name__ == '__main__':
    unittest.main()
