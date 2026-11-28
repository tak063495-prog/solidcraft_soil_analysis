"""Polynomial completeness and signed work checks for the independent C-02 reference."""
from itertools import product
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "interim/2026-11-25-approved-six"))
from c02_geometry import lagrange, simplex_rule
from c02_reference import project


class IndependentReference(unittest.TestCase):
    def test_tetrahedral_values_and_gradients_reproduce_polynomials(self):
        bary, weights = simplex_rule(3, 4)
        self.assertAlmostEqual(sum(weights), 1/6, places=14)
        for degree in (2, 3):
            alpha = np.array([a for a in product(range(degree+1),repeat=4) if sum(a)==degree])
            nodes = alpha/degree
            basis, gradient = lagrange(bary, alpha, degree)
            for exponent in (1, degree):
                nodal = (nodes[:,1] - 2*nodes[:,2] + 3*nodes[:,3])**exponent
                field = bary[:,1] - 2*bary[:,2] + 3*bary[:,3]
                np.testing.assert_allclose(basis@nodal, field**exponent, atol=3e-14)
                physical = gradient[:,:,1:] - gradient[:,:,:1]
                computed = np.einsum("qij,i->qj", physical, nodal)
                expected = (exponent*field**(exponent-1))[:,None]*np.array([1,-2,3])
                np.testing.assert_allclose(computed, expected, atol=3e-14)

    def test_stick_slide_reverse_open_and_free_flight(self):
        law = SimpleNamespace(normal=np.array([[0.,0.,1.]]),kn=np.array([10.]),
                              kt=np.array([4.]),mu=np.array([0.5]))
        old = np.array([[0.,0.,-1.]])
        traction = np.array([[0.,0.,-10.]])
        steps = [
            ([0.5,0,-1], [2,0,-10], [0,0,0], 1),
            ([2,0,-1], [5,0,-10], [0.75,0,0], 2),
            ([-1,0,-1], [-5,0,-10], [-0.5,0,0], 2),
            ([-3,0,1], [0,0,0], [-2.25,0,0], 0),
            ([3,2,1], [0,0,0], [0,0,0], 0),
        ]
        for target, expected, flow, branch in steps:
            target = np.array(target,dtype=float)[None,None,:]
            q, ds, status = project(law, old, traction, target, np.array([4]))
            np.testing.assert_allclose(q[0,0],expected,atol=1e-14)
            np.testing.assert_allclose(ds[0,0],flow,atol=1e-14)
            self.assertEqual(int(status[0,0]),branch)
            self.assertGreaterEqual(float(q[0,0]@ds[0,0]),0.)
            old, traction = target[:,0], q[:,0]
