"""Refuse unsupported branch/holding-load claims; matrices here are algebra only."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
path=ROOT/'interim/2026-11-10-nc-postfold-branches/branch_math.py'
spec=importlib.util.spec_from_file_location('nc_branch_math',path)
branch=importlib.util.module_from_spec(spec)
spec.loader.exec_module(branch)


class BranchFollowup(unittest.TestCase):
    def test_signed_original_component_bounds_refuse_extra_load_and_nonfinite(self):
        h=np.array([2.,-3.,0.])
        branch.component_bounds(h,[1.,-2.,0.])
        for r in ([2.1,-2.,0.],[1.,.1,0.],[1.,-2.,1e-16],[1.,np.nan,0.]):
            with self.assertRaises(ValueError):branch.component_bounds(h,r)

    def test_turning_sign_is_not_a_bifurcation_or_stability_certificate(self):
        result=branch.diagnose(np.diag([1.,1e-12]),np.array([0.,1.]),np.array([0.,1.,-1e-6]),np.array([0.,1.,1e-6]))
        self.assertTrue(result['turning_candidate'])
        self.assertFalse(result['bifurcation_confirmed'])
        self.assertGreater(result['augmented_smallest_singular_kn_m'],.9)
        self.assertGreater(result['weak_left_load_alignment'],.9)

    def test_augmented_rank_loss_still_does_not_confirm_a_branch(self):
        result=branch.diagnose(np.diag([1.,0.]),np.array([1.,0.]),np.array([0.,1.,0.]),np.array([0.,1.,0.]))
        self.assertEqual(result['augmented_smallest_singular_kn_m'],0.)
        self.assertFalse(result['bifurcation_confirmed'])

    def test_parallel_mode_cannot_be_called_a_transverse_branch_seed(self):
        with self.assertRaisesRegex(ValueError,'parallel'):
            branch.probe_tangent([1.,0.,0.],[1.,0.],1.)
        for amplitude in [-1.,-.25,.25,1.]:
            seed=branch.probe_tangent([1.,0.,0.],[0.,1.],amplitude)
            self.assertAlmostEqual(np.linalg.norm(seed),1.)
            self.assertEqual(np.sign(seed[1]),np.sign(amplitude))


if __name__=='__main__':unittest.main()
