"""The independent gate must reject a corrupted solution despite a zero solver report."""
import copy
import unittest
from audit_linear import audit

class IndependentLinearAudit(unittest.TestCase):
    def values(self):
        req=[dict(label='nonsymmetric',matrix_columns=[[4.,1.],[2.,3.]],reduced_residual=[1.,2.])]
        res=[dict(label='nonsymmetric',baseline_x=[-.1,.7],baseline_relative_residual=0.,
                  candidate_x=[-.1,.7],accepted=True,original_relative_residual=0.,
                  compensated_relative_residual=0.,corrections=0)]
        return req,res

    def test_nonsymmetric_orientation_and_acceptance_scope(self):
        req,res=self.values();r=audit(req,res)
        self.assertTrue(r['all_saved_matrices_pass_independent_linear_gate'])
        self.assertFalse(r['passed']);self.assertFalse(r['history_converged'])

    def test_corrupted_solution_is_not_certified_by_solver_report(self):
        req,res=self.values();res[0]['candidate_x'][0]+=.01
        self.assertFalse(audit(req,res)['all_saved_matrices_pass_independent_linear_gate'])

    def test_failed_solve_and_misassociated_matrix_do_not_pass(self):
        req,res=self.values();bad=copy.deepcopy(res);bad[0]['accepted']=False
        self.assertFalse(audit(req,bad)['all_saved_matrices_pass_independent_linear_gate'])
        res[0]['label']='another-matrix'
        with self.assertRaises(AssertionError):audit(req,res)

if __name__=='__main__':unittest.main()
