"""Software and evidence checks, never physical acceptance thresholds."""
import json,os,sys,unittest
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));sys.path.insert(0,str(HERE.parents[1]/'tools'))
from native_arc import tangent
from independent_branches import StrainRate,expected_inputs
import mcc_rate_fem as rate

class EquationChecks(unittest.TestCase):
    def setup_material(self):
        mesh,ground,case=expected_inputs();state=rate.seed(mesh,ground)[0,:1].copy();p=ground['materials'][0]['mcc']
        parameters=np.array([[p['critical_state_m']**2,p['lambda'],p['kappa'],p['poisson']]])
        return state,StrainRate(parameters)
    def test_augmented_tangent_crosses_singular_load_control_fold(self):
        K=np.diag([0.,2.]);L=np.array([1.,0.]);previous=np.array([1.,0.,0.]);t=tangent(K,L,previous)
        np.testing.assert_array_equal(t,[1.,0.,0.]);np.testing.assert_array_equal(K@t[:-1]-L*t[-1],[0.,0.])
    def test_independent_elastic_hydrostatic_ode_closed_form(self):
        state,material=self.setup_material();old=state.copy();ev=-1e-4;strain=np.array([[ev/3]*3+[0.]*3]);new=material.integrate(state,strain,8)
        v=np.exp(state[0,7]);kap=material.parameters[0,2];p=state[0,:3].mean()*np.exp((v-v*np.exp(-ev))/kap)
        np.testing.assert_allclose(new[0,:3],p,rtol=1e-9,atol=0.);self.assertAlmostEqual(new[0,7],state[0,7]-ev,places=14)
        np.testing.assert_array_equal(new[:,8:],state[:,8:]);np.testing.assert_array_equal(state,old)
    def test_associated_nc_plastic_power_is_nonnegative(self):
        state,material=self.setup_material();strain=np.array([[1e-5/3]*3+[0.]*3]);_,dy,_,_=material.rate(state,strain,np.ones(1,dtype=bool))
        _,_,normal,_,_=rate.terms(state,material.parameters);gamma=float(dy[:,8:].ravel()@normal.ravel()/(normal.ravel()@normal.ravel()));power=float(state[:,:6].ravel()@dy[:,8:].ravel())
        self.assertGreater(gamma,0);self.assertGreater(power,0);p=state[0,:3].mean();expected=gamma*material.parameters[0,0]*p*np.exp(state[0,6]);self.assertAlmostEqual(power/expected,1.,places=13)
    def test_zero_strain_retains_all_independent_history(self):
        state,material=self.setup_material();new=material.integrate(state,np.zeros((1,6)),4);np.testing.assert_array_equal(new,state)

class EvidenceChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        p=Path(os.environ['NC_ARC_MEASUREMENTS']);cls.m=json.loads(p.read_text())
    def test_accepted_native_progress_and_false_acceptance_flags_preserved(self):
        m=self.m;self.assertEqual(m['native_accepted_increments'],266);self.assertEqual(m['native_highest_accepted_progress'],.9514000415802002);self.assertEqual(m['native_new_accepted_increments'],0)
        for key in ['passed','history_converged','production_adopted','native_full_removal_reached','full_release_freeze_verified','restart_to_completion_verified']:self.assertIs(m[key],False)
        self.assertEqual(m['windows_ci_for_new_head'],'unconfirmed')
    def test_force_reassembly_uses_existing_force_gate(self):
        self.assertLess(self.m['remaining_region_lambda_zero_native']['relative_force_residual'],1e-8)
        for row in self.m['paths'].values():
            if row.get('maximum_reassembled_relative_force_residual') is not None:self.assertLess(row['maximum_reassembled_relative_force_residual'],1e-8)
    def test_both_material_maps_have_local_turning_corrector_paths(self):
        for name in ['native-coarse','native-fine','independent-map-arc-coarse','independent-map-arc-fine']:
            row=self.m['paths'][name];self.assertIsNone(row['stop_reason']);self.assertGreater(row['positive_progress_increments'],0);self.assertGreater(row['negative_progress_increments'],0)
    def test_local_roots_are_distinct_and_not_stability_acceptance(self):
        for name in ['native-branch-search','independent-branch-search4','independent-branch-search8']:
            s=self.m['local_root_searches'][name];self.assertEqual(s['distinct_local_root_clusters_at_1e_8_m'],2);r=s['roots'][-1]
            self.assertGreater(r['increment_norm_m'],1e-4);self.assertLess(r['audit']['symmetric_eigen_min_kn_m'],0);self.assertGreater(r['audit']['plastic_stress_work_kn_m'],0);self.assertFalse(s['passed'])
    def test_diagnostic_restarts_reproduce_every_saved_field(self):
        for r in self.m['restarts'].values():self.assertTrue(r['all_saved_fields_bitwise_equal']);self.assertEqual(r['whole_samples'],r['resumed_samples']);self.assertFalse(r['restart_to_full_removal_verified'])
    def test_partial_histories_keep_removed_material_and_contact_frozen(self):
        for row in self.m['paths'].values():
            if row.get('equilibrium_samples',0):self.assertTrue(row['removed_material_contact_frozen'])
        h=self.m['partial_history'];self.assertEqual(h['accepted_rows'],h['recounted_rows']);self.assertTrue(h['terminal_and_completed_frames_exact']);self.assertTrue(h['partial_excavation_removed_material_frozen']);self.assertTrue(h['partial_excavation_contact_history_frozen'])
    def test_failed_refinement_and_remeasured_linear_exceedance_remain_visible(self):
        r=self.m['paths']['independent-arc-history400-half'];self.assertEqual(r['equilibrium_samples'],0);self.assertIsNotNone(r['stop_reason']);self.assertFalse(r['passed'])
        self.assertEqual(self.m['linear_system_audit']['independent-arc-extended400']['recomputed_above_existing_1e_10'],2)
    def test_saved_corrector_linear_systems_are_independently_audited(self):
        for name in ['native-coarse','native-fine','independent-map-arc-coarse','independent-map-arc-fine','independent-map-arc-material8']:
            r=self.m['linear_system_audit'][name];self.assertGreater(r['saved_systems'],0);self.assertLess(r['maximum_recomputed_true_relative_residual'],1e-10)

if __name__=='__main__':unittest.main()
