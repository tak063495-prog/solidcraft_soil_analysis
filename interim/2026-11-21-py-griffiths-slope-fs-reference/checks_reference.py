"""Explicitly invoked tests, outside existing tests/tools discovery and CI counts."""
import copy
import importlib
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import lem
import srm
import compare


class EquilibriumTests(unittest.TestCase):
    def test_ground_geometry(self):
        np.testing.assert_allclose(lem.ground([0,12,22,32]),[10,10,5,0])

    def test_arc_endpoints(self):
        b,_,_=lem.arc(8,30,2)
        np.testing.assert_allclose(b([8,30]),lem.ground([8,30]),atol=1e-12)

    def test_forbid_below_fixed_base(self):
        with self.assertRaises(ValueError): lem.slices(8,32,32,lambda x:np.zeros_like(x)-1)

    def test_slice_weight_exact_polygon(self):
        bot=lambda x: 10-(np.asarray(x)-4)/2
        # Whole sloping wedge between x4 crest and x24 surface.
        # Must meet both endpoints, so use a lower V-shaped surface instead.
        bot=lambda x: lem.ground(x)-2*np.sin(np.pi*(np.asarray(x)-4)/24)
        sl=lem.slices(4,28,100,bot)
        x=sl['x'];area=np.trapezoid(lem.ground(x)-bot(x),x)
        # split at corner makes terrain area exact; corner aligns with slice grid here only approximately.
        self.assertLess(abs(sum(sl['W'])/20-area),.01)

    def test_planar_force_analytic(self):
        alpha=.4;W=np.array([100.,100.,100.]);length=np.array([2.,2.,2.])
        # Parallel planar slice bases; independent translational wedge solution.
        sl=dict(W=W,l=length,s=np.full(3,np.sin(alpha)),c=np.full(3,np.cos(alpha)),
                x=np.arange(4)*2*np.cos(alpha),y=-np.arange(4)*2*np.sin(alpha),
                xb=(np.arange(3)+.5)*2*np.cos(alpha),yb=-(np.arange(3)+.5)*2*np.sin(alpha),
                xg=(np.arange(3)+.5)*2*np.cos(alpha))
        F=(10*sum(length)+sum(W)*np.cos(alpha)*np.tan(np.deg2rad(20)))/(sum(W)*np.sin(alpha))
        r,E,N,_=lem.propagate(sl,F,0)
        self.assertLess(abs(E[-1]),1e-10)
        np.testing.assert_allclose(N,W*np.cos(alpha),atol=1e-10)

    def test_complete_closure_without_clipping_tension(self):
        bot,_,_=lem.arc(8.3333333333,31.5,3)
        sl=lem.slices(8.3333333333,31.5,32,bot)
        for fun in ['constant','sine']:
            v=lem.complete(sl,fun)
            r,E,N,_=lem.propagate(sl,v['fs'],v['lambda_value'],fun)
            self.assertLess(max(abs(r)),1e-8)
            self.assertFalse(v['compression_only'])
            self.assertTrue(min(E)<0 or min(N)<0)

    def test_grid_counts_and_no_acceptance(self):
        r=lem.search(3,16)
        self.assertEqual(r['candidates'],27)
        self.assertIsNone(r['accepted'])
        self.assertIsNotNone(r['minima']['bishop'])

    def test_non_circular_count(self):
        r=lem.search(3,16,True)
        self.assertEqual(r['candidates'],81)
        self.assertNotIn('bishop',r['minima'])


class MaterialAndFEMTests(unittest.TestCase):
    def test_mesh_area_gravity_boundary(self):
        m=srm.mesh(3)
        self.assertAlmostEqual(sum(m['area']),220)
        np.testing.assert_allclose(m['load'].reshape(-1,2).sum(axis=0),[0,-4400],atol=1e-10)
        self.assertNotIn(0,m['free']);self.assertNotIn(1,m['free'])

    def test_affine_strain_patch(self):
        m=srm.mesh(2);xy=m['xy'];u=np.c_[.002*xy[:,0]+.003*xy[:,1],-.001*xy[:,1]].ravel()
        strain=np.einsum('eij,ej->ei',m['B'],u[m['dofs']])
        np.testing.assert_allclose(strain,np.tile([.002,-.001,0,.003],(len(strain),1)),atol=1e-15)

    def test_rigid_rotation(self):
        m=srm.mesh(2);xy=m['xy'];u=np.c_[-xy[:,1],xy[:,0]].ravel()
        np.testing.assert_allclose(np.einsum('eij,ej->ei',m['B'],u[m['dofs']]),0,atol=1e-14)

    def test_elastic_energy_and_shear(self):
        D=srm.elastic();self.assertTrue(min(np.linalg.eigvalsh(D))>0)
        self.assertAlmostEqual(D[3,3],100000/2.6)

    def test_hydrostatic_elastic(self):
        mat=srm.MohrCoulomb();s=np.array([[-100.,-100.,-100.,0.]])
        out,dp=mat.update(s);np.testing.assert_allclose(out,s);np.testing.assert_allclose(dp,0,atol=1e-15)

    def test_mc_return_admissibility_zero_dilation(self):
        mat=srm.MohrCoulomb();trial=np.array([[-500.,-20.,-100.,17.]])
        out,dp=mat.update(trial);s=out[0]
        ev=np.linalg.eigvalsh([[s[0],s[3],0],[s[3],s[1],0],[0,0,s[2]]])
        self.assertLessEqual(max(ev@mat.A.T-mat.k),1e-7)
        self.assertAlmostEqual(sum(dp[0,:3]),0,places=10)
        self.assertGreater(float(np.dot(out[0],dp[0])),0)

    def test_rotated_elastic_branch_has_exact_zero_plastic_increment(self):
        trial=np.array([[-100.,-110.,-105.,2.]])
        out,dp=srm.MohrCoulomb().update(trial)
        np.testing.assert_array_equal(out,trial)
        np.testing.assert_array_equal(dp,np.zeros_like(trial))

    def test_tresca_uniaxial_yield(self):
        mat=srm.MohrCoulomb(phi=0);out,dp=mat.update([[-20.,0.,0.,0.]])
        np.testing.assert_allclose(out,[[-20,0,0,0]],atol=1e-12)
        np.testing.assert_allclose(dp,0,atol=1e-12)

    def test_rotational_objectivity(self):
        mat=srm.MohrCoulomb();angle=.47;q=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
        a=np.diag([-500.,-20.]);b=q@a@q.T
        one=mat.update([[-500,-20,-100,0]])[0][0]
        two=mat.update([[b[0,0],b[1,1],-100,b[0,1]]])[0][0]
        expect=q@np.array([[one[0],one[3]],[one[3],one[1]]])@q.T
        np.testing.assert_allclose(two[[0,1,3]],expect[[0,1,0],[0,1,1]],atol=1e-10)

    def test_strength_reduction_tangent_not_angle_division(self):
        a=srm.MohrCoulomb(factor=2)
        f=np.arctan(np.tan(np.deg2rad(20))/2)
        self.assertAlmostEqual(a.k,10*np.cos(f))
        self.assertAlmostEqual(a.A[0,0],1+np.sin(f))

    def test_small_gravity_real_solve(self):
        r=srm.trial_solve(srm.mesh(1),.8,steps=2,seconds=3)
        self.assertEqual(r['status'],'equilibrated');self.assertLess(r['relative_residual'],1e-7)

    def test_timeout_is_unknown(self):
        r=srm.trial_solve(srm.mesh(1),1,seconds=1e-12)
        self.assertEqual(r['status'],'timeout');self.assertFalse(r['failure_certified'])

    def test_budget_overrun_during_evaluation_is_unknown(self):
        m=srm.mesh(1);m['load']*=0
        with patch.object(srm.time,'monotonic',side_effect=[0.,0.,2.,2.]):
            r=srm.trial_solve(m,1,steps=1,seconds=1.)
        self.assertEqual(r['status'],'timeout');self.assertFalse(r['failure_certified'])

    def test_iteration_cap_is_unknown(self):
        r=srm.trial_solve(srm.mesh(1),1,iterations=0)
        self.assertEqual(r['status'],'iteration_limit');self.assertFalse(r['failure_certified'])

    def test_bisection_does_not_certify_unknown(self):
        def fake(m,factor,**kw):return dict(factor=factor,status='equilibrated' if factor<=1.3 else 'timeout')
        with patch.object(srm,'trial_solve',fake):r=srm.bisect_srm(1,max_trials=4)
        self.assertEqual(len(r['trials']),4);self.assertTrue(r['trial_cap_reached'])
        self.assertIsNone(r['certified_fs_interval']);self.assertIsNone(r['accepted'])
        self.assertLessEqual(r['equilibrium_lower'],1.3);self.assertGreater(r['unresolved_upper'],1.3)


class ContractTests(unittest.TestCase):
    def setUp(self):
        self.m=dict(schema_version=1,model_id=compare.MODEL,case_sha256=compare.case_fingerprint(),units=compare.UNITS.copy(),
                    conventions=compare.CONVENTIONS.copy(),reduction=compare.REDUCTION.copy(),accepted=None,fs_tolerance=None,
                    trial_count=2,trials=[dict(factor=1.,status='equilibrated',gravity_fraction=1.,relative_residual=1e-9,residual_tolerance=1e-7),
                                          dict(factor=1.5,status='timeout')],factor_interval=[1.,1.5],
                    endpoint_statuses=['equilibrated','timeout'],interval_kind='computational_unresolved',field_state=None)

    def test_valid_but_not_physical_fs(self):
        compare.validate(self.m);r=compare.compare(self.m,self.m)
        self.assertIsNone(r['accepted']);self.assertEqual(r['fs_relation'],'判定不能')

    def test_pinned_case_required(self):
        self.m['case_sha256']='b'*64
        with self.assertRaises(ValueError):compare.validate(self.m)

    def test_unknown_point_alias_rejected(self):
        self.m['best_factor']=1.25
        with self.assertRaises(ValueError):compare.validate(self.m)

    def test_units_rejected(self):
        self.m['units']['stress']='Pa'
        with self.assertRaises(ValueError):compare.validate(self.m)

    def test_sign_rejected(self):
        self.m['conventions']['y']='down'
        with self.assertRaises(ValueError):compare.validate(self.m)

    def test_reversed_interval_rejected(self):
        self.m['factor_interval']=[1.5,1.]
        with self.assertRaises(ValueError):compare.validate(self.m)

    def test_collapsed_interval_rejected(self):
        self.m['factor_interval']=[1.,1.]
        with self.assertRaises(ValueError):compare.validate(self.m)

    def test_point_values_rejected_recursively(self):
        for key in compare.FORBIDDEN:
            with self.subTest(key=key):
                m=copy.deepcopy(self.m);m['nested']={key:1.25}
                with self.assertRaises(ValueError):compare.validate(m)

    def test_implicit_failed_upper_rejected(self):
        self.m['interval_kind']='certified'
        with self.assertRaises(ValueError):compare.validate(self.m)

    def test_omitted_trial_rejected(self):
        self.m['trial_count']=3
        with self.assertRaises(ValueError):compare.validate(self.m)

    def test_uncertified_failure_rejected(self):
        self.m['trials'][1]['status']='certified_failure'
        with self.assertRaises(ValueError):compare.validate(self.m)

    def test_nonfinite_and_boolean_rejected(self):
        for v in [float('nan'),float('inf'),True]:
            with self.subTest(value=v):
                m=copy.deepcopy(self.m);m['trials'][0]['factor']=v
                with self.assertRaises(ValueError):compare.validate(m)

    def test_false_equilibrium_rejected(self):
        self.m['trials'][0]['relative_residual']=1e-3
        with self.assertRaises(ValueError):compare.validate(self.m)

    def test_new_acceptance_rejected(self):
        self.m['accepted']=True
        with self.assertRaises(ValueError):compare.validate(self.m)

    def test_different_input_rejected(self):
        r=copy.deepcopy(self.m);r['case_sha256']='b'*64
        with self.assertRaises(ValueError):compare.compare(self.m,r)

    def test_json_fields(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'nodes.json';p.write_text(json.dumps(dict(columns=compare.NODE,rows=[[1,0,0,.001,-.002]])))
            np.testing.assert_allclose(compare.table(p,compare.NODE),[[1,0,0,.001,-.002]])

    def test_boolean_json_field_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'nodes.json';p.write_text(json.dumps(dict(columns=compare.NODE,rows=[[True,0,0,0,0]])))
            with self.assertRaises(ValueError):compare.table(p,compare.NODE)

    def test_duplicate_csv_ids_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'nodes.csv';p.write_text(','.join(compare.NODE)+'\n1,0,0,0,0\n1,0,0,0,0\n')
            with self.assertRaises(ValueError):compare.table(p,compare.NODE)

    def test_csv_common_point_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder);(p/'a.csv').write_text(','.join(compare.NODE)+'\n1,0,0,0,0\n');(p/'b.csv').write_text(','.join(compare.NODE)+'\n1,1,0,0,0\n')
            f=dict(factor=1.,gravity_fraction=1.,reference='stress_free_gravity_turn_on',plastic_measure='engineering_tensor_components',nodes_csv='a.csv',plastic_csv='p.csv')
            a=copy.deepcopy(self.m);a['field_state']=f;b=copy.deepcopy(a);b['field_state']['nodes_csv']='b.csv'
            with self.assertRaisesRegex(ValueError,'common IDs'):compare.compare(a,b,p,p)

if __name__=='__main__':unittest.main()
