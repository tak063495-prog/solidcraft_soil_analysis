"""Original-equation gates for reused independent tangent factors."""
import os
for _name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[_name]='1'
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
import fem_reference as ref
import mcc_projected_sparse_reference as fresh
import mcc_projected_factor_reference as reused
import mcc_projected_matrixfree_reference as matrixfree
import mcc_projected_volume_diagnostic as initial
import mcc_projected_cycle_diagnostic as cycle


def specimen():
    xyz=np.vstack((np.zeros(3),np.eye(3)))
    xyz=np.vstack((xyz,[(xyz[a]+xyz[b])/2 for a,b in ref.EDGES]))
    angle=.37
    rotation=np.array([[np.cos(angle),-np.sin(angle),0.],[np.sin(angle),np.cos(angle),0.],[0.,0.,1.]])@np.array([[1.,0.,0.],[0.,np.cos(.29),-np.sin(.29)],[0.,np.sin(.29),np.cos(.29)]])
    xyz=xyz@rotation.T
    mesh=dict(nodes=xyz.tolist(),tetrahedra=[[0,1,2,3]],tetrahedra10=[list(range(10))],quadrature_order=3,
        material_ids=[1],node_sets={f'n{i}':[i] for i in range(4)})
    ground=dict(materials=[dict(id=1,mcc=dict(critical_state_m=1.2,**{'lambda':.2,'kappa':.04},poisson=.3))])
    constraints=[dict(node_set=f'n{n}',direction=rotation[:,a].tolist(),value_m=0.) for n,axes in [(0,[0,1,2]),(1,[1,2]),(2,[2])] for a in axes]
    case=dict(settings=dict(kinematics='explicit'),boundary=dict(constraints=constraints))
    return mesh,ground,case


def assert_rates(test,a,b):
    for left,right in zip(a,b,strict=True):np.testing.assert_allclose(left,right,rtol=1e-9,atol=1e-11)


class FactorReference(unittest.TestCase):
    candidate_factory=reused.ReusedFactorProjectedSystem
    def systems(self):
        mesh,ground,case=specimen()
        return fresh.SparseProjectedSystem(mesh,ground,[True],case),self.candidate_factory(mesh,ground,[True],case)

    def state(self,system,p=100.,pc=200.):
        state=np.zeros((system.weight.size,14));state[:,:3]=p;state[:,6:8]=np.log([pc,2.]);return state

    def test_coercive_rates_reuse_lu_and_retain_all_three_original_block_gates(self):
        a,b=self.systems();state=self.state(a);trial=state.copy();trial[:,:6]+=[1.,2.,3.,.4,.3,.2]
        force=a.force(trial)-a.force(state);plastic=np.zeros(len(state),bool)
        for i in range(12):
            candidate=state.copy();candidate[:,:3]+=i*.01
            assert_rates(self,a.rate(candidate,force,plastic),b.rate(candidate,force,plastic))
        self.assertGreater(b.cached_factor_hits,0)
        self.assertLess(b.actual_tangent_factorizations,a.factorizations)
        self.assertLessEqual(max(b.maximum_block_residuals.values()),1e-10)
        self.assertLessEqual(b.maximum_true_linear_residual,1e-10)
        self.assertLessEqual(len(b._factors),2)

    def test_large_coercive_tangent_change_refactors_the_current_original_matrix(self):
        a,b=self.systems();state=self.state(a);trial=state.copy();trial[:,:3]+=[1.,2.,3.]
        force=a.force(trial)-a.force(state);plastic=np.zeros(len(state),bool)
        b.rate(state,force,plastic)
        changed=self.state(a,800.,2000.)
        assert_rates(self,a.rate(changed,force,plastic),b.rate(changed,force,plastic))
        self.assertGreater(b.fresh_fallbacks+getattr(b,'matrix_free_fallbacks',0),0)
        self.assertLessEqual(max(b.maximum_block_residuals.values()),1e-10)

    def test_noncoercive_mcc_branch_uses_fresh_lu_instead_of_a_cached_factor(self):
        a,b=self.systems();state=self.state(a,80.,200.);q=np.sqrt(1.2**2*80.*120.)
        state[:,0]+=2*q/3;state[:,1:3]-=q/3
        trial=state.copy();trial[:,:6]+=[1.,2.,3.,.4,.3,.2]
        force=a.force(trial)-a.force(state);plastic=np.ones(len(state),bool)
        assert_rates(self,a.rate(state,force,plastic),b.rate(state,force,plastic))
        self.assertEqual(b.cached_factor_hits,0)
        self.assertEqual(b.noncoercive_bypasses,1)
        self.assertEqual(b.actual_tangent_factorizations,1)

    def test_complete_plastic_k0_cycle_matches_fresh_reference_without_changing_physical_gates(self):
        points=[[float(x),float(y),float(z)] for z in (0,1) for y in (0,1) for x in (0,1)]
        cells=[[0,1,3,7],[0,3,2,7],[0,2,6,7],[0,6,4,7],[0,4,5,7],[0,5,1,7]]
        middle={};quadratic=[]
        for cell in cells:
            row=cell.copy()
            for a,b in ref.EDGES:
                edge=tuple(sorted((cell[a],cell[b])))
                if edge not in middle:middle[edge]=len(points);points.append(((np.asarray(points[edge[0]])+points[edge[1]])/2).tolist())
                row.append(middle[edge])
            quadratic.append(row)
        mesh=dict(nodes=points,tetrahedra=cells,tetrahedra10=quadratic,quadrature_order=3,material_ids=[1]*6,
            region_ids=[1]*6,regions=[dict(id=1,body='Clay',material_id=1,bounds_m=[[0.,0.,0.],[1.,1.,1.]])],node_sets={})
        ground=dict(schema='soil-fem-ground/1',units='m-kN-kPa',materials=[dict(id=1,unit_weight_kn_m3=0.,mcc=dict(
            critical_state_m=1.2,**{'lambda':.2,'kappa':.04},poisson=.3,specific_volume_ncl_at_reference=2.7,reference_pressure_kpa=1.))],
            layers=[dict(body='Clay',material_id=1,k0_x=1.,k0_y=1.,consolidation=dict(definition='mcc_pc_over_mean_effective',ratio=1.1),void_ratio=dict(mode='from_mcc_reference'))],
            initialization=dict(method='horizontal_k0',surface_surcharge_kpa=100.,water_level_m=None,relative_equilibrium_tolerance=1e-9))
        case=dict(schema='soil-fem-projected-mcc-stages/1',settings=dict(kinematics='side_rollers'),pressure_on_original_surface_only=True,
            stages=[dict(kind='pressure',name='preload',target_pressure_kpa=140.),dict(kind='pressure',name='unload',target_pressure_kpa=110.)])
        with tempfile.TemporaryDirectory() as tmp,contextlib.redirect_stdout(io.StringIO()):
            folder=Path(tmp);files=[folder/name for name in ('mesh.json','ground.json','case.json','mesh-identity.json')]
            for path,data in zip(files,(mesh,ground,case,{}),strict=True):path.write_text(json.dumps(data),encoding='utf-8')
            outputs=[]
            for name,factory in [('fresh',fresh.SparseProjectedSystem),('reused',self.candidate_factory)]:
                start=folder/(name+'-initial.json');end=folder/(name+'-cycle.json')
                extra=[Path(fresh.__file__),Path(reused.__file__),Path(matrixfree.__file__)]
                initial.solve(*files,start,1.,40,factory,extra)
                cycle.solve(*files,start,end,40,factory,extra)
                outputs.append(json.loads(end.read_text()))
            for left,right in zip(outputs[0]['frames'],outputs[1]['frames'],strict=True):
                for key in ('displacement_m','material_state'):
                    np.testing.assert_allclose(left[key],right[key],rtol=1e-8,atol=1e-10)
                self.assertLess(right['relative_work_conjugate_force_residual'],1e-8)
                self.assertLess(right['relative_frozen_volume_constraint_residual'],1e-8)
            self.assertGreater(outputs[1]['frames'][1]['plastic_points'],0)
            stats=outputs[1]['frames'][-1]['sparse_solver_diagnostics']
            self.assertLess(stats['actual_tangent_lu_factorizations'],stats['linear_rate_systems'])
            self.assertGreater(stats['cached_factor_hits'],0)


class MatrixFreeReference(FactorReference):
    candidate_factory=matrixfree.MatrixFreeFactorProjectedSystem

    def test_matrix_free_action_matches_independently_assembled_three_field_matrix(self):
        from scipy.sparse import bmat
        import mcc_rate_fem as rate
        a,b=self.systems();state=self.state(a,100.,110.)
        _,d,_,dn,den=rate.terms(state,a.parameters)
        active=np.arange(len(state))%2==0;ep=d.copy()
        ep[active]-=np.einsum('qa,qb,q->qab',dn[active],dn[active],1/den[active])
        dd,dt,td,tt=a.tangent_blocks(ep);m,div=a.mass,a.reduced_divergence
        matrix=bmat([[tt,-td,-m],[-dt,dd,-div.T],[-m,-div,None]],format='csc')
        vector=np.random.default_rng(917).normal(size=matrix.shape[0])
        expected=matrix@vector;actual,_=b.original_action(ep,vector)
        self.assertLess(np.linalg.norm(actual-expected)/np.linalg.norm(expected),1e-12)


if __name__=='__main__':unittest.main()
