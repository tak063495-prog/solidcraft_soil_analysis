"""Check both uncondensed MCC/member rate forms against a full elastic solve."""
import os
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[name]='1'
from pathlib import Path
import sys
import unittest
import contextlib
import io
import json
import tempfile
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from test_mcc_factor_reference import specimen
from beam_reference import System as Members
from mcc_structural_reference import StructuralSystem
import mcc_rate_fem as rate


def configuration(mesh):
    nodes=[dict(id=f'n{i}',position_m=mesh['nodes'][i],soil_node_set=f'n{i}') for i in range(3)]
    nodes.append(dict(id='free',position_m=(np.asarray(mesh['nodes'][1])+[.3,.4,.8]).tolist(),soil_node_set=None))
    return dict(schema='soil-fem-structures/1',nodes=nodes,
        sections=[dict(id='b',young_kpa=2e6,poisson=.2,area_m2=.03,iy_m4=1e-4,iz_m4=2e-4,
            torsion_j_m4=1e-4,unit_weight_kn_m3=24.,formulation=dict(kind='euler_bernoulli'))],
        members=[dict(id='beam',nodes=['n1','n2'],section_id='b',local_y=[0.,0.,1.],initially_active=True,
            rotational_releases=[[False,False,True],[False,False,False]],boundary_line_load_kn_m=[.1,-.2,-.3])],
        plate_sections=[dict(id='p',young_kpa=3e6,poisson=.2,thickness_m=.08,unit_weight_kn_m3=24.,drilling_factor=1e-4)],
        plates=[dict(id='plate',nodes=['n0','n1','free'],section_id='p',initially_active=True,boundary_traction_kpa=[.2,-.1,-.4])],
        rotation_constraints=[dict(node_id='n0',direction=np.eye(3)[a].tolist(),value=0.) for a in range(3)],
        translation_constraints=[],nodal_loads=[],rotation_length_scale_m=1.,maximum_condensed_dofs=100)


class StructuralRate(unittest.TestCase):
    def test_rk4_member_checkpoint_restores_exact_history_and_refuses_other_integrator(self):
        from prepare_mcc_structures import inputs
        from run_mcc_structural_reference import solve
        mesh,ground,case=inputs(Path(__file__).resolve().parents[1]);case['stages']=case['stages'][:3]
        case['stages'][2]['target_factor']=.1
        with tempfile.TemporaryDirectory() as tmp,contextlib.redirect_stdout(io.StringIO()):
            folder=Path(tmp);paths=[folder/(n+'.json') for n in ('mesh','ground','case')]
            for path,value in zip(paths,(mesh,ground,case)):path.write_text(json.dumps(value))
            solve(*paths,folder/'whole.json',16,integration_method='rk4',checkpoint_dir=folder/'checkpoints')
            checkpoint=folder/'checkpoints/phase-002.json'
            solve(*paths,folder/'resumed.json',16,integration_method='rk4',resume_file=checkpoint)
            whole=json.loads((folder/'whole.json').read_text());resumed=json.loads((folder/'resumed.json').read_text())
            self.assertEqual(whole['integration_method'],'rk4')
            self.assertEqual(whole['release_audits'],resumed['release_audits'])
            for a,b in zip(whole['frames'],resumed['frames'],strict=True):
                for key in ('material_state','displacement_m','beams','active_elements','primary_pressure'):
                    self.assertEqual(a[key],b[key],key)
            self.assertTrue(any(k.endswith('mcc_rate_rk4.py') for k in whole['input_and_source_sha256']))
            with self.assertRaisesRegex(ValueError,'binding changed'):
                solve(*paths,folder/'wrong-method.json',16,integration_method='midpoint',resume_file=checkpoint)
            self.assertFalse((folder/'wrong-method.json').exists())

    def test_independent_stage_checkpoints_restore_member_and_material_histories_and_reject_changed_bindings(self):
        from prepare_mcc_structures import inputs
        from run_mcc_structural_reference import solve
        mesh,ground,case=inputs(Path(__file__).resolve().parents[1])
        original=case['stages']
        case['stages']=[original[0],original[1],original[2],next(s for s in original if s['kind']=='excavate'),original[-1]]
        case['stages'][2]['target_factor']=.1
        with tempfile.TemporaryDirectory() as tmp,contextlib.redirect_stdout(io.StringIO()):
            folder=Path(tmp);paths=[folder/(n+'.json') for n in ('mesh','ground','case')]
            for path,value in zip(paths,(mesh,ground,case)):path.write_text(json.dumps(value))
            solve(*paths,folder/'continuous.json',16,checkpoint_dir=folder/'checkpoints')
            continuous=json.loads((folder/'continuous.json').read_text())
            self.assertGreater(continuous['frames'][1]['plastic_points'],0)
            self.assertFalse(all(continuous['frames'][-1]['active_elements']))
            for phase in (2,4):
                output=folder/f'resumed-{phase}.json'
                solve(*paths,output,16,resume_file=folder/'checkpoints'/f'phase-{phase:03}.json')
                resumed=json.loads(output.read_text())
                self.assertEqual(resumed['release_audits'],continuous['release_audits'])
                for a,b in zip(continuous['frames'],resumed['frames'],strict=True):
                    for name in ('material_state','displacement_m','beams','active_elements','primary_pressure','plastic_points'):
                        self.assertEqual(a[name],b[name],(phase,name))
            checkpoint=folder/'checkpoints/phase-001.json';saved=json.loads(checkpoint.read_text())
            with self.assertRaisesRegex(ValueError,'binding changed'):
                solve(*paths,folder/'wrong-steps.json',8,resume_file=checkpoint)
            saved['payload']['surface_pressure_kpa']+=1.;bad=folder/'tampered.json';bad.write_text(json.dumps(saved))
            with self.assertRaisesRegex(ValueError,'binding changed'):
                solve(*paths,folder/'tampered-result.json',16,resume_file=bad)
            changed=json.loads(paths[2].read_text());changed['stages'][2]['target_factor']=.2;paths[2].write_text(json.dumps(changed))
            with self.assertRaisesRegex(ValueError,'binding changed'):
                solve(*paths,folder/'changed-result.json',16,resume_file=checkpoint)

    def check_form(self,projected):
        mesh,ground,case=specimen();system=StructuralSystem(mesh,ground,[True],case,projected=projected)
        members=Members(mesh,configuration(mesh));history=members.initial_state();system.configure_structures(members,history,.7)
        state=np.zeros((system.weight.size,14));state[:,:3]=100.;state[:,6:8]=np.log([200.,2.])
        elastic=rate.terms(state,system.parameters)[1]
        force=np.random.default_rng(60261026).normal(size=system.xyz.size)
        du,de=system.linear_rate(elastic,force)
        identity=np.eye(system.xyz.size)
        transform=np.stack([system.lift(identity[i]) for i in system.free],axis=1)
        ksoil=np.zeros((system.xyz.size,system.xyz.size))
        for e,dofs in enumerate(system.dofs):
            local=np.einsum('qai,qab,qbj,q->ij',system.b[e],elastic.reshape(*system.weight.shape,6,6)[e],system.b[e],system.weight[e])
            ksoil[np.ix_(dofs,dofs)]+=local
        kb,_,f0,*_=members.assemble(history,0.);_,_,f1,*_=members.assemble(history,1.)
        _,binding,q,_=members.space(history,np.zeros_like(system.xyz),transform)
        matrix=np.block([[transform.T@ksoil@transform+binding.T@kb@binding,binding.T@kb@q],
                         [q.T@kb@binding,q.T@kb@q]])
        load=(f1-f0)*.7
        rhs=np.r_[transform.T@force+binding.T@load,q.T@load]
        direct=np.linalg.solve(matrix,rhs);expected=transform@direct[:len(system.free)]
        np.testing.assert_allclose(du,expected,rtol=1e-9,atol=1e-11)
        strain=-np.einsum('eqai,ei->eqa',system.b,expected[system.dofs]).reshape(-1,6)
        np.testing.assert_allclose(de,strain,rtol=1e-9,atol=1e-11)
        self.assertLess(max(system.maximum_block_residuals.values()),1e-10)
        self.assertGreater(system.retained_structural_dofs,0)

    def test_projected_rates_match_full_soil_plate_frame_released_rotation_system(self):self.check_form(True)
    def test_gp_volume_rates_match_full_soil_plate_frame_released_rotation_system(self):self.check_form(False)


if __name__=='__main__':unittest.main()
