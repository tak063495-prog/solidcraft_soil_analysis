"""Bounded rotating-law and intervention checks; heavy replay remains sealed."""
import importlib.util
import json
from math import fsum
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/'interim/2026-12-16-c01-regen-global-control'
sys.path.insert(0,str(HERE))
from c01_regen_rotating import build,request,compare,KN,KT,MU
from c01_regen_schedule import project_forced,open_at


class RegenerationControlTests(unittest.TestCase):
    def test_manufactured_rotating_path_exact_flow_and_disk(self):
        # Exact x - tau/kt - plastic is constant until an open excursion.
        samples=build('arc_variable_pressure',32)
        for r in samples:
            x=np.asarray(r['jump']);p=np.asarray(r['plastic_slip']);q=np.asarray(r['traction'])
            np.testing.assert_allclose(x[:2]-q[:2]/KT-p[:2],0.,rtol=0,atol=1e-18)
            self.assertAlmostEqual(float(np.linalg.norm(q[:2])),MU*-q[2],places=10)
        self.assertAlmostEqual(samples[-1]['dissipation'],.004*(40+120)/2+.002*(120+60)/2)
        self.assertLess(abs(samples[63]['plastic_slip'][0]),1e-18)
        self.assertGreater(samples[63]['plastic_slip'][1],0.)

    def test_manufactured_open_recontact_retains_plastic_memory(self):
        samples=build('polygon_open_recontact',8)
        opening=[r for r in samples if r['segment']==5]
        self.assertGreater(opening[3]['jump'][2],0.)
        for r in opening:
            self.assertEqual(r['traction'],[0.,0.,KN*min(r['jump'][2],0.)])
            self.assertEqual(r['plastic_slip'],opening[0]['plastic_slip'])
            self.assertEqual(r['dissipation'],opening[0]['dissipation'])
        self.assertGreater(samples[-1]['dissipation'],opening[-1]['dissipation'])

    def test_forced_open_removes_penetration_force_and_does_not_create_D(self):
        s=SimpleNamespace(normal=np.array([[0.,0.,1.]]),kn=np.array([100.]),kt=np.array([100.]),mu=np.array([.5]))
        schedule=dict(intervals={'phase':[[(0.,1.)]]},endpoint_open={'phase':[.1]})
        q,dp,acc,pd,_,_=project_forced(s,np.zeros((1,3)),np.zeros((1,3)),np.array([[.2,0.,-.1]]),[2],schedule,'phase',0.,1.)
        np.testing.assert_array_equal(q,0.);np.testing.assert_array_equal(dp,0.);np.testing.assert_array_equal(pd,0.)
        self.assertTrue(open_at(schedule,'phase',1.)[0])

    def test_forced_closed_branch_is_explicitly_modified_problem(self):
        s=SimpleNamespace(normal=np.array([[0.,0.,1.]]),kn=np.array([100.]),kt=np.array([100.]),mu=np.array([.5]))
        schedule=dict(intervals={'phase':[[]]})
        q,*_=project_forced(s,np.zeros((1,3)),np.zeros((1,3)),np.array([[0.,0.,.1]]),[1],schedule,'phase',0.,1.)
        self.assertEqual(float(q[0,0,2]),10.) # tensile force: not original unilateral solution

    def test_default_off_reference_matches_frozen_previous_reference(self):
        import c01_regen_global as g
        old,prior=g.diag.dependencies();new=g.reference_module()
        mesh=old.read(ROOT/'examples/joint-water-column-tet10-mesh.json')
        ground=old.fem.normalize_ground(old.read(ROOT/'examples/ground-joint-water-mc.json'))
        case=old.configuration(old.fem.normalize_stage_case(old.read(ROOT/'examples/stages-joint-water-cycle.json')),2,5e-5)
        a=prior._reference(mesh,ground,case,2,1,7,False,True,maximum_m=5e-5)
        b=new._reference(mesh,ground,case,2,1,7,False,True,maximum_m=5e-5)
        self.assertEqual(json.dumps(prior.serializable(a),sort_keys=True),json.dumps(new.serializable(b),sort_keys=True))

    def test_native_rotating_accuracy_and_nonnegative_D(self):
        subprocess.run(['cargo','build','--locked','--release','--example','wet_joint_local_path_probe'],cwd=ROOT,check=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        exe=ROOT/'target/release/examples'/('wet_joint_local_path_probe.exe' if sys.platform=='win32' else 'wet_joint_local_path_probe')
        req,truth=request((2,32),(1,4))
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'in.json').write_text(json.dumps(req),encoding='utf-8')
            subprocess.run([str(exe),str(p/'in.json'),str(p/'out.json')],check=True)
            ordinary=json.loads((p/'out.json').read_text());result=compare(ordinary,truth)
            from copy import deepcopy
            rotation=np.array([[0.,0.,1.],[1.,0.,0.],[0.,1.,0.]])
            rotated=deepcopy(req);rotated['normal']=(rotation@np.array(req['normal'])).tolist()
            for run in rotated['runs']:
                for step in run['steps']:step['jump_m']=(rotation@np.array(step['jump_m'])).tolist()
            (p/'rotated.json').write_text(json.dumps(rotated),encoding='utf-8')
            subprocess.run([str(exe),str(p/'rotated.json'),str(p/'rotated-out.json')],check=True)
            changed=json.loads((p/'rotated-out.json').read_text())
            for left,right in zip(ordinary,changed,strict=True):
                for a,b in zip(left['samples'],right['samples'],strict=True):
                    for key in ('plastic_slip_m','traction_kpa'):
                        np.testing.assert_array_equal(b['state'][key],rotation@np.array(a['state'][key]))
                    self.assertEqual(a['state']['plastic_dissipation_kpa_m'],b['state']['plastic_dissipation_kpa_m'])

        for row in result['runs'].values():
            self.assertGreaterEqual(row['minimum_increment_D'],0.)
            self.assertLess(row['maximum_errors']['endpoint_implementation'],1e-10)
        for name in ('arc_variable_pressure','arc_reverse','polygon_open_recontact'):
            for k in (1,4):
                for field in ('slip','traction','cumulative_D','direction_radians'):
                    self.assertLess(result['runs'][f'{name}/n32/k{k}']['maximum_errors'][field],result['runs'][f'{name}/n2/k{k}']['maximum_errors'][field])

    def test_norm_aggregation_preserves_fixed_order_and_cancellation(self):
        import c01_regen_global as g
        x=np.array([[1e16,1.],[-1e16,0.]])
        self.assertEqual(g.dot(x,np.ones_like(x)),1.)
        self.assertEqual(g.dot(np.array(x,order='F'),np.ones_like(x)),1.)

    def test_sealed_summary_and_original_failure(self):
        from c01_regen_evidence import verify
        self.assertTrue(verify()['archive_verified'])
        a=json.loads((HERE/'acceptance.json').read_text(encoding='utf-8'))
        self.assertIs(a['passed'],False)
        self.assertEqual(a['c01_decision'],'未達')

    def test_docs_lf_and_real_crlf(self):
        for name in ('README.ja.md','REPRODUCE.ja.md'):
            text=(HERE/name).read_text(encoding='utf-8')
            self.assertTrue(text.endswith('C-01は未達、passed=false。\n'))
            with tempfile.TemporaryDirectory() as tmp:
                p=Path(tmp)/name;p.write_bytes(text.replace('\n','\r\n').encode('utf-8'))
                self.assertEqual(p.read_text(encoding='utf-8'),text)

if __name__=='__main__':unittest.main()
