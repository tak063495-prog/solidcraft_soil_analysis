"""Light A-02 cause observations; global recomputation is sealed outside CI."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/'interim/2026-12-21-a02-joint-birth-cause-diagnosis'
sys.path.insert(0,str(HERE))
from a02cause_metrics import dot,norm,metric,compare,common
from a02cause_binding import study,reference
from a02cause_projection import solve
from a02cause_driver import local_replay


class A02JointBirthCauseDiagnosis(unittest.TestCase):
    def test_fixed_order_norm_and_small_denominator_are_observations(self):
        self.assertEqual(dot([1e16,1.,-1e16],[1.,1.,1.]),1.)
        self.assertEqual(norm([[3.,4.],[0.,12.]]),13.)
        m=metric([1e-18],[0.]);self.assertTrue(m['denominator_floor_active'])
        self.assertEqual(m['denominator'],1e-12)
        with self.assertRaises(ValueError):metric([1.],[1.,2.])
        with self.assertRaises(ValueError):metric([float('nan')],[0.])

    def test_coordinate_matching_never_interpolates_or_drops_fields(self):
        row=lambda p:(dict(stage='birth',progress=p),dict(u=np.array([p])))
        a=[row(.5),row(1.)];b=[row(.25),row(.5),row(.75),row(1.)]
        self.assertEqual(compare(a,common(b,2),set())['fields']['u']['absolute_l2'],0.)
        with self.assertRaises(ValueError):compare(a,b,set())
        with self.assertRaises(ValueError):compare(a,[row(.25),row(1.)],set())
        with self.assertRaises(ValueError):compare(a,[(a[0][0],{}),a[1]],set())

    def test_weighted_projection_range_and_order_preserve_tip_zero_rows(self):
        mesh,ground,case=study.fixture('transfer_force',2);s=reference.Model(mesh,ground,case).system
        a=np.einsum('qag,q->gqa',s.c,np.sqrt(s.weight)).reshape(3*len(mesh['nodes']),-1)
        f=np.einsum('qag,a,q->g',s.c,[12.,5.,-100.],s.weight)
        q,r=solve(a,f,s.weight,s.normal,s.mu);self.assertTrue(r['nodal_gate']);self.assertTrue(r['cone_gate'])
        reverse,rr=solve(a,f,s.weight,s.normal,s.mu,gp_order=np.arange(len(s.weight))[::-1])
        self.assertTrue(rr['nodal_gate']);self.assertLess(metric(reverse,q)['maximum_absolute'],1e-10)
        tip=case['joint_system']['joints'][0]['crack_tip']['edges'][0][0]
        self.assertTrue(np.all(a[3*tip:3*tip+3]==0))
        f[3*tip]+=1.;_,bad=solve(a,f,s.weight,s.normal,s.mu)
        self.assertFalse(bad['nodal_gate']);self.assertEqual(bad['maximum_nodal_error_kn'],1.)

    def test_cone_rejection_is_not_an_initial_return_mapping(self):
        mesh,ground,case=study.fixture('transfer_force',2);s=reference.Model(mesh,ground,case).system
        a=np.einsum('qag,q->gqa',s.c,np.sqrt(s.weight)).reshape(3*len(mesh['nodes']),-1)
        for target in ([0.,0.,100.],[40.,0.,-100.]):
            f=np.einsum('qag,a,q->g',s.c,target,s.weight);before=f.copy()
            q,r=solve(a,f,s.weight,s.normal,s.mu)
            self.assertTrue(r['nodal_gate']);self.assertFalse(r['cone_gate']);np.testing.assert_array_equal(before,f)
            if target[2]>0:self.assertGreater(max(r['normal_traction_kpa']),0.)
            else:self.assertGreater(max(r['friction_excess_kpa']),0.)

    def test_generation_reset_in_saved_path_replay(self):
        from types import SimpleNamespace
        s=SimpleNamespace(points=[{}],normal=np.array([[0.,0.,1.]]),kn=np.array([100.]),kt=np.array([100.]),mu=np.array([.5]))
        def f(g,j):return dict(generation=np.array([g]),reference_jump=np.zeros((1,3)),initial_traction=np.zeros((1,3)),
            elastic_offset=np.zeros((1,3)),joint_active=np.array([True]),jump=np.array([[0.,0.,j]]))
        rows=[(dict(stage='birth',progress=1.),f(1,-.1)),(dict(stage='rebirth',progress=1.),f(2,0.))]
        work=local_replay(rows,s,1);self.assertGreater(work[0,0],0.);self.assertEqual(work[1,0],0.)

    def test_document_terminal_status_with_real_crlf(self):
        text=(HERE/'README.ja.md').read_text(encoding='utf-8');self.assertTrue(text.endswith('A-02は未達、passed=false。\n'))
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'readme';p.write_bytes(text.replace('\n','\r\n').encode('utf-8'))
            self.assertEqual(p.read_text(encoding='utf-8'),text)

    def test_sealed_evidence_parts_without_recomputing_equilibrium(self):
        a=json.loads((HERE/'acceptance.json').read_text(encoding='utf-8'))
        index=HERE/'archive-index.json';self.assertEqual(hashlib.sha256(index.read_bytes()).hexdigest(),a['archive_index_sha256'])
        data=json.loads(index.read_text(encoding='utf-8'))
        for part in data['parts']:
            path=HERE/part['name'];self.assertEqual(path.stat().st_size,part['size'])
            with path.open('rb') as f:self.assertEqual(hashlib.file_digest(f,'sha256').hexdigest(),part['sha256'])
        self.assertTrue(any('projection-observations/diagnosis.json' in s for s in data['entries']))
        self.assertTrue(any('n8/stress_free-cycle/native-fields.jsonl' in s for s in data['entries']))

    def test_original_acceptance_and_cloud_remaining_status(self):
        a=json.loads((HERE/'acceptance.json').read_text(encoding='utf-8'))
        for key in ('passed','a02_complete','original_acceptance_all_met','new_limited_acceptance','pc_executed',
                    'a01_complete','a03_a04_complete','b_items_complete','c01_c04_complete','nc01_complete','p6_complete'):
            self.assertIs(a[key],False,key)
        for key in ('production_defaults_unchanged','refusal_scope_unchanged','cloud_diagnostics_remaining'):
            self.assertIs(a[key],True,key)


if __name__=='__main__':unittest.main()
