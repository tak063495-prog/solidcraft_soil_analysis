"""Light A-01 diagnostics; heavy global histories remain sealed, not rerun by CI."""
import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
STUDY=ROOT/'interim/2026-12-20-a01-joint-combo-cause-diagnosis'
sys.path.insert(0,str(STUDY))
from a01cause_metrics import dot,norm,metric
from a01cause_work import identity,replay,geometry
from a01cause_reference import Model,run
from a01cause_binding import fixtures,original_reference
from a01cause_linear import assess

class A01JointComboCauseDiagnosis(unittest.TestCase):
    def test_fixed_order_reductions_and_invalid_correspondence(self):
        self.assertEqual(dot([1e16,1.,-1e16],[1.,1.,1.]),1.)
        self.assertEqual(norm([[3.,4.],[0.,12.]]),13.)
        with self.assertRaises(ValueError):metric([1.],[1.,2.])
        with self.assertRaises(ValueError):metric([float('nan')],[0.])

    def test_endpoint_identity_keeps_unexplained_substep_work(self):
        f=dict(jump=np.array([[[0.,0.,-.1]],[[.2,0.,-.05]]]),traction=np.array([[[0.,0.,-10.]],[[2.5,0.,-5.]]]),
            reference_jump=np.zeros((2,1,3)),initial_traction=np.zeros((2,1,3)),elastic_offset=np.zeros((2,1,3)),
            joint_active=np.ones((2,1),bool),generation=np.ones((2,1),int),endpoint_work=np.array([[1.],[1.8]]))
        b=copy.deepcopy(f);b['endpoint_work'][1,0]=1.7
        r=identity([dict(progress=.5),dict(progress=1.)],f,b)
        self.assertEqual(norm(r['traction_term']),0.);self.assertEqual(norm(r['jump_term']),0.)
        self.assertGreater(norm(r['within_increment_remainder']),.09)
        self.assertEqual(r['identity_residual_max'],0.)
        s=SimpleNamespace(normal=np.array([[0.,0.,1.]]),kn=np.array([100.]),kt=np.array([100.]),mu=np.array([.5]))
        self.assertGreater(norm(replay(f,s,1)-replay(f,s,4)),.01)

    def test_reference_default_retains_original_small_input_only_solution(self):
        original=[];observed=[]
        values=fixtures.shared(steps=2)
        a=original_reference.run(original_reference.Model(*values),2,lambda r,f:original.append(f))
        b=run(Model(*values),2,lambda r,f:observed.append(f))
        self.assertEqual(a,b)
        for x,y in zip(original,observed):
            for key in x:np.testing.assert_array_equal(x[key],y[key])
        for count in (0,-1,65,1.5,True):
            with self.assertRaises(ValueError):run(None,2,None,fixed_parts=count)

    def test_duplicate_gauss_points_are_detected(self):
        system=Model(*fixtures.shared(steps=2)).system
        self.assertEqual(geometry(system)['unique_paired_faces'],2)
        bad=copy.copy(system);bad.points=system.points+[system.points[0]]
        with self.assertRaises(ValueError):geometry(bad)

    def test_same_system_svd_distinguishes_range_incompatibility(self):
        result=assess([[1.,0.],[0.,0.]],[0.,1.])
        self.assertEqual(result['numerical_rank'],1)
        self.assertEqual(result['rhs_outside_numerical_range_relative'],1.)
        self.assertEqual(result['svd_range_solve']['true_relative_residual'],1.)
        self.assertIn('error',result['direct_lu'])

    def test_document_terminal_status_with_real_crlf(self):
        text=(STUDY/'README.ja.md').read_text(encoding='utf-8')
        self.assertTrue(text.endswith('A-01は未達、passed=false。\n'))
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder)/'README.ja.md';p.write_bytes(text.replace('\n','\r\n').encode('utf-8'))
            self.assertEqual(p.read_text(encoding='utf-8'),text)

    def test_small_fixture_and_archive_parts_sha_without_resolving_equilibrium(self):
        a=json.loads((STUDY/'acceptance.json').read_text(encoding='utf-8'))
        for file,key in [('diagnostic-fixture.json','diagnostic_fixture_sha256'),('archive-index.json','archive_index_sha256')]:
            self.assertEqual(hashlib.sha256((STUDY/file).read_bytes()).hexdigest(),a[key])
        fixture=json.loads((STUDY/'diagnostic-fixture.json').read_text(encoding='utf-8'))
        for row in fixture['metric_examples']:self.assertEqual(metric(row['a'],row['b']),row['expected'])
        index=json.loads((STUDY/'archive-index.json').read_text(encoding='utf-8'))
        for part in index['parts']:
            path=STUDY/part['name'];self.assertEqual(path.stat().st_size,part['size'])
            with path.open('rb') as f:self.assertEqual(hashlib.file_digest(f,'sha256').hexdigest(),part['sha256'])
        self.assertTrue(any('linear/replay.json' in s for s in index['entries']))
        self.assertTrue(any('stacked/native/trace.jsonl' in s for s in index['entries']))

    def test_original_acceptance_refusals_and_remaining_items_are_retained(self):
        a=json.loads((STUDY/'acceptance.json').read_text(encoding='utf-8'))
        for key in ('passed','a01_complete','original_acceptance_all_met','new_limited_acceptance','pc_executed','nc01_complete','p6_complete','a02_a04_complete','b_items_complete','c01_c04_complete'):
            self.assertIs(a[key],False,key)
        self.assertTrue(a['production_defaults_unchanged']);self.assertTrue(a['refusal_scope_unchanged'])
        self.assertTrue(a['cloud_diagnostics_remaining'])

if __name__=='__main__':unittest.main()
