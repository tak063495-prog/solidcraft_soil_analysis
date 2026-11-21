"""Small contracts and sealed-byte audits only; never rerun the FEM matrix."""
import importlib,json,sys,tempfile,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/'interim/2026-12-29-a02-joint-birth-cloud-followup'
sys.path.insert(0,str(HERE))
from a02follow_binding import norm,dot
from a02follow_cone import ulp,orders,ordered_sum
from a02follow_report import switches,switchdiff
from a02follow_evidence import verify

def read(n):return json.loads((HERE/n).read_text(encoding='utf-8'))
class A02CloudFollowup(unittest.TestCase):
    def test_scope_is_incomplete_and_refusals_retained(self):
        v=read('acceptance.json')
        for k in ('passed','a02_complete','new_limited_acceptance','production_rust_changed','nc_code_or_configuration_changed','acceptance_tolerance_changed','old_values_rewritten','physical_instability_proved','other_items_complete'):
            self.assertIs(v[k],False,k)
        self.assertIs(v['refusals_retained'],True)
    def test_crlf_document_tail(self):
        tail='A-02は未達、passed=false。\n'
        for n in ('README.ja.md','REPRODUCE.ja.md'):
            s=(HERE/n).read_text(encoding='utf-8');self.assertTrue(s.endswith(tail))
            with tempfile.TemporaryDirectory() as d:
                p=Path(d)/n;p.write_bytes(s.replace('\n','\r\n').encode('utf-8'))
                self.assertTrue(p.read_text(encoding='utf-8').endswith(tail))
    def test_fixed_order_aggregation(self):
        self.assertEqual(norm([[3.,0.],[0.,4.]]),5.)
        self.assertEqual(dot([1e16,1.,-1e16],[1.,1.,1.]),1.)
    def test_ulp_input_and_order_are_separate(self):
        self.assertEqual(ordered_sum([1e16,1.,-1e16]),0.)
        self.assertEqual(ordered_sum([1e16,-1e16,1.]),1.)
        x=np.array([30.,0.,-100.]);np.testing.assert_array_equal(ulp(x,0),x)
        self.assertLess(ulp(x[0],-1),x[0]);self.assertGreater(ulp(x[0],1),x[0])
        self.assertTrue(orders(x,[0.,0.,1.],.3)['fixed_fsum']['cone'])
        self.assertFalse(orders([float(ulp(30.,1)),0.,-100.],[0.,0.,1.],.3)['fixed_fsum']['cone'])
    def test_generation_reset_is_not_a_physical_switch(self):
        def row(progress,status,generation):
            return ({'stage':'open','progress':progress},{'joint_active':np.array([True]),'generation':np.array([generation]),'contact_status':np.array([status]),'accumulated':np.array([0.]),'contact_open':np.array([status==0])})
        rows=[row(.25,1,1),row(.5,0,1),row(.75,1,2)]
        s=switches(rows);self.assertEqual(len(s),2);self.assertTrue(all(r['upper']==.5 for r in s))
        changed=[dict(r,upper=.75) for r in s]
        d=switchdiff(s,changed);self.assertEqual(d['matched'],2);self.assertEqual(d['maximum_position_difference'],.25)
        self.assertEqual(d['only_left']+d['only_right'],0)
        same_branch=[row(.25,0,1),row(.5,0,1)]
        same_branch[1][1]['accumulated'][0]=.125
        activity=switches(same_branch);self.assertEqual(len(activity),1)
        self.assertEqual(activity[0]['kind'],'slip_activity');self.assertEqual(activity[0]['endpoint_status'],0)
        self.assertIs(activity[0]['exact_internal_switch_root_known'],False)
    def test_complete_cross_and_bidirectional_exchange(self):
        m=read('measurements.json');rows=m['paths'];self.assertEqual(len(rows),78)
        cross=[r for r in rows if r['path'].startswith('cross/')];ex=[r for r in rows if r['path'].startswith('exchange/')]
        self.assertEqual(len({(r['n'],r['minimum_parts'],r['case'],r['engine']) for r in cross}),54)
        self.assertEqual(len({(r['case'],r['engine'],r['predictor'],r['repartition']) for r in ex}),24)
        self.assertTrue(all(r['complete'] and r['exit_code']==0 for r in rows));self.assertIs(m['finest_grid_is_exact_reference'],False)
        addon=m['additional_actual_contact_exchange'];self.assertEqual(len(addon['paths']),8)
        self.assertTrue(all(r['accepted_contact_parts_maximum']>1 for r in addon['path_observations'].values()))
        self.assertEqual(m['total_fresh_global_paths_including_compatibility'],129)
    def test_boundary_builds_and_linear_paths_are_recorded(self):
        v=read('cone-summary.json');self.assertEqual(v['local_requests'],36);self.assertEqual(v['projection_requests'],44)
        self.assertIs(v['unchanged_formula_tolerance_clip'],True);self.assertEqual(v['build_profiles'],['release','debug'])
        self.assertTrue(v['rejected_work_is_null']);self.assertFalse(v['windows_reproduction_confirmed'])
    def test_production_defaults_and_import_namespace(self):
        for name in ('stages.rs','joint_law.rs','joints_birth.rs'):
            self.assertNotIn('A02FOLLOW_', (ROOT/'src'/name).read_text(encoding='utf-8'))
        names=('reference','fields','study','projection','geometry_audit');before={k:sys.modules.get(k) for k in names}
        importlib.reload(importlib.import_module('a02follow_binding'));self.assertEqual(before,{k:sys.modules.get(k) for k in names})
    def test_evidence_sha(self):
        v=verify(HERE);self.assertTrue(v['all_sha256_verified']);self.assertEqual(v['entries'],len(read('archive-index.json')['entries']))
    def test_source_access_limits(self):
        v=read('primary-source-review.json');self.assertFalse(v['external_reference_values_introduced']);self.assertFalse(v['book_body_claims_invented'])
        for title in ('Wriggers','Laursen'):
            self.assertTrue(any(title in s['title'] and '本文未取得' in s['access'] for s in v['sources']))
if __name__=='__main__':unittest.main()
