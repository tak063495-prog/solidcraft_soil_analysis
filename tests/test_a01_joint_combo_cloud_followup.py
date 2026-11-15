"""Light contracts and sealed evidence; never rerun the diagnostic FEM matrix."""
import importlib, json
from pathlib import Path
import sys, tempfile, unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/'interim/2026-12-28-a01-joint-combo-cloud-followup'
sys.path.insert(0,str(HERE))
from a01follow_contact import point
from a01follow_evidence import verify
from a01follow_binding import norm, dot

def read(name):return json.loads((HERE/name).read_text(encoding='utf-8'))

class A01CloudFollowup(unittest.TestCase):
    def test_whole_scope_stays_incomplete(self):
        v=read('acceptance.json')
        for key in ('passed','a01_complete','new_limited_acceptance','production_rust_changed','nc_code_or_configuration_changed','acceptance_tolerance_changed','old_values_rewritten','physical_instability_proved','a02_a04_b01_b03_c01_c04_nc01_p6_complete'):
            self.assertIs(v[key],False,key)
        self.assertIs(v['refusals_retained'],True)

    def test_documents_use_universal_newlines(self):
        tail='A-01は未達、passed=false。\n'
        for name in ('README.ja.md','REPRODUCE.ja.md'):
            text=(HERE/name).read_text(encoding='utf-8')
            self.assertTrue(text.endswith(tail))
            with tempfile.TemporaryDirectory() as d:
                p=Path(d)/name;p.write_bytes(text.replace('\n','\r\n').encode('utf-8'))
                self.assertTrue(p.read_text(encoding='utf-8').endswith(tail))

    def test_unavailable_primary_bodies_are_not_claimed(self):
        v=read('primary-source-review.json')
        self.assertIs(v['external_reference_values_introduced'],False)
        self.assertIs(v['book_body_claims_invented'],False)
        books=[s for s in v['sources'] if any(k in s['title'] for k in ('Wriggers','Laursen','Nocedal','Dennis'))]
        self.assertEqual(len(books),4)
        self.assertTrue(all('未取得' in s['access'] or '目次だけ' in s['access'] for s in books))

    def test_fixed_order_fsum_does_not_absorb_order_error(self):
        self.assertEqual(norm([[3.,0.],[0.,4.]]),5.)
        self.assertEqual(dot([1e16,1.,-1e16],[1.,1.,1.]),1.)

    def test_analytic_sliding_derivative_has_normal_coupling(self):
        q,d,branch=point([0,0,1],2.,1.,.5,[0,0,-.25],[0,0,-.5],[0,0,0],[.5,0,-.25],1)
        np.testing.assert_array_equal(q,[.25,0,-.5])
        np.testing.assert_array_equal(d,[[0,0,-1],[0,.5,0],[0,0,2]])
        self.assertEqual(branch['path'],['slip'])

    def test_exact_boundary_is_reported(self):
        _,_,b=point([0,0,1],2.,1.,.5,[0,0,0],[0,0,0],[0,0,0],[0,0,0],1)
        self.assertIn('gap=0',b['boundaries'])
        self.assertIn('friction=cap',b['boundaries'])

    def test_full_bidirectional_inventory_is_retained(self):
        m=read('measurements.json');rows=m['exchange_paths']
        self.assertEqual(len(rows),24)
        self.assertEqual(len({(r['case'],r['engine'],r['predictor'],r['repartition']) for r in rows}),24)
        self.assertTrue(all(r['complete'] is True and r['accepted_states']==96 for r in rows))
        self.assertIs(m['passed'],False)
        for case in m['cases'].values():
            self.assertEqual(len(case['paired_work_relative_l2']),4)
            self.assertTrue(all(v==0 for v in case['paired_open_mismatches'].values()))

    def test_two_contact_neutral_modes_are_separated(self):
        m=read('mode-separation.json')
        self.assertEqual(m['branch_gp_indices'],dict(open=[0,1,3,4,6,7,8,9,10,11],stick=[],slip=[2,5]))
        self.assertEqual(m['linearizations']['analytic']['middle_rigid_contact_rank'],4)
        self.assertEqual(m['linearizations']['native_central_fd']['middle_rigid_contact_rank'],5)
        self.assertEqual(m['analytic_traction_vs_native_maximum'],0.)
        self.assertIs(m['new_support_or_regularization'],False)

    def test_local_descent_and_finite_candidate_exhaustion_are_distinct(self):
        v=read('descent-decomposition.json');r=v['records'][-1]
        self.assertIs(v['strict_decrease_and_18_candidates_unchanged'],True)
        self.assertLess(r['source_actual_merit_slope'],0.)
        baseline=r['factorial'][0]
        self.assertEqual(len(baseline['candidate_norms']),18)
        self.assertEqual(baseline['accepted_cuts'],[])
        self.assertGreater(baseline['best_norm_change'],0.)
        extra=r['outside_candidate_grid_observations'][0]
        self.assertEqual(extra['cut'],18);self.assertLess(extra['decrease'],0.);self.assertIs(extra['adopted'],False)
        self.assertEqual(r['initial_projected_force_norm'],0.)

    def test_production_has_no_new_policy_switch(self):
        self.assertNotIn('A01FOLLOW_', (ROOT/'src/stages.rs').read_text(encoding='utf-8'))
        v=read('base-preservation.json')
        self.assertIs(v['all_materialized_old_bytes_equal_git'],True)
        self.assertEqual(v['unexpected_missing'],[]);self.assertEqual(v['changed_old_paths'],[])

    def test_import_does_not_replace_generic_legacy_modules(self):
        names=('fixtures','reference','fields','study')
        before={k:sys.modules.get(k) for k in names}
        importlib.reload(importlib.import_module('a01follow_binding'))
        self.assertEqual(before,{k:sys.modules.get(k) for k in names})

    def test_new_raw_evidence_is_fully_sha_verified(self):
        result=verify(HERE)
        self.assertIs(result['all_sha_verified'],True)
        self.assertEqual(result['entries'],len(read('archive-index.json')['entries']))

if __name__=='__main__':unittest.main()
