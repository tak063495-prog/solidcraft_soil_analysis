"""Bounded C-04 observation checks; archived global controls are not rerun in CI."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
STUDY=ROOT/'interim/2026-12-19-c04-transient-water-cause-diagnosis'
sys.path.insert(0,str(STUDY))
from c04cause_metrics import norm,dot,metric,compare
from c04cause_observe import dissipation_identity,margin_bin,replay
from c04cause_reference import run


class C04TransientWaterCauseDiagnosis(unittest.TestCase):
    def test_fixed_reduction_and_invalid_correspondence(self):
        self.assertEqual(dot([1e16,1.,-1e16],[1.,1.,1.]),1.)
        self.assertEqual(norm([[3.,4.],[0.,12.]]),13.)
        with self.assertRaises(ValueError):metric(np.zeros(2),np.zeros(3))
        with self.assertRaises(ValueError):metric([float('nan')],[0.])
        a=([dict(stage='water',progress=.5)],dict(u=np.zeros((1,1,3)),traction=np.zeros((1,1,3))))
        b=([dict(stage='water',progress=1.)],a[1])
        with self.assertRaises(ValueError):compare(a,b)

    def test_endpoint_products_do_not_determine_substep_dissipation(self):
        system=SimpleNamespace(normal=np.array([[0.,0.,1.]]),kn=np.array([100.]),kt=np.array([100.]),mu=np.array([.5]))
        # The shear force is saturated while the normal pressure falls along a straight jump segment.
        jumps=np.array([[[0.,0.,-.1]],[[.2,0.,-.05]]])
        coarse=replay(jumps,system,1);fine=replay(jumps,system,4)
        np.testing.assert_allclose(coarse['traction'],fine['traction'],rtol=0,atol=1e-14)
        np.testing.assert_allclose(coarse['plastic_slip'],fine['plastic_slip'],rtol=0,atol=1e-14)
        self.assertNotEqual(coarse['plastic_dissipation'][-1,0],fine['plastic_dissipation'][-1,0])
        a=([],coarse);b=([],fine)
        result=dissipation_identity(a,b,system)
        self.assertEqual(result['identity_residual_max'],0.)
        self.assertGreater(norm(result['within_increment_remainder']),.01)

    def test_margin_bins_are_observation_labels_and_keep_boundary_points(self):
        self.assertEqual(margin_bin(0.),'0..1e-14')
        self.assertEqual(margin_bin(-1e-15),'0..1e-14')
        self.assertEqual(margin_bin(1e-7),'>1e-8')
        fixture=json.loads((STUDY/'diagnostic-fixture.json').read_text(encoding='utf-8'))
        self.assertEqual(len(fixture['finest_mismatch_details']),133)
        self.assertTrue(all(r['margin_bin']=='0..1e-14' for r in fixture['finest_mismatch_details']))
        self.assertTrue(all(r['zero_saved_plastic_increment_both'] for r in fixture['finest_mismatch_details']))
        self.assertEqual(sum(fixture['finest_margin_bins'].values()),133)

    def test_diagnostic_fixed_parts_reject_out_of_budget_inputs(self):
        for count in (0,65,-1,1.5):
            with self.assertRaises(ValueError):run(None,fixed_contact_parts=count)

    def test_sealed_small_metrics_exact_without_float_tolerance(self):
        fixture=json.loads((STUDY/'diagnostic-fixture.json').read_text(encoding='utf-8'))
        for row in fixture['metric_examples']:
            self.assertEqual(metric(row['a'],row['b']),row['expected'])
        self.assertFalse(fixture['local_replay_is_global_equilibrium'])

    def test_document_terminal_status_with_actual_crlf(self):
        text=(STUDY/'README.ja.md').read_text(encoding='utf-8')
        self.assertTrue(text.endswith('C-04は未達、passed=false。\n'))
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'README.ja.md'
            path.write_bytes(text.replace('\n','\r\n').encode('utf-8'))
            self.assertEqual(path.read_text(encoding='utf-8'),text)

    def test_parts_and_fixture_are_sha_fixed_without_equilibrium_recalculation(self):
        acceptance=json.loads((STUDY/'acceptance.json').read_text(encoding='utf-8'))
        index_bytes=(STUDY/'archive-index.json').read_bytes()
        self.assertEqual(hashlib.sha256(index_bytes).hexdigest(),acceptance['archive_index_sha256'])
        self.assertEqual(hashlib.sha256((STUDY/'diagnostic-fixture.json').read_bytes()).hexdigest(),acceptance['diagnostic_fixture_sha256'])
        index=json.loads(index_bytes)
        for part in index['parts']:
            path=STUDY/part['name'];self.assertEqual(path.stat().st_size,part['size'])
            with path.open('rb') as stream:self.assertEqual(hashlib.file_digest(stream,'sha256').hexdigest(),part['sha256'])
        self.assertTrue(any('fixed4-w32/fields.jsonl' in name for name in index['entries']))
        self.assertTrue(any('native-w32/raw/accepted.jsonl' in name for name in index['entries']))

    def test_original_acceptance_and_richards_scope_remain_unachieved(self):
        a=json.loads((STUDY/'acceptance.json').read_text(encoding='utf-8'))
        for key in ('passed','c04_complete','a04_complete','c01_complete','c02_complete','c03_complete',
                    'nc01_complete','p6_complete','ab_items_completed','pc_executed','original_acceptance_all_met',
                    'richards_wet_joint_connected'):
            self.assertIs(a[key],False,key)
        self.assertTrue(a['production_defaults_unchanged'])
        self.assertTrue(a['old_evidence_unchanged'])
        self.assertTrue(a['cloud_diagnostics_remaining'])

if __name__=='__main__':unittest.main()
