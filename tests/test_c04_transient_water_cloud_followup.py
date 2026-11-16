"""Lightweight frozen evidence and observation contracts; no physical reruns."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/'interim/2027-01-03-c04-transient-water-cloud-followup'
sys.path.insert(0,str(HERE))
import c04follow_common as common
import c04follow_run as run
import c04follow_report as report
import c04follow_evidence as evidence

class C04TransientWaterCloudFollowup(unittest.TestCase):
    def test_matrix_is_only_mechanics16_64_water32(self):
        plan=run.specs()
        self.assertEqual(len(plan),16)
        self.assertEqual({p['m'] for p in plan},{16,64})
        self.assertEqual({p['water'] for p in plan},{32})
        self.assertEqual(sum(p['expected_states'] for p in plan),7552)
        self.assertEqual(len({(p['m'],p['engine'],p['p'],p['r']) for p in plan}),16)

    def test_fixed_order_and_raw_denominator(self):
        self.assertEqual(report.dot([1e16,1,-1e16],[1,1,1]),1.)
        value=report.metric([3.,4.],[0.,0.])
        self.assertEqual((value['absolute_l2'],value['reference_l2'],value['denominator'],value['relative_l2']),(5.,0.,1e-12,5e12))

    def test_boundary_and_saved_zero_are_separate(self):
        zero=dict(returned_in_observation_band=True,saved_plastic_increment_exact_zero=True)
        nonzero=dict(returned_in_observation_band=True,saved_plastic_increment_exact_zero=False)
        outside=dict(returned_in_observation_band=False,saved_plastic_increment_exact_zero=True)
        self.assertEqual(report.classify(zero,zero),'returned_boundary_and_both_saved_increments_exact_zero')
        self.assertIn('nonzero',report.classify(zero,nonzero))
        self.assertIn('not_ruled_out',report.classify(outside,zero))

    def test_documents_with_actual_crlf(self):
        for name in ('README.ja.md','REPRODUCE.ja.md'):
            text=(HERE/name).read_text(encoding='utf-8')
            self.assertTrue(text.endswith('C-04は未達、passed=false。\n'))
            with tempfile.TemporaryDirectory() as temp:
                p=Path(temp)/name;p.write_bytes(text.replace('\n','\r\n').encode())
                self.assertEqual(p.read_text(encoding='utf-8'),text)

    def test_no_production_switch(self):
        for path in (ROOT/'src').rglob('*.rs'):
            self.assertNotIn('C04FOLLOW_',path.read_text(encoding='utf-8'))
            self.assertNotIn('c04follow_decision',path.read_text(encoding='utf-8'))

    def test_false_status_and_retained_refusals(self):
        for name in ('acceptance.json','measurements.json'):
            v=common.read(HERE/name)
            self.assertIs(v['passed'],False);self.assertIs(v['c04_complete'],False)
            self.assertIs(v['diagnostic_only'],True)
        v=common.read(HERE/'acceptance.json')
        self.assertIs(v['one_way_accepted'],False);self.assertIs(v['bidirectional_accepted'],False)
        self.assertIsNone(v['joint_mass_balance'])

    def test_retained_category_does_not_fabricate_current_flow(self):
        v=common.read(HERE/'category-mechanism-review.json')
        self.assertFalse(v['mechanical32_reclassified'])
        for step,data in v['levels'].items():
            self.assertTrue(data['native_current_shortcut_all'])
            self.assertEqual(data['independent_current_lhs_values'],[0.])
            self.assertGreater(data['inherited_native_last_flow_range_m'][0],0.)
            self.assertTrue(data['all_both_saved_increment_exact_zero'])
            self.assertEqual(data['status_pairs'],{'(2, 1)':data['comparison_occurrences']})

    def test_sealed_sha_not_recomputation(self):
        m=common.read(HERE/'archive-index.json');combined=hashlib.sha256();size=0
        for part in m['parts']:
            p=HERE/part['name'];self.assertEqual(p.stat().st_size,part['size'])
            self.assertEqual(evidence.sha(p),part['sha256'])
            with p.open('rb') as f:
                for raw in iter(lambda:f.read(1024*1024),b''):combined.update(raw);size+=len(raw)
        self.assertEqual((combined.hexdigest(),size),(m['archive_sha256'],m['archive_bytes']))
        self.assertIn('c04follow-formal/result.json',m['entries'])
        self.assertIn('c04follow-observations/measurements.json',m['entries'])

    def test_formal_completeness_and_full_update_counts(self):
        m=common.read(HERE/'measurements.json')
        self.assertEqual((m['fresh_equilibrium_paths'],m['fresh_accepted_states']),(16,7552))
        self.assertFalse(m['rates_calculated']);self.assertFalse(m['observation_band_is_acceptance'])
        for step,level in m['levels'].items():
            self.assertEqual(len(level['all_cross_policy_pairs']),16)
            for engine,policies in level['contact_partition_counts'].items():
                self.assertEqual(set(policies),{'NN','NI','IN','II'})
                for counts in policies.values():self.assertEqual(sum(counts.values()),6*(7*int(step)+192))
            for pair in level['all_cross_policy_pairs'].values():
                self.assertEqual(sum(pair['classifications'].values()),pair['category_mismatches'])
                q=pair['fields']['plastic_dissipation']
                self.assertEqual(q['relative_l2'],q['absolute_l2']/q['denominator'])

if __name__=='__main__':unittest.main()
