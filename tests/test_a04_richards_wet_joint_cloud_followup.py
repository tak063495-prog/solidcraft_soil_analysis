"""Light contracts and sealed raw evidence only; no equilibrium recomputation."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/'interim/2026-12-31-a04-richards-wet-joint-cloud-followup'
sys.path.insert(0,str(HERE))
import a04follow_common as common
import a04follow_run as run
import a04follow_report as report
import a04follow_evidence as evidence

class A04RichardsWetJointCloudFollowup(unittest.TestCase):
    def test_bounded_matrix_and_false_status(self):
        p=run.specs()
        self.assertEqual((len(p['soil_cross']),len(p['order']),len(p['exchange'])),(36,48,8))
        self.assertFalse(p['a04_complete'])
        self.assertLessEqual(max(v['layers']*6 for v in p['soil_cross']),48)
        self.assertTrue(p['diagnostic_only'])

    def test_fsum_fixed_order(self):
        self.assertEqual(report.dot([1e16,1,-1e16],[1,1,1]),1.)
        m=common.compare([3.,4.],[0.,0.])
        self.assertEqual(m['absolute_l2'],5.)
        self.assertEqual(m['reference_l2'],0.)
        self.assertEqual(m['relative_l2'],5e12)

    def test_rate_without_extrapolation(self):
        r=report.rates([[1.],[.5],[.25]])
        self.assertEqual(r['observed_log2_rate'],1.)
        self.assertFalse(r['extrapolation'])
        self.assertIsNone(report.rates([[0.],[0.],[0.]])['observed_log2_rate'])

    def test_physical_contact_canonicalization(self):
        layout=dict(element_keys=[[4,5,6,7],[0,1,2,3]],contact_keys=[[[4],[5],[.8,.2]],[[0],[1],[.2,.8]]])
        fields=dict(u=[[1.,2.,3.]],traction=[[5.,0.,0.],[2.,0.,0.]],stress=[[5.]*6,[2.]*6])
        out=report.canonical_joint(fields,layout)
        self.assertEqual(out['traction'][:,0].tolist(),[2.,5.])
        self.assertEqual(out['stress'][:,0].tolist(),[2.,5.])
        self.assertEqual(out['u'].tolist(),fields['u'])

    def test_production_default_is_untouched(self):
        for p in (ROOT/'src').rglob('*.rs'):
            self.assertNotIn('A04FOLLOW_',p.read_text(encoding='utf-8'))

    def test_document_end_with_real_crlf(self):
        for name in ('README.ja.md','REPRODUCE.ja.md'):
            text=(HERE/name).read_text(encoding='utf-8')
            self.assertTrue(text.endswith('A-04は未達、passed=false。\n'))
            with tempfile.TemporaryDirectory() as tmp:
                p=Path(tmp)/name;p.write_bytes(text.replace('\n','\r\n').encode('utf-8'))
                self.assertEqual(p.read_text(encoding='utf-8'),text)

    def test_scientific_status_is_not_software_status(self):
        for name in ('acceptance.json','measurements.json'):
            v=common.read(HERE/name)
            self.assertIs(v['passed'],False)
            self.assertIs(v['a04_complete'],False)
            self.assertIs(v['diagnostic_only'],True)
        v=common.read(HERE/'acceptance.json')
        self.assertIsNone(v['joint_mass_balance'])
        self.assertIs(v['one_way_accepted'],False)
        self.assertIs(v['bidirectional_accepted'],False)

    def test_order_failure_not_relaxed(self):
        accepted=common.read(HERE/'acceptance.json')
        self.assertIs(accepted['ordering_nodal_force_invariance_passed'],False)
        measurements=common.read(HERE/'measurements.json')
        pairs=measurements['ordering']['dry-native']['effects']
        maximum=max(v['fields']['reaction']['maximum_absolute'] for p in pairs.values() for v in p)
        self.assertGreater(maximum,1e-10)
        stop=common.read(HERE/'order-newton-stop-observation.json')
        self.assertEqual(stop['original_force_tolerance'],1e-11)
        self.assertEqual(stop['records']['e0f1q0']['attempt']['newton_iterations'],5)

    def test_sealed_parts_and_index(self):
        index=common.read(HERE/'archive-index.json')
        combined=hashlib.sha256();size=0
        for part in index['parts']:
            digest=hashlib.sha256();count=0
            with (HERE/part['name']).open('rb') as f:
                for raw in iter(lambda:f.read(1024*1024),b''):
                    digest.update(raw);combined.update(raw);count+=len(raw)
            self.assertEqual((count,digest.hexdigest()),(part['size'],part['sha256']))
            size+=count
        self.assertEqual((size,combined.hexdigest()),(index['archive_bytes'],index['archive_sha256']))
        self.assertIn('formal/solve-ledger.json',index['entries'])
        self.assertGreater(len(index['entries']),0)

if __name__=='__main__':unittest.main()
