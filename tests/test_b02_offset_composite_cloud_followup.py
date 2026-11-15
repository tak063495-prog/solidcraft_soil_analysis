"""Small factorial/portability/evidence contracts; no coupled solve in CI."""
import hashlib
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
HERE=ROOT/'interim/2027-01-02-b02-offset-composite-cloud-followup'
saved=sys.path[:]
try:
    sys.path.insert(0,str(HERE))
    from b02follow_common import dot, metric, quadrature, legacy
    from b02follow_run import input_case
    from b02follow_evidence import verify
finally:
    sys.path[:]=saved

class B02OffsetCompositeCloudFollowup(unittest.TestCase):
    def test_physical_input_inverse_and_bounded_thickness_series(self):
        for scale in (.5,1.,2.):
            hashes=[]
            for mode in ('default','reverse','shuffle'):
                m,g,c,a=input_case(0,scale,mode)
                self.assertEqual(len(m['tetrahedra']),12)
                self.assertEqual([v[:4] for v in m['tetrahedra10']],m['tetrahedra'])
                self.assertTrue(a['inverse_exact'])
                self.assertEqual(a['offset_m']/a['thickness_m'],.5)
                self.assertEqual(a['seed'],20270102)
                hashes.append(a['canonical_input_sha256'])
            self.assertEqual(len(set(hashes)),1)
        with self.assertRaises(ValueError):input_case(0,1.,'invalid')
        with self.assertRaises(ValueError):input_case(2,1.,'default')
        with self.assertRaises(ValueError):input_case(0,4.,'default')

    def test_rule_is_plate_only_and_restored_after_exception(self):
        import plate6_reference
        original=plate6_reference.np
        m,g,c,a=input_case(0,1.,'default')
        soil_count=len(legacy['cases'].ref.volume(m,0,independent=True)['weight'])
        for order in (4,5,6):
            with self.assertRaisesRegex(RuntimeError,'retained'):
                with quadrature(order):
                    kernel=legacy['reference'].System(m,c['beam_system']).plate_data[0][2]
                    self.assertEqual(len(kernel.neutral.quadrature),order*order)
                    self.assertEqual(len(legacy['cases'].ref.volume(m,0,independent=True)['weight']),soil_count)
                    raise RuntimeError('retained')
            self.assertIs(plate6_reference.np,original)
        with self.assertRaises(ValueError):
            with quadrature(2):pass

    def test_fixed_sum_and_raw_denominator(self):
        self.assertEqual(dot([1e16,1.,-1e16],[1.,1.,1.]),1.)
        self.assertEqual(metric([1e-18],[0.])['relative_l2'],1e-18/1e-12)
        with self.assertRaises(ValueError):metric([math.nan],[0.])

    def test_crlf_document_suffix(self):
        for name in ('README.ja.md','REPRODUCE.ja.md'):
            text=(HERE/name).read_text(encoding='utf-8')
            self.assertTrue(text.endswith('B-02は未達、passed=false。\n'))
            with tempfile.TemporaryDirectory() as tmp:
                p=Path(tmp)/name;p.write_bytes(text.replace('\n','\r\n').encode())
                self.assertEqual(p.read_text(encoding='utf-8'),text)

    def test_saved_status_and_fixed_matrix_controls(self):
        for name in ('acceptance','measurements'):
            v=json.loads((HERE/(name+'.json')).read_text(encoding='utf-8'))
            self.assertIs(v['passed'],False);self.assertIs(v['b02_complete'],False)
            self.assertIs(v['diagnostic_only'],True);self.assertIsNone(v['unsupported_combined_equilibrium'])
        v=json.loads((HERE/'measurements.json').read_text())
        self.assertEqual((v['native_histories'],v['independent_histories']),(18,54))
        self.assertEqual(len(v['fixed_order_matrix_controls']),54)
        self.assertTrue(all(all(r['fixed_matrix_exact'].values()) and r['canonical_input_same'] for r in v['fixed_order_matrix_controls']))
        self.assertFalse(v['single_cause_identified'])

    def test_old_records_byte_bound_and_production_untouched(self):
        binding=json.loads((HERE/'prior-record-binding.json').read_text())
        for path,digest in binding.items():
            self.assertEqual(hashlib.sha256((ROOT/path).read_bytes()).hexdigest(),digest,path)
        self.assertNotIn('b02follow',(ROOT/'src/plates/mindlin.rs').read_text())

    def test_seal_and_ci_allocation(self):
        # CI hashes sealed parts only; full ZIP-entry verification is a retained manual check.
        index=json.loads((HERE/'archive-index.json').read_text())
        combined=hashlib.sha256()
        for part in index['parts']:
            path=HERE/part['name']
            self.assertEqual(path.stat().st_size,part['size'])
            with path.open('rb') as stream:
                self.assertEqual(hashlib.file_digest(stream,'sha256').hexdigest(),part['sha256'])
            with path.open('rb') as stream:
                for block in iter(lambda:stream.read(1024*1024),b''):combined.update(block)
        self.assertEqual(combined.hexdigest(),index['archive_sha256'])
        workflow=(ROOT/'.github/workflows/fem-ci.yml').read_text()
        self.assertNotIn('!/interim/2027-01-02-b02-offset-composite-cloud-followup/',workflow)
        index=json.loads((HERE/'archive-index.json').read_text())
        self.assertIn('b02follow-formal/plan.json',index['entries'])
        self.assertIn('b02follow-formal/signed-factorial-effects.json',index['entries'])

if __name__=='__main__':unittest.main()
