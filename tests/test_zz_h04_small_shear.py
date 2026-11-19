"""Check sealed inputs and remeasure one signed response without a path solve."""
import hashlib
import json
import math
from pathlib import Path
import tempfile
import unittest
import zipfile


HERE = Path(__file__).resolve().parents[1]/'interim/2027-01-16-a03-h04-small-explicit-shear'


class SmallShearAtUnequalThickness(unittest.TestCase):
    def test_record_and_fixed_matrix(self):
        acceptance = json.loads((HERE/'acceptance.json').read_bytes())
        raw = (HERE/'measurements.json').read_bytes()
        self.assertEqual(hashlib.sha256(raw).hexdigest(), acceptance['measurements_sha256'])
        measured = json.loads(raw)
        self.assertIs(acceptance['passed'], False)
        self.assertIs(acceptance['a03_complete'], False)
        self.assertIs(acceptance['remaining_same_direction_cloud_diagnostics'], False)
        paths = measured['paths']
        self.assertEqual(len(paths), 21)
        self.assertEqual(sum(len(path['rows']) for path in paths.values()), 768)
        for stress in ('zero','xz1m','xz1p'):
            gl = measured['comparisons'][stress+'/gl5-gl7/common16']
            self.assertEqual(gl['states'], 16)
            for n in (4,8,16):
                for e in ('native','independent'):
                    self.assertEqual(len(paths[f'h4-{stress}-n{n}-gl5/{e}']['rows']), 4*n)

    def test_raw_horizontal_and_vertical_displacements(self):
        index = json.loads((HERE/'archive-index.json').read_bytes())
        acceptance = json.loads((HERE/'acceptance.json').read_bytes())
        digest = hashlib.sha256()
        with tempfile.TemporaryFile() as stream:
            for part in index['parts']:
                data = (HERE/part['name']).read_bytes()
                self.assertEqual((len(data), hashlib.sha256(data).hexdigest()), (part['size'],part['sha256']))
                digest.update(data)
                stream.write(data)
            self.assertEqual(digest.hexdigest(), acceptance['archive_sha256'])
            self.assertEqual(digest.hexdigest(), index['archive_sha256'])
            with zipfile.ZipFile(stream) as archive:
                states = {}
                for s, tau in (('zero',0.),('xz1m',-.05),('xz1p',.05)):
                    prefix = f'formal/h4-{s}-n8-gl5/'
                    audit = json.loads(archive.read(prefix+'input-audit.json'))
                    self.assertEqual(audit['height_m'], .4)
                    self.assertEqual(audit['fresh']['effective_stress_compression_positive_kpa'],
                                     [2.5,2.5,5.,0.,0.,tau])
                    with archive.open(prefix+'independent-fields.jsonl') as rows:
                        records = [json.loads(line) for line in rows]
                    states[s] = [(r['row'],r['fields']['u']) for r in records if (4*r['row']['progress']).is_integer()]
        observed = json.loads((HERE/'measurements.json').read_bytes())['directional']['n8/independent']['u']
        for s, key in (('xz1p','plus_l2'),('xz1m','minus_l2')):
            self.assertEqual([r for r,_ in states[s]], [r for r,_ in states['zero']])
            self.assertEqual(len(states[s]), 16)
            for axis in range(3):
                diff = (a[axis]-b[axis] for (_,left),(_,right) in zip(states[s],states['zero'],strict=True)
                        for a,b in zip(left,right,strict=True))
                self.assertEqual(math.sqrt(math.fsum(x*x for x in diff)), observed[axis][key])
