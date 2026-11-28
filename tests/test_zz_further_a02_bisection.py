"""Verify the archived 24-Tet10 observations without running a nonlinear solve."""
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import tempfile
import unittest
import zipfile


HERE = Path(__file__).resolve().parents[1]/'interim/2027-01-15-a02-open-interval-further-bisection'


class A02FurtherBisection(unittest.TestCase):
    def test_first_event_intervals(self):
        measured = json.loads((HERE/'measurements.json').read_bytes())
        data = (HERE/'first-events.csv').read_bytes()
        self.assertEqual(hashlib.sha256(data).hexdigest(), measured['observed_sha256']['first-events.csv'])
        rows = list(csv.DictReader(io.StringIO(data.decode('utf-8'))))
        self.assertEqual(len(rows), 4*2*32)
        for level in range(4):
            equivalents = []
            for engine in ('native', 'independent'):
                events = [r for r in rows if r['engine'] == engine and int(r['level']) == level]
                self.assertEqual([int(row['gp']) for row in events], list(range(32)))
                equivalents.append([{k: v for k, v in row.items() if k != 'engine'} for row in events])
                for row in events:
                    for key in ('open', 'slip_endpoint', 'slip_activity'):
                        lo, hi = row[key+'_lower'], row[key+'_upper']
                        if not hi:
                            self.assertEqual((key, lo), ('slip_endpoint', ''))
                            continue
                        self.assertEqual(float(hi)-float(lo), 1/(64*2**level))
                        self.assertLessEqual(1/32, float(lo))
                        self.assertLessEqual(float(hi), 1/16)
            self.assertEqual(*equivalents)

    def test_seal_and_one_interval(self):
        manifest = json.loads((HERE/'archive-index.json').read_bytes())
        acceptance = json.loads((HERE/'acceptance.json').read_bytes())
        self.assertIs(acceptance['passed'], False)
        with tempfile.TemporaryFile() as stream:
            digest = hashlib.sha256()
            for part in manifest['parts']:
                data = (HERE/part['name']).read_bytes()
                self.assertEqual(len(data), part['size'])
                self.assertEqual(hashlib.sha256(data).hexdigest(), part['sha256'])
                digest.update(data)
                stream.write(data)
            self.assertEqual(digest.hexdigest(), manifest['archive_sha256'])
            self.assertEqual(digest.hexdigest(), acceptance['archive_sha256'])
            with zipfile.ZipFile(stream) as archive:
                key = 'report/accepted-contact.jsonl'
                data = archive.read(key)
                self.assertEqual(hashlib.sha256(data).hexdigest(), manifest['entries'][key]['sha256'])
        raw = [json.loads(line) for line in data.splitlines()]
        states = {(r['path'], r['upper']): r for r in raw}
        lo, hi = 20/512, 21/512
        a, b = [states[f'L3/{engine}', hi] for engine in ('native', 'independent')]
        a0, b0 = [states[f'L3/{engine}', lo] for engine in ('native', 'independent')]
        vectors = {
            'traction_kpa': [x-y for ga, gb in zip(a['traction'], b['traction']) for x, y in zip(ga, gb)],
            'increment_work_kpa_m': [(x-x0)-(y-y0) for x, x0, y, y0 in zip(
                a['endpoint_work'], a0['endpoint_work'], b['endpoint_work'], b0['endpoint_work'])],
            'cumulative_work_kpa_m': [x-y for x, y in zip(a['endpoint_work'], b['endpoint_work'])],
        }
        table = (HERE/'interval-differences.csv').read_bytes()
        measured = json.loads((HERE/'measurements.json').read_bytes())
        self.assertEqual(hashlib.sha256(table).hexdigest(), measured['observed_sha256']['interval-differences.csv'])
        candidates = [r for r in csv.DictReader(io.StringIO(table.decode('utf-8')))
                      if r['comparison'] == 'L3 native-independent' and float(r['upper']) == hi]
        self.assertEqual(len(candidates), 1)
        for key, vec in vectors.items():
            self.assertEqual(math.sqrt(math.fsum(x*x for x in vec)), float(candidates[0][key+'_l2']))
