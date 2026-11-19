"""Read sealed binary64 components and repeat small fixed-order Decimal comparisons."""
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path
import unittest
import zipfile

HERE = Path(__file__).resolve().parents[1]/'interim/2027-01-17-a04-iteration4-tangent-rhs-correction'


class FourthDryCorrection(unittest.TestCase):
    def test_seal_and_scope(self):
        status = json.loads((HERE/'acceptance.json').read_bytes())
        raw = (HERE/'measurements.json').read_bytes()
        self.assertEqual(hashlib.sha256(raw).hexdigest(),status['measurements_sha256'])
        self.assertIs(status['passed'],False)
        self.assertIs(status['a04_complete'],False)
        self.assertIs(status['remaining_same_direction_cloud_diagnostics'],False)
        index = json.loads((HERE/'archive-index.json').read_bytes())
        part, = index['parts']
        data = (HERE/part['name']).read_bytes()
        self.assertEqual(len(data),part['size'])
        self.assertEqual(hashlib.sha256(data).hexdigest(),status['archive_sha256'])
        self.assertEqual(part['sha256'],index['archive_sha256'])
        self.assertEqual(part['sha256'],status['archive_sha256'])

    def test_measured_tangent_component_and_endpoints(self):
        result = json.loads((HERE/'measurements.json').read_bytes())['orders']
        with zipfile.ZipFile(HERE/'evidence.zip.part000') as z, localcontext() as ctx:
            ctx.prec = 40
            endpoints = []
            for order in ('000','010'):
                history = [json.loads(line) for line in z.read(f'formal/{order}/native.jsonl').splitlines()]
                matrix = next(v for v in history if v['kind']=='linear')
                exact = json.loads(z.read(f'retained/a04it4-final40/{order}.json'))
                i,j = result[order]['native']['tangent_vs_elastic']['index']
                gi,gj = (int(exact['ids'][n]) for n in (i,j))
                (ri,ci), = matrix['rows'][gi]
                (rj,cj), = matrix['rows'][gj]
                self.assertEqual((ci,cj),(1.,1.))
                diff = +Decimal.from_float(matrix['matrix'][ri][rj])-Decimal(exact['K'][i][j])
                self.assertEqual(diff,Decimal(result[order]['native']['tangent_vs_elastic']['signed']))
                last = [v for v in history if v['kind']=='iterate'][-1]
                self.assertEqual(last['i'],result[order]['native_stop'])
                _,_,terminal = json.loads(z.read(f'formal/{order}/independent-path.json'))
                self.assertEqual(terminal['i'],2)
                endpoints.append([Decimal.from_float(f)-Decimal.from_float(p)
                                  for f,p in zip(last['internal'],last['target'],strict=True)])
            measured = max(abs(b-a) for a,b in zip(*endpoints,strict=True))
            saved = result['010_minus_000']['observed_endpoint_reaction_kn']
            self.assertEqual(float(measured),saved['max_abs'])
