"""Verify actual-solve byte reconstruction, including asymmetric systems."""
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from linear_evidence_codec import encode, read_rows, json_lines


class LinearEvidenceCodec(unittest.TestCase):
    def test_roundtrip_reconstructs_every_original_json_byte(self):
        rows = [dict(label='step-7', iteration=3, matrix=[[1., -0.], [float.fromhex('0x1.0000000000001p-34'), 3.]],
                     rhs=[np.finfo(float).tiny, 1e100], solution=[-1e-100, .3], true_relative_residual=1e-16)]
        rows.append(dict(label='step-8', iteration=0, matrix=[[1.00000000001, 0.], [-3., 3.00000000001]],
                         rhs=[-np.finfo(float).tiny, 1e100], solution=[-1e-100, -.3], true_relative_residual=2e-16))
        original = ''.join(json.dumps(x, allow_nan=False)+'\n' for x in rows).encode()
        with tempfile.TemporaryDirectory() as temp:
            a=Path(temp)/'a.gz'; b=Path(temp)/'b.gz'
            a.write_bytes(gzip.compress(original)); encode(a,b)
            self.assertEqual(b''.join(json_lines(b)),original)
            self.assertEqual(list(read_rows(b)),rows)
            # The original non-differential float64 framing remains readable.
            from linear_evidence_codec import MAGIC, ARRAYS
            import struct
            header=dict(keys=list(rows[0]),metadata={k:v for k,v in rows[0].items() if k not in ARRAYS},
                        shapes={k:np.array(rows[0][k]).shape for k in ARRAYS})
            h=json.dumps(header).encode(); raw=MAGIC+struct.pack('<I',len(h))+h
            raw+=b''.join(np.array(rows[0][k],dtype='<f8').tobytes() for k in ARRAYS)
            old=Path(temp)/'old.gz';old.write_bytes(gzip.compress(raw))
            self.assertEqual(list(read_rows(old)),rows[:1])

    def test_truncated_payload_is_rejected(self):
        from linear_evidence_codec import MAGIC
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'bad.gz'; p.write_bytes(gzip.compress(MAGIC+b'\x01\x00'))
            with self.assertRaisesRegex(ValueError,'truncated'):
                list(read_rows(p))

    def test_nonfinite_values_cannot_enter_encoded_evidence(self):
        with tempfile.TemporaryDirectory() as temp:
            a=Path(temp)/'a.gz'; b=Path(temp)/'b.gz'
            a.write_bytes(gzip.compress((json.dumps(dict(matrix=[[float('nan')]],rhs=[0.],solution=[0.]))+'\n').encode()))
            with self.assertRaisesRegex(ValueError,'nonfinite'):
                encode(a,b)


if __name__ == '__main__':
    unittest.main()
