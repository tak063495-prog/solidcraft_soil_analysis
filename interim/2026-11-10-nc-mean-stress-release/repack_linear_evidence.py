"""Apply bitwise XOR compression; verify identical metadata and float64 bits."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools'))
from linear_evidence_codec import read_rows, encode_rows, ARRAYS


def signature(path):
    digest=hashlib.sha256();count=0
    for row in read_rows(path):
        header=dict(keys=list(row),metadata={k:v for k,v in row.items() if k not in ARRAYS},
                    shapes={k:np.asarray(row[k]).shape for k in ARRAYS})
        digest.update((json.dumps(header,allow_nan=False)+'\n').encode())
        for key in ARRAYS:digest.update(np.asarray(row[key],dtype='<f8').tobytes())
        count+=1
    return digest.hexdigest(),count


def run(evidence):
    manifest=evidence/'compressed-jsonl.json';records=json.loads(manifest.read_text())
    for name,record in records.items():
        if record.get('encoding')!='float64-lossless' or record.get('bitwise_delta'):
            continue
        source=evidence/(name+'.f64.gz');target=source.with_name(source.name+'.new')
        before,count=signature(source);encode_rows(read_rows(source),target);after,other_count=signature(target)
        if before!=after or count!=other_count:raise ValueError('bitwise compression changed actual linear evidence')
        record.update(previous_encoded_sha256=record['encoded_sha256'],encoded_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
                      encoded_size=target.stat().st_size,bitwise_delta=True,metadata_and_float64_bits_sha256=after,actual_systems=count)
        target.replace(source);manifest.write_text(json.dumps(records,indent=2)+'\n')
        print(name,record['encoded_size'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--evidence',type=Path,required=True)
    run(p.parse_args().evidence.resolve())
