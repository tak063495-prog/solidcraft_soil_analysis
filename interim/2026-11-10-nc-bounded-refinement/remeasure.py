"""Verify every evidence byte and remeasure with the sealed source snapshot."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tarfile

HERE = Path(__file__).resolve().parent


class PartsReader:
    def __init__(self, paths):
        self.paths=iter(paths); self.current=None
    def read(self, size):
        data=bytearray()
        while len(data)<size:
            if self.current is None:
                try:self.current=next(self.paths).open('rb')
                except StopIteration:break
            chunk=self.current.read(size-len(data))
            if chunk:data.extend(chunk)
            else:self.current.close();self.current=None
        return bytes(data)
    def close(self):
        if self.current is not None:self.current.close()


def run(output):
    output.mkdir(parents=True,exist_ok=False)
    manifest=json.loads((HERE/'archive-sha256.json').read_text()); total=hashlib.sha256(); size=0
    for part in manifest['parts']:
        path=HERE/part['name']; checksum=hashlib.sha256(); count=0
        with path.open('rb') as stream:
            for raw in iter(lambda:stream.read(1024*1024),b''):
                checksum.update(raw);total.update(raw);count+=len(raw);size+=len(raw)
        if checksum.hexdigest()!=part['sha256'] or count!=part['size']:
            raise ValueError('archive part changed: '+part['name'])
    if total.hexdigest()!=manifest['archive_sha256'] or size!=manifest['archive_size']:
        raise ValueError('complete archive changed')
    evidence=output/'evidence'; seen=set(); reader=PartsReader([HERE/x['name'] for x in manifest['parts']])
    try:
        with tarfile.open(fileobj=reader,mode='r|gz') as archive:
            for entry in archive:
                path=Path(entry.name)
                if not entry.isfile() or path.is_absolute() or '..' in path.parts or entry.name in seen or entry.name not in manifest['entries']:
                    raise ValueError('unsafe/unlisted archive entry')
                expected=manifest['entries'][entry.name]
                if entry.size!=expected['size']:raise ValueError('entry size changed')
                destination=evidence/path;destination.parent.mkdir(parents=True,exist_ok=True)
                checksum=hashlib.sha256();count=0
                with archive.extractfile(entry) as source,destination.open('xb') as target:
                    for raw in iter(lambda:source.read(1024*1024),b''):
                        checksum.update(raw);count+=len(raw);target.write(raw)
                if checksum.hexdigest()!=expected['sha256'] or count!=expected['size']:
                    raise ValueError('archive entry changed: '+entry.name)
                seen.add(entry.name)
    finally:
        reader.close()
    if seen!=set(manifest['entries']):raise ValueError('archive inventory changed')
    snapshot=evidence/'sources'
    script=snapshot/'interim/2026-11-10-nc-bounded-refinement/measure.py'
    subprocess.run([sys.executable,'-I',str(script),'--evidence',str(evidence),'--output',str(output/'measurements.json')],cwd=snapshot,check=True)
    identical=(output/'measurements.json').read_bytes()==(HERE/'measurements.json').read_bytes()
    if not identical:raise ValueError('remeasured JSON differs; no numerical acceptance inferred')
    result=dict(archive_verified=True,verified_entries=len(seen),measurements_byte_identical=True,
        numerical_acceptance=False,all_new_linear_matrices_archived=False,measurement_source='sealed isolated source snapshot')
    (output/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    print(json.dumps(run(a.output.resolve()),indent=2))
