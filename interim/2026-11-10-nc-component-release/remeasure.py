"""Verify every distributed byte, extract independently and rerun measurements."""
import argparse
import hashlib
import json
import subprocess
import sys
import tarfile
from pathlib import Path

HERE=Path(__file__).resolve().parent


class Parts:
    def __init__(self,paths):self.paths=iter(paths);self.f=None
    def read(self,n=-1):
        if n<0:raise ValueError('bounded read required')
        result=bytearray()
        while len(result)<n:
            if self.f is None:
                try:self.f=next(self.paths).open('rb')
                except StopIteration:break
            b=self.f.read(n-len(result))
            if b:result.extend(b)
            else:self.f.close();self.f=None
        return bytes(result)


def run(out):
    out.mkdir(parents=True,exist_ok=False);m=json.loads((HERE/'archive-sha256.json').read_text());h=hashlib.sha256()
    for part in m['parts']:
        ph=hashlib.sha256();size=0
        with (HERE/part['name']).open('rb') as f:
            while b:=f.read(1024*1024):h.update(b);ph.update(b);size+=len(b)
        if size!=part['bytes'] or ph.hexdigest()!=part['sha256']:raise ValueError('archive part differs')
    if h.hexdigest()!=m['archive_sha256']:raise ValueError('archive digest differs')
    seen=set();root=out/'evidence'
    with tarfile.open(fileobj=Parts([HERE/p['name'] for p in m['parts']]),mode='r|gz') as z:
        for ent in z:
            name=Path(ent.name)
            if not ent.isfile() or name.is_absolute() or '..' in name.parts or ent.name in seen:raise ValueError('unsafe/duplicate entry')
            dst=root/name;dst.parent.mkdir(parents=True,exist_ok=True);eh=hashlib.sha256();size=0
            with z.extractfile(ent) as f,dst.open('xb') as target:
                while b:=f.read(1024*1024):eh.update(b);target.write(b);size+=len(b)
            if dict(sha256=eh.hexdigest(),bytes=size)!=m['entries'][ent.name]:raise ValueError('entry differs')
            seen.add(ent.name)
    if seen!=set(m['entries']):raise ValueError('archive inventory differs')
    for p in HERE.glob('*.py'):
        archived=root/'diagnostic-source'/p.name
        if not archived.is_file() or archived.read_bytes()!=p.read_bytes():raise ValueError('diagnostic source differs')
    result=out/'measurements.json'
    subprocess.run([sys.executable,str(HERE/'measure.py'),'--evidence',str(root),'--output',str(result)],check=True)
    same=result.read_bytes()==(HERE/'measurements.json').read_bytes()
    (out/'remeasurement.json').write_text(json.dumps(dict(verified_entries=len(seen),archive_sha256=h.hexdigest(),
        measurements_byte_identical=same,numerical_solver_executed=False,passed=False),indent=2)+'\n')
    if not same:raise ValueError('remeasurements differ')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    run(p.parse_args().output.resolve())
