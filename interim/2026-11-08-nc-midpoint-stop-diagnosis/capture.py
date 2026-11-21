"""Capture two recorded increments for diagnosis, not a continuation run."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time
import zipfile
from prepare_diagnostic import prepare, ROOT
from measure import measure

HERE=Path(__file__).resolve().parent
PRIOR=ROOT/'interim/2026-11-08-nc-midpoint-resume-and-cut-history'

def digest(data):return hashlib.sha256(data).hexdigest()

def extract(out):
    manifest=json.loads((PRIOR/'archive-sha256.json').read_text())
    archive=PRIOR/'evidence.zip'
    assert digest(archive.read_bytes())==manifest['archive_sha256'],'prior archive changed'
    out.mkdir(parents=True,exist_ok=False)
    checks={}
    with zipfile.ZipFile(archive) as z:
        for name,expected in manifest['entries'].items():
            if name.startswith('nc-'):
                data=z.read(name);assert digest(data)==expected,('prior entry changed',name)
                p=out/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data);checks[name]=expected
    return dict(archive=str(archive.relative_to(ROOT)),archive_sha256=manifest['archive_sha256'],entries=checks)

def command(args,log,records):
    start=time.monotonic()
    with log.open('x') as f:
        r=subprocess.run([str(v) for v in args],stdout=f,stderr=subprocess.STDOUT,timeout=300)
    records.append(dict(command=[str(v) for v in args],exit_code=r.returncode,elapsed_seconds=time.monotonic()-start))
    if r.returncode:raise RuntimeError('diagnostic process failed: '+str(log))


def capture(out):
    out.mkdir(parents=True,exist_ok=False)
    binding=extract(out/'prior');(out/'prior-binding.json').write_text(json.dumps(binding,indent=2))
    prepare(out/'source');records=[]
    command(['cargo','build','--locked','--release','--manifest-path',out/'source/Cargo.toml',
             '--example','nc_stop_probe','--example','nc_midpoint_map_probe'],out/'build.log',records)
    suffix='.exe' if __import__('os').name=='nt' else ''
    exe=out/('source/target/release/examples/nc_stop_probe'+suffix)
    probe=out/('source/target/release/examples/nc_midpoint_map_probe'+suffix)
    dst=out/'diagnostic';dst.mkdir()
    for name,cp,target in [('last-accepted','recovery-checkpoint-1.json','0.9514000415802002'),
                           ('terminal-rejected','recovery-checkpoint-0.json','0.9514002799987793')]:
        command([exe,out/'prior/nc-resumed/input/mesh.json',out/'prior/nc-resumed/job'/cp,target,dst/name],dst/(name+'.log'),records)
    (dst/'process.json').write_text(json.dumps(records,indent=2))
    (dst/'binary-sha256.json').write_text(json.dumps({str(p):digest(p.read_bytes()) for p in (exe,probe)},indent=2))
    return measure(out/'prior',dst,probe,out/'measurements')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    r=capture(a.output.resolve());print(json.dumps({k:r[k] for k in ('passed','history_converged','checkpoint_progress')},indent=2))
