"""Verify existing evidence; prepare unchanged observer and input-only history."""
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run(out):
    out.mkdir(parents=True,exist_ok=False)
    prior=HERE.parent/'2026-11-09-nc-turning-path-diagnosis'
    manifest=json.loads((prior/'archive-sha256.json').read_text())['archives']['evidence.tar.xz']
    archive=prior/'evidence.tar.xz'
    if sha(archive)!=manifest['archive_sha256']:
        raise ValueError('existing evidence archive differs')
    seen=set()
    with tarfile.open(archive) as z:
        for ent in z:
            name=Path(ent.name)
            if not ent.isfile() or name.is_absolute() or '..' in name.parts or ent.name in seen:
                raise ValueError('unsafe or duplicate evidence entry')
            b=z.extractfile(ent).read()
            if hashlib.sha256(b).hexdigest()!=manifest['entries'][ent.name]:
                raise ValueError('existing evidence entry differs')
            seen.add(ent.name)
            if ent.name in ['field-records.json','assembled-matrices.json'] or ent.name.startswith('prior/fresh-head-midpoint/input/') or ent.name=='prior/fresh-head-midpoint/job/recovery-checkpoint-0.json':
                q=out/'origin'/ent.name;q.parent.mkdir(parents=True,exist_ok=True);q.write_bytes(b)
    if seen!=set(manifest['entries']):raise ValueError('existing evidence inventory differs')
    source=out/'observer';source.mkdir();shutil.copytree(ROOT/'src',source/'src')
    for n in ['Cargo.toml','Cargo.lock','rust-toolchain.toml']:shutil.copyfile(ROOT/n,source/n)
    for n,a in [('mcc.rs','mcc-diagnostic-addition.rs'),('mcc_midpoint.rs','midpoint-diagnostic-addition.rs')]:
        p=source/'src'/n;p.write_bytes(p.read_bytes().replace(b'\r\n',b'\n')+(prior/a).read_bytes())
    (source/'examples').mkdir();shutil.copyfile(prior/'material_diagnostic.rs',source/'examples/material_diagnostic.rs')
    subprocess.run(['cargo','fmt','--manifest-path',str(source/'Cargo.toml')],check=True)
    expected=json.loads((HERE.parent/'2026-11-09-nc-pseudo-arclength/archive-sha256.json').read_text())['entries']
    for p in source.rglob('*'):
        if p.is_file() and sha(p)!=expected['isolated/'+p.relative_to(source).as_posix()]:
            raise ValueError('observer source differs: '+str(p))
    subprocess.run(['cargo','build','--locked','--offline','--release','--manifest-path',str(source/'Cargo.toml'),
                    '--target-dir',str(out/'observer-target'),'--example','material_diagnostic'],check=True)
    subprocess.run([sys.executable,str(HERE.parent/'2026-11-09-nc-pseudo-arclength/capture_reference_history.py'),'--output',str(out/'reference400'),
                    '--steps','400'],check=True)
    if not (out/'reference400/unload.json').is_file():raise ValueError('input-only prehistory incomplete')
    (out/'preparation.json').write_text(json.dumps(dict(archive_sha256=manifest['archive_sha256'],verified_entries=len(seen),
        observer_sources_matched=True,passed=False,prehistory_is_not_full_release=True),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    run(p.parse_args().output.resolve())
