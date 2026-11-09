"""Seal full accepted fields, real linear systems, seeds and source snapshots.

Native trial request/response bytes remain in work; their individual SHA/size
manifest is distributed. All accepted material/contact fields are distributed.
"""
import argparse
import hashlib
import json
import tarfile
from pathlib import Path

HERE=Path(__file__).resolve().parent


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        while b:=f.read(1024*1024):h.update(b)
    return h.hexdigest()


def run(w):
    files={}
    paths=sorted(p for p in w.iterdir() if p.is_dir() and (p/'controls.json').exists())
    if len(paths)!=10 or any(not (p/'report.json').exists() for p in paths):
        raise ValueError('ten completed diagnostic processes required')
    for folder in paths:
        io=folder/'material'
        if io.exists():
            entries={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(io.iterdir())}
            (folder/'material-io-sha256.json').write_text(json.dumps(dict(trial_bytes_included=False,entries=entries),indent=2)+'\n')
        for p in folder.iterdir():
            if p.is_file():files[p.relative_to(w).as_posix()]=p
    for directory in ['origin','reference400','observer']:
        for p in (w/directory).rglob('*'):
            if p.is_file() and '__pycache__' not in p.parts:files[p.relative_to(w).as_posix()]=p
    for p in w.glob('*.log'):files[p.name]=p
    for name in ['numerical-processes.json','preparation.json']:
        if (w/name).exists():files[name]=w/name
    for p in HERE.glob('*.py'):files['diagnostic-source/'+p.name]=p
    entries={name:dict(sha256=sha(p),bytes=p.stat().st_size) for name,p in sorted(files.items())}
    archive=w/'new-evidence.tar.gz'
    with tarfile.open(archive,'w:gz',compresslevel=6) as z:
        for name,p in sorted(files.items()):
            if sha(p)!=entries[name]['sha256']:raise ValueError('source mutated during seal')
            z.add(p,arcname=name,recursive=False)
    manifest=dict(entries=entries,archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,parts=[],
                  all_accepted_fields_included=True,native_trial_io_included=False,binary_included=False)
    with archive.open('rb') as f:
        i=0
        while b:=f.read(16*1024*1024):
            i+=1;p=HERE/f'evidence.tar.gz.part{i:02d}';p.write_bytes(b)
            manifest['parts'].append(dict(name=p.name,bytes=len(b),sha256=hashlib.sha256(b).hexdigest()))
    (HERE/'archive-sha256.json').write_text(json.dumps(manifest,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--evidence',type=Path,required=True)
    run(p.parse_args().evidence.resolve())
