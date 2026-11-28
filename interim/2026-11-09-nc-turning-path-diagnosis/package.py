"""Package final and exploratory raw evidence without including executables."""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile

HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def package(raw):
    final_roots=['material-requests.json','field-records.json','material-results-final.json',
        'material-final-run.log','material-final-process.json','material-build-final.log',
        'final-source-sha256.json','final-binary-sha256.json','environment.json',
        'source-reproduction.json','source-reproduction-final.log','source-after.json',
        'data-reproduction.json','derivative-measurements.json','assembled-matrices.json',
        'derivative-guarded.log','audit-published.log']
    final={n:raw/n for n in final_roots}
    for folder in ['isolated','start','checks','path-diagnostic-coarse','path-diagnostic-fine','audit-published']:
        final.update({p.relative_to(raw).as_posix():p for p in (raw/folder).rglob('*') if p.is_file()})
    for folder in ['fresh-head-midpoint/input','fresh-head-midpoint/job']:
        files=(raw/'prior'/folder).glob('*.json') if folder.endswith('input') else [raw/'prior'/folder/'recovery-checkpoint-0.json']
        final.update({p.relative_to(raw).as_posix():p for p in files})
    for n in ['fresh-head-midpoint/experiment.json','fresh-head-midpoint/process-observation.json']:
        final['prior/'+n]=raw/'prior'/n
    for folder in ['last-accepted','terminal-rejected']:
        for p in (raw/'prior/experiment/prior/evidence/diagnostic'/folder).glob('*'):
            if p.is_file():final[p.relative_to(raw).as_posix()]=p
    prior_manifest=HERE.parent/'2026-11-09-nc-strict-removal/archive-sha256.json'
    final['prior-archive-sha256.json']=prior_manifest
    for p in HERE.iterdir():
        if p.suffix in ('.py','.rs'):final['study-source/'+p.name]=p
    development={p.name:p for p in raw.iterdir() if p.is_file() and p.name not in final and not p.name.startswith('solid-package') and not p.name.endswith(('.tar.xz','.zip'))}
    paths=sorted(p for p in raw.glob('path-*') if p.is_dir() and p.name not in ['path-diagnostic-coarse','path-diagnostic-fine'])
    for folder in paths:
        development.update({p.relative_to(raw).as_posix():p for p in folder.rglob('*') if p.is_file()})
    groups={'evidence.tar.xz':final,'development-evidence.tar.xz':development}
    manifests={}
    for name,files in groups.items():
        destination=HERE/name
        entries={}
        if destination.exists():
            with tarfile.open(destination,'r:xz') as z:
                if set(z.getnames())!=set(files):raise ValueError('existing archive inventory differs')
                for n,p in sorted(files.items()):
                    entries[n]=sha(p)
                    if hashlib.sha256(z.extractfile(n).read()).hexdigest()!=entries[n]:raise ValueError('existing archive differs: '+n)
        else:
            with tarfile.open(destination,'x:xz',preset=3) as z:
                for n,p in sorted(files.items()):
                    if p.suffix in ['.exe','.rlib','.rmeta'] or not p.is_file():raise ValueError('unexpected evidence file')
                    entries[n]=sha(p);z.add(p,arcname=n,recursive=False)
        if destination.stat().st_size>=100*1024*1024:raise ValueError('split this archive before committing')
        manifests[name]=dict(archive_sha256=sha(destination),bytes=destination.stat().st_size,entries=entries)
        print(name,len(entries),destination.stat().st_size,manifests[name]['archive_sha256'],flush=True)
    (HERE/'archive-sha256.json').write_text(json.dumps(dict(archives=manifests),indent=2)+'\n')
    return manifests
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--raw',type=Path,required=True);a=p.parse_args();package(a.raw.resolve())
