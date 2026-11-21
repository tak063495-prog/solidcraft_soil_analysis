"""Verify old archive bytes and extract only local-fold replay evidence.

No solver is run. Original archives, checkpoints, manifests and snapshots are
read-only. Output paths must be new. All members are hashed, not just extracts.
"""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile

ROOT=Path(__file__).resolve().parents[2]
PRIOR=ROOT/'interim/2026-11-09-nc-pseudo-arclength'
PATHS=('native-coarse','native-fine','independent-map-arc-coarse','independent-map-arc-fine','independent-map-arc-material8')


def sha(path):
    result=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):result.update(chunk)
    return result.hexdigest()


def selected(name):
    if name.startswith('prior/prior/fresh-head-midpoint/input/') or name in (
            'prior/field-records.json','prior/assembled-matrices.json',
            'prior/prior/fresh-head-midpoint/job/recovery-checkpoint-0.json',
            'independent-input-rk4-fine/unload.json','independent-input-rk4-fine/reference.json',
            'independent-arc-history400/samples.json'):return True
    return any(name==p+'/'+f for p in PATHS for f in ('samples.json','linear-systems.json','events.json','controls.json','seed.json','report.json'))


def run(output):
    output.mkdir(parents=True,exist_ok=False)
    manifest=json.loads((PRIOR/'archive-sha256.json').read_text(encoding='utf-8'))
    archive=output/'original.tar.xz';combined=hashlib.sha256()
    with archive.open('xb') as target:
        for part in manifest['parts']:
            path=PRIOR/part['name']
            if path.stat().st_size!=part['bytes'] or sha(path)!=part['sha256']:raise ValueError('old part checksum differs')
            with path.open('rb') as source:
                for chunk in iter(lambda:source.read(1024*1024),b''):target.write(chunk);combined.update(chunk)
    if combined.hexdigest()!=manifest['archive_sha256']:raise ValueError('old archive checksum differs')
    destination=output/'evidence';destination.mkdir();names=set();extracted={}
    with tarfile.open(archive,'r|xz') as stream:
        for member in stream:
            if not member.isfile() or member.name in names or member.name not in manifest['entries']:
                raise ValueError('unexpected old archive member')
            names.add(member.name);hasher=hashlib.sha256();file=None
            if selected(member.name):
                path=destination/member.name;path.parent.mkdir(parents=True,exist_ok=True);file=path.open('xb')
            try:
                source=stream.extractfile(member)
                for chunk in iter(lambda:source.read(1024*1024),b''):
                    hasher.update(chunk)
                    if file:file.write(chunk)
            finally:
                if file:file.close()
            if hasher.hexdigest()!=manifest['entries'][member.name]:raise ValueError('old member checksum differs')
            if selected(member.name):extracted[member.name]=hasher.hexdigest()
    if names!=set(manifest['entries']):raise ValueError('old member set differs')
    # Only the accepted native material evaluations adjacent to the turning point.
    wanted=set()
    for path in PATHS[:2]:
        samples=json.loads((destination/path/'samples.json').read_text())
        peak=max(range(len(samples)),key=lambda i:samples[i]['progress'])
        steps={r['step'] for r in samples[max(0,peak-2):peak+4]}
        events=json.loads((destination/path/'events.json').read_text())
        for step in steps:
            event=[e for e in events if e['step']==step][-1]
            wanted.add(f'{path}/call-{event["call"]:05d}-result.json')
    with tarfile.open(archive,'r|xz') as stream:
        for member in stream:
            if member.name in wanted:
                raw=stream.extractfile(member).read()
                if hashlib.sha256(raw).hexdigest()!=manifest['entries'][member.name]:raise ValueError('native replay member differs')
                path=destination/member.name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(raw)
                extracted[member.name]=manifest['entries'][member.name]
    if not wanted<=set(extracted):raise ValueError('accepted native GP tangents missing')
    record=dict(archive_sha256=manifest['archive_sha256'],all_old_members_verified=len(names),extracted_sha256=extracted,
                original_evidence_modified=False,solver_run=False)
    (output/'extraction.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();print(json.dumps(run(args.output),indent=2))
