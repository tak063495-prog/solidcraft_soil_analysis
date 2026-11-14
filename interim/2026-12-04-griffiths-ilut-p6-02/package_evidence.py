"""Split and verify exact raw evidence. Never counts verification as a solve."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import zipfile


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def write(p,v):
    with p.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')


def create(work,out):
    archive=work/'p602-evidence.zip';members=[]
    groups=['base','final','frozen','restarts','source-before','source-final',
            'validation-rust','validation-python','python-repairs','preliminary','pc-smoke','observer-build']
    paths=[]
    for name in groups:
        paths += [(p,name+'/'+p.relative_to(work/name).as_posix()) for p in (work/name).rglob('*')
                  if p.is_file() and '__pycache__' not in p.parts and 'target' not in p.relative_to(work/name).parts]
    paths += [(p,'logs-and-contracts/'+p.name) for p in work.iterdir() if p.is_file()
              and p.suffix in ('.log','.json','.txt','.py')
              and p.name not in ('package.log','artifact-verify.log','collect-archive.log')]
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
        for p,name in sorted(paths,key=lambda x:x[1]):
            h=sha(p);size=p.stat().st_size;z.write(p,name)
            assert sha(p)==h and p.stat().st_size==size
            members.append(dict(path=name,bytes=size,sha256=h))
    parts=[]
    with archive.open('rb') as f:
        while data:=f.read(8*1024*1024):
            p=out/f'evidence.zip.part{len(parts):03d}'
            with p.open('xb') as target:target.write(data)
            parts.append(dict(path=p.name,bytes=len(data),sha256=sha(p)))
    manifest=dict(schema='p602-ilut-evidence/1',archive_bytes=archive.stat().st_size,
        archive_sha256=sha(archive),parts=parts,members=members,
        verification_is_not_a_new_solve=True)
    write(out/'evidence-manifest.json',manifest)
    return manifest


def verify(manifest,scratch):
    d=json.loads(manifest.read_text());scratch.mkdir(parents=True,exist_ok=False)
    archive=scratch/'combined.zip'
    with archive.open('xb') as f:
        for part in d['parts']:
            p=manifest.parent/part['path']
            assert p.stat().st_size==part['bytes'] and sha(p)==part['sha256']
            with p.open('rb') as src:shutil.copyfileobj(src,f)
    assert archive.stat().st_size==d['archive_bytes'] and sha(archive)==d['archive_sha256']
    members={v['path']:v for v in d['members']}
    with zipfile.ZipFile(archive) as z:
        assert len(z.namelist())==len(members) and set(z.namelist())==set(members)
        for name,expected in members.items():
            p=Path(name)
            assert not p.is_absolute() and '..' not in p.parts
            h=hashlib.sha256();size=0
            with z.open(name) as f:
                for b in iter(lambda:f.read(1024*1024),b''):h.update(b);size+=len(b)
            assert size==expected['bytes'] and h.hexdigest()==expected['sha256'],name
        z.extractall(scratch/'raw')
    result=dict(passed=True,parts_verified=len(d['parts']),members_verified=len(members),
                archive_sha256=d['archive_sha256'],verification_is_not_a_new_solve=True)
    write(scratch/'verification.json',result)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--work',type=Path);p.add_argument('--output',type=Path)
    p.add_argument('--verify',type=Path);p.add_argument('--scratch',type=Path);a=p.parse_args()
    print(json.dumps(verify(a.verify,a.scratch) if a.verify else create(a.work,a.output)))
