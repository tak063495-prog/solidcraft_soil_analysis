"""Lossless split ZIP with SHA for every part/member; verification is not solving."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import zipfile
from cases import write


def digest(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def verify(manifest,scratch):
    info=json.loads(manifest.read_text());scratch.mkdir(parents=True,exist_ok=False);archive=scratch/'evidence.zip'
    with archive.open('xb') as dst:
        for part in info['parts']:
            p=manifest.parent/part['name'];assert p.stat().st_size==part['bytes'] and digest(p)==part['sha256']
            with p.open('rb') as src:shutil.copyfileobj(src,dst)
    assert archive.stat().st_size==info['zip_bytes'] and digest(archive)==info['zip_sha256']
    with zipfile.ZipFile(archive) as z:
        assert set(z.namelist())==set(info['entries']) and len(z.namelist())==len(info['entries'])
        for name,expected in info['entries'].items():
            h=hashlib.sha256();size=0
            with z.open(name) as src:
                for block in iter(lambda:src.read(1024*1024),b''):h.update(block);size+=len(block)
            assert size==expected['bytes'] and h.hexdigest()==expected['sha256'],name
    write(scratch/'verification.json',dict(parts=len(info['parts']),entries=len(info['entries']),all_part_zip_and_member_sha256_verified=True,new_solve=False))


def package(work,destination,scratch,labels):
    scratch.mkdir(parents=True,exist_ok=False);archive=scratch/'evidence.zip';selected=[]
    for label in labels:
        for p in (work/label).rglob('*'):
            if p.is_file() and '__pycache__' not in p.parts:selected.append((p,'runs/'+p.relative_to(work).as_posix()))
    for p in Path(__file__).parent.iterdir():
        if p.suffix in ('.py','.rs') or p.name=='SPEC.ja.md':selected.append((p,'supplemental-sources/'+p.name))
    entries={}
    with zipfile.ZipFile(archive,'x',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p,name in sorted(selected,key=lambda v:v[1]):
            assert name not in entries
            entries[name]=dict(bytes=p.stat().st_size,sha256=digest(p));z.write(p,name)
    parts=[]
    with archive.open('rb') as src:
        i=0
        while block:=src.read(24*1024*1024):
            name=f'evidence.zip.part{i:03d}';p=destination/name
            with p.open('xb') as dst:dst.write(block)
            parts.append(dict(name=name,bytes=len(block),sha256=hashlib.sha256(block).hexdigest()));i+=1
    write(destination/'evidence-manifest.json',dict(zip_bytes=archive.stat().st_size,zip_sha256=digest(archive),parts=parts,entries=entries))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--work',type=Path);p.add_argument('--destination',type=Path);p.add_argument('--scratch',type=Path,required=True);p.add_argument('--labels',nargs='+');p.add_argument('--verify',type=Path);a=p.parse_args()
    if a.verify:verify(a.verify.resolve(),a.scratch.resolve())
    else:package(a.work.resolve(),a.destination.resolve(),a.scratch.resolve(),a.labels)
