"""Retain raw fresh evidence in Git-sized parts with byte/member SHA manifests."""
import argparse
import hashlib
import json
import zipfile
from pathlib import Path
from common import write

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        while b:=f.read(8*1024*1024):h.update(b)
    return h.hexdigest()

def package(source,destination,scratch):
    scratch.mkdir(parents=True,exist_ok=False)
    archive=scratch/'evidence.zip';entries=[]
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
        for p in sorted(source.rglob('*')):
            if not p.is_file():continue
            name=str(p.relative_to(source)).replace('\\','/');z.write(p,name)
            entries.append(dict(name=name,bytes=p.stat().st_size,sha256=digest(p)))
    parts=[]
    with archive.open('rb') as f:
        index=0
        while b:=f.read(48*1024*1024):
            p=destination/f'evidence.zip.part{index:03d}'
            with p.open('xb') as o:o.write(b)
            parts.append(dict(name=p.name,bytes=len(b),sha256=hashlib.sha256(b).hexdigest()));index+=1
    manifest=dict(passed=False,numerical_acceptance=False,archive_sha256=digest(archive),archive_bytes=archive.stat().st_size,parts=parts,members=entries)
    write(destination/'evidence-manifest.json',manifest)
    verify(destination/'evidence-manifest.json',scratch/'verified')

def verify(manifest_path,scratch):
    manifest=json.loads(manifest_path.read_text());scratch.mkdir(parents=True,exist_ok=False);archive=scratch/'evidence.zip'
    with archive.open('xb') as out:
        for part in manifest['parts']:
            p=manifest_path.parent/part['name'];assert p.stat().st_size==part['bytes'] and digest(p)==part['sha256']
            with p.open('rb') as f:
                while b:=f.read(8*1024*1024):out.write(b)
    assert archive.stat().st_size==manifest['archive_bytes'] and digest(archive)==manifest['archive_sha256']
    with zipfile.ZipFile(archive) as z:
        assert z.namelist()==[e['name'] for e in manifest['members']]
        for entry in manifest['members']:
            h=hashlib.sha256();size=0
            with z.open(entry['name']) as f:
                while b:=f.read(8*1024*1024):size+=len(b);h.update(b)
            assert size==entry['bytes'] and h.hexdigest()==entry['sha256']
    write(scratch/'verification.json',dict(parts=len(manifest['parts']),members=len(manifest['members']),all_sha256_equal=True,new_solve=False))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path);p.add_argument('--destination',type=Path);p.add_argument('--scratch',type=Path,required=True);p.add_argument('--verify',type=Path);a=p.parse_args()
    if a.verify:verify(a.verify.resolve(),a.scratch.resolve())
    else:package(a.source.resolve(),a.destination.resolve(),a.scratch.resolve())
