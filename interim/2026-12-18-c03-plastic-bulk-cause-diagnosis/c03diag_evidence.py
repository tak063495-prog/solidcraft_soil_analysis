"""Seal and verify raw C-03 evidence. Verification never starts a solve."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import zipfile


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def write(path,value):
    with path.open('x',encoding='utf-8',newline='\n') as stream:
        json.dump(value,stream,indent=2,ensure_ascii=False,allow_nan=False);stream.write('\n')


def pack(root,output):
    entries={};parts=[]
    with tempfile.TemporaryDirectory(dir=output.parent) as tmp:
        zp=Path(tmp)/'evidence.zip'
        with zipfile.ZipFile(zp,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as archive:
            for p in sorted(root.rglob('*')):
                if not p.is_file() or '__pycache__' in p.parts:continue
                name=p.relative_to(root).as_posix();before=sha(p);archive.write(p,name)
                if sha(p)!=before:raise ValueError('evidence changed during sealing: '+name)
                entries[name]=dict(size=p.stat().st_size,sha256=before)
        with zp.open('rb') as stream:
            for i in range((zp.stat().st_size+19_999_999)//20_000_000):
                p=output/f'evidence.zip.part{i:03}'
                with p.open('xb') as f:f.write(stream.read(20_000_000))
                parts.append(dict(name=p.name,size=p.stat().st_size,sha256=sha(p)))
        write(output/'archive-index.json',dict(schema='c03-cause-sealed-evidence/1',zip_sha256=sha(zp),
              zip_size=zp.stat().st_size,parts=parts,entries=entries,raw_evidence=True,new_solve_on_verify=False))


def verify(source,output):
    output.mkdir();index=json.loads((source/'archive-index.json').read_text(encoding='utf-8'));zp=output/'evidence.zip'
    with zp.open('xb') as f:
        for part in index['parts']:
            p=source/part['name']
            if p.stat().st_size!=part['size'] or sha(p)!=part['sha256']:raise ValueError('sealed part differs: '+p.name)
            with p.open('rb') as stream:shutil.copyfileobj(stream,f)
    if zp.stat().st_size!=index['zip_size'] or sha(zp)!=index['zip_sha256']:raise ValueError('ZIP differs')
    with zipfile.ZipFile(zp) as archive:
        if sorted(archive.namelist())!=sorted(index['entries']):raise ValueError('ZIP entry inventory differs')
        for name,item in index['entries'].items():
            p=Path(name)
            if p.is_absolute() or '..' in p.parts:raise ValueError('unsafe member')
            with archive.open(name) as f:
                h=hashlib.sha256();size=0
                for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk);size+=len(chunk)
            if size!=item['size'] or h.hexdigest()!=item['sha256']:raise ValueError('entry differs: '+name)
    write(output/'verification.json',dict(parts=len(index['parts']),entries=len(index['entries']),all_sha256_verified=True,
          new_solve=False,checkpoint_restart=False,passed=False))


def main():
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=['pack','verify'],required=True)
    p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.mode=='pack':pack(a.source,a.output)
    else:verify(a.source,a.output)

if __name__=='__main__':main()
