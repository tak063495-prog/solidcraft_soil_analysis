"""Verify the original Git tree without lazy-fetching or rewriting any blob."""
import argparse, hashlib, json, os
from pathlib import Path
import subprocess

BASE='b383aa007afd99e74f5988060335b2e07c8d341f'
ROOT=Path(__file__).resolve().parents[2]

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    env=dict(os.environ,GIT_NO_LAZY_FETCH='1');rows=[];missing=[];errors=[]
    for line in subprocess.check_output(['git','ls-tree','-r',BASE],cwd=ROOT,text=True,env=env).splitlines():
        meta,name=line.split('\t');mode,kind,expected=meta.split()
        if kind!='blob':continue
        path=ROOT/name
        if not path.is_file():missing.append(name);continue
        size=path.stat().st_size;s=hashlib.sha256();g=hashlib.sha1(('blob '+str(size)+'\0').encode())
        with path.open('rb') as stream:
            while data:=stream.read(1024*1024):s.update(data);g.update(data)
        if g.hexdigest()!=expected:errors.append(name)
        rows.append(dict(path=name,sha256=s.hexdigest(),git_blob_sha1=g.hexdigest(),size=size,mode=mode))
    unexpected_missing=[name for name in missing if not name.startswith('interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part')]
    report=dict(base=BASE,files=rows,excluded_postfold=missing,unexpected_missing=unexpected_missing,changed_old_paths=errors,
        all_materialized_old_bytes_equal_git=not errors and not unexpected_missing)
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    if errors or unexpected_missing:raise ValueError((errors,unexpected_missing))

if __name__=='__main__':main()
