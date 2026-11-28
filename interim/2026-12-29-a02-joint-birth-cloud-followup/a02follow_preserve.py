"""Verify the original Git tree without lazy-fetching or rewriting any blob."""
import argparse, hashlib, json, os
from pathlib import Path
import subprocess

BASE='11df919ad053fc742fa40cb986953d8f53de28db'
ROOT=Path(__file__).resolve().parents[2]

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--after-registration',action='store_true');args=p.parse_args()
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
    expected_registration_changes=False
    if args.after_registration:
        inv='interim/2026-12-10-ci-p601-inventory-fix/inventory.toml'
        auditor='interim/2026-12-06-p6-01-status-reconciliation/audit_status.py'
        original_inventory=subprocess.check_output(['git','show',BASE+':'+inv],cwd=ROOT,env=env)
        original_auditor=subprocess.check_output(['git','show',BASE+':'+auditor],cwd=ROOT,env=env)
        actual_inventory=(ROOT/inv).read_bytes();actual_auditor=(ROOT/auditor).read_bytes()
        before=hashlib.sha256(original_inventory).hexdigest().encode();after=hashlib.sha256(actual_inventory).hexdigest().encode()
        expected_registration_changes=(set(errors)=={inv,auditor} and actual_inventory.startswith(original_inventory)
            and actual_auditor.replace(after,before)==original_auditor)
    report=dict(base=BASE,expected_registration_changes=expected_registration_changes,files=rows,excluded_postfold=missing,unexpected_missing=unexpected_missing,changed_old_paths=errors,
        all_materialized_old_bytes_equal_git=not errors and not unexpected_missing)
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    if (errors and not expected_registration_changes) or unexpected_missing:raise ValueError((errors,unexpected_missing))

if __name__=='__main__':main()
