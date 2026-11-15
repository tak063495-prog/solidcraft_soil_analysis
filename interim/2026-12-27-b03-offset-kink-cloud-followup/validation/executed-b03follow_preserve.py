import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
BASE = '672824272e2baa1fe6ea776c8fc5889eacfc9a15'
ALLOWED = {'interim/2026-12-10-ci-p601-inventory-fix/inventory.toml',
           'interim/2026-12-06-p6-01-status-reconciliation/audit_status.py'}
p = argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args()
started=time.monotonic();rows={};missing=[];mismatches=[];changed=[]
for entry in subprocess.check_output(['git','ls-tree','-rz',BASE],cwd=ROOT).split(b'\0'):
    if not entry:continue
    meta,name=entry.split(b'\t',1);mode,kind,blob=meta.decode().split();name=name.decode();path=ROOT/name
    if kind!='blob':raise ValueError('unexpected tree entry')
    if not path.exists() and not path.is_symlink():
        expected=name.startswith('interim/2026-11-10-nc-postfold-branches/evidence.tar.gz.part')
        missing.append(dict(path=name,git_blob_sha1=blob,existing_ci_exclusion=expected))
        if not expected:mismatches.append(name)
        continue
    if path.is_symlink():
        raw=os.readlink(path).encode();size=len(raw);h1=hashlib.sha1(b'blob '+str(size).encode()+b'\0'+raw);h256=hashlib.sha256(raw)
    else:
        size=path.stat().st_size;h1=hashlib.sha1(b'blob '+str(size).encode()+b'\0');h256=hashlib.sha256()
        with path.open('rb') as stream:
            for chunk in iter(lambda:stream.read(4*1024*1024),b''):h1.update(chunk);h256.update(chunk)
    executable_matches=path.is_symlink() or bool(path.stat().st_mode&0o111)==(mode=='100755')
    equal=h1.hexdigest()==blob and executable_matches
    rows[name]=dict(bytes=size,sha256=h256.hexdigest(),git_blob_sha1=h1.hexdigest(),baseline_blob=blob,executable_matches=executable_matches,unchanged=equal)
    if not equal:
        changed.append(name)
        if name not in ALLOWED:mismatches.append(name)
result=dict(base_commit=BASE,baseline_paths=len(rows)+len(missing),present_checked=len(rows),
            changed_original_paths=changed,allowed_original_changes=sorted(ALLOWED),
            unexpected_changes=mismatches,absent_legacy_parts=missing,
            legacy_git_tree_entries_retained=True,all_other_present_original_files_byte_exact=not mismatches,
            integrity_ok=not mismatches,elapsed_seconds=time.monotonic()-started,files=rows)
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('files','absent_legacy_parts')}),flush=True)
