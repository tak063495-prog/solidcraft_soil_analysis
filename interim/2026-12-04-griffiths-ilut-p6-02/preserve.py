"""Compare all present old files with Git blobs; never fetch omitted evidence."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

BASE='2a716084dc66754b681f29e062e5d9155d505647'
PREFIX={'docs/P6_STATUS.ja.md','interim/2026-11-16-issue-triage/todo.ja.md'}
CHANGED={'src/linalg.rs','src/checkpoint.rs'}


def preserve(root):
    tree=subprocess.check_output(['git','ls-tree','-r',BASE],cwd=root,text=True)
    rows=[];omitted=[];prefixes=[];changes=[]
    for line in tree.splitlines():
        meta,name=line.split('\t');mode,kind,blob=meta.split();p=root/name
        if not p.is_file():omitted.append(dict(path=name,mode=mode,blob=blob));continue
        actual_mode='100755' if p.stat().st_mode&0o111 else '100644';assert actual_mode==mode,name
        if name in PREFIX|CHANGED:
            original=subprocess.check_output(['git','show',BASE+':'+name],cwd=root)
            current=p.read_bytes()
            r=dict(path=name,original_sha256=hashlib.sha256(original).hexdigest(),current_sha256=hashlib.sha256(current).hexdigest(),original_bytes=len(original),current_bytes=len(current))
            if name in PREFIX:
                assert current.startswith(original),name
                r['all_original_bytes_preserved_as_prefix']=True;prefixes.append(r)
            else:changes.append(r)
        else:
            size=p.stat().st_size;h=hashlib.sha1(f'blob {size}\0'.encode())
            with p.open('rb') as f:
                for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
            assert h.hexdigest()==blob,name
            rows.append(dict(path=name,bytes=size,blob=blob,mode=mode))
    assert {x['path'] for x in changes}==CHANGED
    return dict(passed=True,baseline=BASE,baseline_paths=len(tree.splitlines()),
        present_unchanged_files=len(rows),present_unchanged_bytes=sum(x['bytes'] for x in rows),
        unchanged=rows,prefix_only_documents=prefixes,explicit_production_changes=changes,
        omitted_count=len(omitted),omitted_original_mode_blob=omitted,omitted_contents_not_reaudited=True,
        old_hash_snapshot_failure_and_reference_files_unchanged=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    with a.output.open('x',encoding='utf-8') as f:json.dump(preserve(a.root),f,indent=2);f.write('\n')
