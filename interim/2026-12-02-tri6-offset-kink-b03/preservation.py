"""Check existing Git content by streaming object hashing, without rewriting old snapshots."""
import argparse,hashlib,os,subprocess
from pathlib import Path
from b03_cases import ROOT,write
BASE='0dad058b00306b35467fa22b3dbaffcf52a0d851'
DOCS={'docs/P6_STATUS.ja.md','docs/STRUCTURE_WATER_JOINT_REMAINDERS.ja.md','interim/2026-11-16-issue-triage/todo.ja.md'}


def run(output):
    tree=subprocess.check_output(['git','ls-tree','-rz',BASE],cwd=ROOT).split(b'\0');same=[];omitted=[];prefix=[]
    for line in tree:
        if not line:continue
        meta,rawpath=line.split(b'\t');mode,kind,oid=meta.decode().split();relative=rawpath.decode();p=ROOT/relative
        if not p.exists():omitted.append(dict(path=relative,git_blob_sha1=oid,mode=mode));continue
        if relative in DOCS:
            before=subprocess.check_output(['git','show',BASE+':'+relative],cwd=ROOT);after=p.read_bytes();assert after.startswith(before)
            prefix.append(dict(path=relative,base_bytes=len(before),appended_bytes=len(after)-len(before),original_prefix_sha256=hashlib.sha256(before).hexdigest()));continue
        size=p.stat().st_size;h=hashlib.sha1(f'blob {size}\0'.encode())
        with p.open('rb') as f:
            for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
        assert h.hexdigest()==oid,relative
        assert bool(p.stat().st_mode&0o111)==(mode=='100755'),relative
        same.append(dict(path=relative,git_blob_sha1=oid,bytes=size,mode=mode))
    result=dict(passed=True,base=BASE,existing_tree_paths=len(same)+len(omitted)+len(prefix),
                present_unchanged_paths=len(same),present_unchanged_bytes=sum(r['bytes'] for r in same),unchanged=same,
                document_prefixes=prefix,sparse_omitted_paths=omitted,sparse_omitted_count=len(omitted),
                note='Omitted contents were not reverified. Existing Git blob/mode is preserved; old hash/snapshot files and failures are not rewritten.')
    write(output,result)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args();run(a.output)
