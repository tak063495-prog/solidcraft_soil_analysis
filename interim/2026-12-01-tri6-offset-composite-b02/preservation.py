"""Original tree identities + present-byte audit; omitted files are not certified."""
import argparse
import hashlib
from pathlib import Path
import subprocess
from cases import ROOT,write
BASE='e940bad4015652a6cc713d0bdcde35f062c8f546'
DOCUMENTS={'docs/P6_STATUS.ja.md','docs/STRUCTURE_WATER_JOINT_REMAINDERS.ja.md','interim/2026-11-16-issue-triage/todo.ja.md'}


def audit(output):
    tree=subprocess.check_output(['git','ls-tree','-r',BASE],cwd=ROOT,text=True)
    present=[];omitted=[];changed=[];documents=[]
    for line in tree.splitlines():
        meta,name=line.split('\t');mode,kind,expected=meta.split();p=ROOT/name
        if not p.is_file():omitted.append(name);continue
        if name in DOCUMENTS:
            original=subprocess.check_output(['git','show',BASE+':'+name],cwd=ROOT);current=p.read_bytes()
            documents.append(dict(path=name,original_git_blob=expected,original_bytes=len(original),original_sha256=hashlib.sha256(original).hexdigest(),
                                  current_bytes=len(current),current_sha256=hashlib.sha256(current).hexdigest(),original_prefix_preserved=current.startswith(original)))
            continue
        size=p.stat().st_size;h=hashlib.sha1(f'blob {size}\0'.encode());sha=hashlib.sha256()
        with p.open('rb') as f:
            for block in iter(lambda:f.read(1024*1024),b''):h.update(block);sha.update(block)
        same=h.hexdigest()==expected
        mode_same=bool(p.stat().st_mode&0o111)==(mode=='100755')
        record=dict(path=name,bytes=size,git_blob_sha1=expected,working_blob_sha1=h.hexdigest(),sha256=sha.hexdigest(),bytes_match=same,executable_mode_match=mode_same)
        present.append(record)
        if not same or not mode_same:changed.append(record)
    write(output,dict(base=BASE,old_git_paths=len(tree.splitlines()),present_protected_files=len(present),present_protected_bytes=sum(x['bytes'] for x in present),
                      all_present_protected_bytes_and_modes_match=not changed,changed_protected=changed,append_only_documents=documents,
                      all_original_document_prefixes_preserved=all(x['original_prefix_preserved'] for x in documents),
                      omitted_paths=omitted,omitted_count=len(omitted),omitted_contents_revalidated=False,present_files=present))
    assert not changed and all(x['original_prefix_preserved'] for x in documents)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args();audit(a.output)
