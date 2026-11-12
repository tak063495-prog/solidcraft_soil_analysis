"""Original Git paths plus a byte/mode check for locally present files only."""
import argparse
import hashlib
import os
from pathlib import Path
import subprocess
from common import ROOT,read,write
from freeze_sources import manifest
BASE='1885f1f374ca417f6ceede181db2621015ccc0b8'

def audit(output):
    entries=subprocess.check_output(['git','ls-tree','-r','-z',BASE],cwd=ROOT)
    present=[];omitted=[];changed=[]
    for entry in entries.split(b'\0'):
        if not entry:continue
        meta,name=entry.split(b'\t',1);mode,kind,blob=meta.decode().split();name=name.decode();p=ROOT/name
        if not p.exists():omitted.append(name);continue
        h=hashlib.sha1();s=hashlib.sha256()
        size=len(os.readlink(p).encode()) if p.is_symlink() else p.stat().st_size
        h.update(f'blob {size}\0'.encode())
        if p.is_symlink():
            data=os.readlink(p).encode();h.update(data);s.update(data)
        else:
            with p.open('rb') as f:
                while data:=f.read(8*1024*1024):h.update(data);s.update(data)
        same=h.hexdigest()==blob;mode_same=p.is_symlink() if mode=='120000' else bool(p.stat().st_mode&0o111)==(mode=='100755')
        record=dict(path=name,bytes=size,git_blob=blob,working_blob=h.hexdigest(),sha256=s.hexdigest(),bytes_equal=same,mode_equal=mode_same)
        present.append(record)
        if not same or not mode_same:changed.append(record)
    before=read(output.parent/'validation/source-snapshot.json');after=manifest()
    result=dict(base_commit=BASE,original_paths=len(present)+len(omitted),present_count=len(present),present_bytes=sum(r['bytes'] for r in present),
                all_present_bytes_and_modes_equal=not changed,changed=changed,omitted_count=len(omitted),omitted_paths=omitted,omitted_bytes_revalidated=False,
                protected_source_snapshot_count=len(before),protected_source_snapshot_before_after_equal=before==after,source_snapshot_after=after,present_files=present)
    write(output,result);assert not changed and before==after

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args();audit(a.output.resolve())
