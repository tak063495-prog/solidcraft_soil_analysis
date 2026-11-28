"""Seal raw fresh-solve evidence; replay is explicitly not another fresh solve."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import zipfile
from study import read,write,sha


def pack(runs,source,output,retained):
    output.mkdir(parents=True,exist_ok=True);entries={};files=[]
    for p in sorted(runs.rglob('*')):
        if p.is_file():
            relative=p.relative_to(runs)
            name=('reports/'+p.name if relative.as_posix() in ('measurements.json','completion.json') else 'data/final/'+relative.as_posix())
            files.append((p,name))
    binaries={}
    for root in retained:
        for p in sorted(root.rglob('*')):
            if p.is_file() and p.name!='probe' and p.name!='probe.exe':files.append((p,'data/retained/'+root.name+'/'+p.relative_to(root).as_posix()))
        for binary in [root/'probe',root/'probe.exe']:
            if binary.is_file():binaries.setdefault(sha(binary),binary)
    for digest,p in binaries.items():files.append((p,'data/retained/binaries/'+digest))
    for directory in ['src','tools','tests','examples','interim/2026-11-20-wet-joint-plastic-bulk-integration']:
        for p in sorted((source/directory).rglob('*')):
            if p.is_file() and '__pycache__' not in p.parts and p.suffix in ('.py','.rs','.json','.toml','.txt'):files.append((p,'source/'+p.relative_to(source).as_posix()))
    for name in ['Cargo.toml','Cargo.lock','rust-toolchain.toml','requirements-reference.lock.txt','requirements.lock.txt']:
        p=source/name;files.append((p,'source/'+name))
    assert len(files)==len({name for _,name in files})
    with tempfile.TemporaryDirectory(dir=output.parent) as temporary:
        zp=Path(temporary)/'evidence.zip'
        with zipfile.ZipFile(zp,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as archive:
            for p,name in files:
                before=sha(p);archive.write(p,name);assert sha(p)==before
                entries[name]=dict(sha256=before,size=p.stat().st_size)
        parts=[]
        with zp.open('rb') as stream:
            for index in range((zp.stat().st_size+20_000_000-1)//20_000_000):
                target=output/f'evidence.zip.part{index:03}'
                with target.open('xb') as f:f.write(stream.read(20_000_000))
                parts.append(dict(name=target.name,sha256=sha(target),size=target.stat().st_size))
        write(output/'archive-sha256.json',dict(zip_sha256=sha(zp),zip_size=zp.stat().st_size,parts=parts,entries=entries,
              raw_evidence_not_replaced_by_reaudit=True,source_origin='specified Git blobs plus new observation code'))


def replay(archive_folder,target):
    assert not target.exists();target.mkdir(parents=True)
    index=read(archive_folder/'archive-sha256.json');zp=target/'evidence.zip'
    with zp.open('xb') as out:
        for item in index['parts']:
            p=archive_folder/item['name'];assert p.stat().st_size==item['size'] and sha(p)==item['sha256']
            with p.open('rb') as stream:shutil.copyfileobj(stream,out)
    assert zp.stat().st_size==index['zip_size'] and sha(zp)==index['zip_sha256']
    root=target/'unpacked';root.mkdir()
    with zipfile.ZipFile(zp) as archive:
        assert set(archive.namelist())==set(index['entries']) and len(archive.namelist())==len(index['entries'])
        for name,item in index['entries'].items():
            relative=Path(name);assert not relative.is_absolute() and '..' not in relative.parts
            p=root/relative;p.parent.mkdir(parents=True,exist_ok=True)
            with archive.open(name) as src,p.open('xb') as dst:shutil.copyfileobj(src,dst)
            assert p.stat().st_size==item['size'] and sha(p)==item['sha256']
    script=root/'source/interim/2026-11-20-wet-joint-plastic-bulk-integration/measure.py';result=target/'remeasured.json'
    routes=['baseline','material-half','material-quarter','contact-half','contact-quarter','single-partition-control']
    bootstrap="import pathlib,runpy,sys; script=sys.argv[1]; sys.path.insert(0,str(pathlib.Path(script).parent)); sys.argv=sys.argv[1:]; runpy.run_path(script,run_name='__main__')"
    cmd=[sys.executable,'-I','-X','utf8','-c',bootstrap,str(script),'--output',str(root/'data/final'),'--steps','32','--routes',*routes,'--report',str(result)]
    with (target/'replay.log').open('x') as log:process=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
    assert process.returncode==0,(process.returncode,target/'replay.log')
    assert read(result)==read(root/'reports/measurements.json'),'sealed numerical remeasurement differs'
    write(target/'verification.json',dict(archive_parts_verified=len(index['parts']),entries_verified=len(index['entries']),measurements_object_exact=True,
              new_solve=False,real_checkpoint_resume_executed_by_replay=False,source='sealed raw evidence in isolated subprocess'))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--mode',choices=['pack','replay'],required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--runs',type=Path);p.add_argument('--source',type=Path);p.add_argument('--retained',type=Path,nargs='*',default=[]);p.add_argument('--archive',type=Path,default=Path(__file__).resolve().parent);a=p.parse_args()
    if a.mode=='pack':pack(a.runs.resolve(),a.source.resolve(),a.output.resolve(),[x.resolve() for x in a.retained])
    else:replay(a.archive.resolve(),a.output.resolve())
if __name__=='__main__':main()
