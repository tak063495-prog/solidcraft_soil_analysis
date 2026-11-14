"""Seal new raw evidence and retained failures; integrity verification is not solving."""
import argparse
import hashlib
from pathlib import Path
import shutil
import zipfile
from study import read,write


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def verify(manifest_path,scratch):
    m=read(manifest_path);scratch.mkdir(parents=True,exist_ok=False);archive=scratch/'joined-evidence.zip'
    with archive.open('xb') as joined:
        for part in m['parts']:
            p=manifest_path.parent/part['name']
            assert p.stat().st_size==part['size'] and sha(p)==part['sha256']
            with p.open('rb') as source:shutil.copyfileobj(source,joined)
    assert archive.stat().st_size==m['archive_bytes'] and sha(archive)==m['archive_sha256']
    with zipfile.ZipFile(archive) as z:
        assert len(z.namelist())==len(set(z.namelist())) and set(z.namelist())==set(m['entries'])
        for name,expected in m['entries'].items():
            h=hashlib.sha256();size=0
            with z.open(name) as f:
                for block in iter(lambda:f.read(1024*1024),b''):size+=len(block);h.update(block)
            assert size==expected['size'] and h.hexdigest()==expected['sha256'],name
    result=dict(archive_verified=True,entries_verified=len(m['entries']),integrity_only_not_fresh_solve=True)
    write(scratch/'verification.json',result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path);p.add_argument('--validation',type=Path)
    p.add_argument('--retained',type=Path,action='append',default=[]);p.add_argument('--destination',type=Path)
    p.add_argument('--scratch',type=Path,required=True);p.add_argument('--verify',type=Path);a=p.parse_args()
    if a.verify:print(verify(a.verify,a.scratch));return
    assert a.run and a.destination
    a.destination.mkdir(parents=True,exist_ok=True);a.scratch.mkdir(parents=True,exist_ok=False)
    folders={'fresh':a.run}
    if a.validation:folders['validation']=a.validation
    for folder in a.retained:folders['retained-failures/'+folder.name]=folder
    entries={};archive=a.scratch/'evidence.zip'
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for prefix,folder in sorted(folders.items()):
            for path in sorted(folder.rglob('*')):
                if path.is_file() and 'bin' not in path.relative_to(folder).parts and '__pycache__' not in path.parts:
                    name=prefix+'/'+path.relative_to(folder).as_posix()
                    assert name not in entries
                    z.write(path,name);entries[name]=dict(size=path.stat().st_size,sha256=sha(path))
    parts=[]
    with archive.open('rb') as f:
        while block:=f.read(24*1024*1024):
            part=a.destination/f'evidence.zip.part{len(parts):03}'
            with part.open('xb') as stream:stream.write(block)
            parts.append(dict(name=part.name,size=len(block),sha256=sha(part)))
    manifest=a.destination/'evidence-manifest.json'
    write(manifest,dict(archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,parts=parts,entries=entries,
        solver_binaries_archived=False,all_selected_raw_states_inputs_checkpoints_and_failures_retained=True))
    print(verify(manifest,a.scratch/'verification'))
    print('sealed',len(entries),'entries',archive.stat().st_size,'bytes',len(parts),'parts')


if __name__=='__main__':main()
