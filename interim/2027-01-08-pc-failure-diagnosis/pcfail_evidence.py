"""Sparse-aware evidence verification; missing bytes are incomplete, exit 2."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
import tempfile
import zipfile
ROOT=Path(__file__).resolve().parents[2]
PROFILES={
    'b01':('2026-12-25-b01-curved-tri6-cause-diagnosis','archive-index.json'),
    'b01-cli':('2026-12-03-tri6-curved-b01','evidence-manifest.json'),
    'a01':('2026-12-20-a01-joint-combo-cause-diagnosis','archive-index.json'),
    'a01-old':('2026-11-23-joint-combo-extension','evidence-manifest.json'),
}


def sha(path):
    with path.open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()


def safe(name):
    p=PurePosixPath(name)
    if p.is_absolute() or '..' in p.parts or '\\' in name or ':' in name:
        raise ValueError('unsafe evidence path: '+name)
    return p


def manifest(profile):
    directory,name=PROFILES[profile];folder=ROOT/'interim'/directory
    return folder,json.loads((folder/name).read_text(encoding='utf-8'))


def parts(folder,index):
    result=[]
    for item in index['parts']:
        name=item.get('name',item.get('path'));safe(name)
        if len(PurePosixPath(name).parts)!=1: raise ValueError('part must be a basename')
        result.append((folder/name,item.get('size',item.get('bytes')),item['sha256']))
    return result


def availability(folder,index):
    rows=[]
    for path,size,digest in parts(folder,index):
        state='missing' if not path.is_file() else ('verified' if path.stat().st_size==size and sha(path)==digest else 'corrupt')
        rows.append(dict(path=path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else path.as_posix(),
                         status=state, expected_bytes=size, expected_sha256=digest))
    status='failed' if any(r['status']=='corrupt' for r in rows) else ('incomplete' if any(r['status']=='missing' for r in rows) else 'available')
    return dict(status=status, passed=False, numerical_acceptance=False, parts=rows)


def verify(folder,index,extract=None,profile=None):
    result=availability(folder,index)
    if result['status']!='available': return result
    members=index.get('entries')
    if members is None: members={r['path']:dict(size=r['bytes'],sha256=r['sha256']) for r in index['members']}
    # Never create a misleading unsealed directory when required parts are missing.
    if extract is not None and extract.exists(): raise FileExistsError(extract)
    digest=hashlib.sha256();size=0
    with tempfile.TemporaryFile() as stream:
        for path,_,_ in parts(folder,index):
            with path.open('rb') as f:
                for raw in iter(lambda:f.read(1024*1024),b''):
                    stream.write(raw);digest.update(raw);size+=len(raw)
        if size!=index['archive_bytes'] or digest.hexdigest()!=index['archive_sha256']: raise ValueError('combined archive SHA/size mismatch')
        stream.seek(0)
        with zipfile.ZipFile(stream) as z:
            names=z.namelist()
            if len(names)!=len(set(names)) or set(names)!=set(members): raise ValueError('member inventory mismatch')
            # Verify everything before extraction. Two passes keep memory bounded.
            for name in names:
                safe(name);h=hashlib.sha256();n=0
                with z.open(name) as f:
                    for raw in iter(lambda:f.read(1024*1024),b''): h.update(raw);n+=len(raw)
                if n!=members[name]['size'] or h.hexdigest()!=members[name]['sha256']: raise ValueError('member SHA/size mismatch: '+name)
            if extract is not None:
                extract.mkdir(parents=True,exist_ok=False)
                for name in names:
                    if profile=='a01-old' and (not name.startswith(('selected/','interacting-stacked-attempt/')) or 'restarts/' in name or 'checkpoint-' in name): continue
                    target=extract.joinpath(*safe(name).parts);target.parent.mkdir(parents=True,exist_ok=True)
                    with z.open(name) as source,target.open('xb') as dest:
                        for raw in iter(lambda:source.read(1024*1024),b''):dest.write(raw)
    result.update(status='verified',all_bytes_verified=True,entries=len(members),archive_bytes=size,archive_sha256=digest.hexdigest())
    return result


def fetch(folder,index,commit):
    if not re.fullmatch('[0-9a-f]{40}',commit): raise ValueError('fixed full commit SHA required')
    # Only listed parts, no checkout, history traversal or archive expansion.
    result=availability(folder,index)
    if result['status']=='failed': raise ValueError('existing corrupt evidence retained; use a separate worktree')
    for path,size,digest in parts(folder,index):
        if path.exists(): continue
        rel=path.relative_to(ROOT).as_posix()
        with tempfile.TemporaryFile() as stream:
            subprocess.run(['git','show',commit+':'+rel],cwd=ROOT,stdout=stream,check=True)
            stream.seek(0);h=hashlib.sha256();count=0
            for raw in iter(lambda:stream.read(1024*1024),b''):h.update(raw);count+=len(raw)
            if count!=size or h.hexdigest()!=digest:raise ValueError('fixed blob differs from manifest: '+rel)
            stream.seek(0)
            path.parent.mkdir(parents=True,exist_ok=True)
            with path.open('xb') as dest:
                for raw in iter(lambda:stream.read(1024*1024),b''):dest.write(raw)
    return availability(folder,index)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--profile',choices=PROFILES,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--extract',type=Path)
    p.add_argument('--fetch-commit');a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    try:
        folder,index=manifest(a.profile)
        if a.fetch_commit:fetch(folder,index,a.fetch_commit)
        result=verify(folder,index,a.extract,a.profile)
        code=0 if result['status']=='verified' else (2 if result['status']=='incomplete' else 1)
    except (OSError,ValueError,KeyError,subprocess.CalledProcessError,zipfile.BadZipFile) as exc:
        result=dict(status='failed',passed=False,numerical_acceptance=False,error=repr(exc));code=1
    (a.output/'evidence-status.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False));return code

if __name__=='__main__':raise SystemExit(main())
