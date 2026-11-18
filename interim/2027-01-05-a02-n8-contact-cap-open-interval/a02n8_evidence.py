"""Seal only new observations; verify every part/entry byte without a solver."""
import argparse, hashlib, json
from pathlib import Path
import tempfile, zipfile

def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def seal(sources,destination):
    entries={};parts=[]
    with tempfile.TemporaryFile() as stream:
        with zipfile.ZipFile(stream,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
            for prefix,folder in sources:
                for p in sorted(folder.rglob('*')):
                    if not p.is_file() or '__pycache__' in p.parts:continue
                    name=prefix+'/'+p.relative_to(folder).as_posix()
                    if name in entries:raise ValueError('duplicate evidence entry')
                    z.write(p,name);entries[name]=dict(size=p.stat().st_size,sha256=sha(p))
        size=stream.tell();stream.seek(0);digest=hashlib.sha256()
        while data:=stream.read(24*1024*1024):
            name='evidence.zip.part'+str(len(parts)).zfill(3);p=destination/name
            with p.open('xb') as out:out.write(data)
            digest.update(data);parts.append(dict(name=name,size=len(data),sha256=sha(p)))
    index=dict(archive_bytes=size,archive_sha256=digest.hexdigest(),parts=parts,entries=entries,
        raw_failures_retained=True,solver_binaries_archived=False)
    with (destination/'archive-index.json').open('x',encoding='utf-8') as f:json.dump(index,f,indent=2);f.write('\n')
    return index

def verify(folder,extract=None):
    m=json.loads((folder/'archive-index.json').read_text(encoding='utf-8'));digest=hashlib.sha256();size=0
    if extract:extract.mkdir(parents=True,exist_ok=False)
    with tempfile.TemporaryFile() as stream:
        for part in m['parts']:
            p=folder/part['name']
            if p.stat().st_size!=part['size'] or sha(p)!=part['sha256']:raise ValueError('part SHA differs: '+str(p))
            with p.open('rb') as f:
                while data:=f.read(1024*1024):stream.write(data);digest.update(data);size+=len(data)
        if size!=m['archive_bytes'] or digest.hexdigest()!=m['archive_sha256']:raise ValueError('archive SHA differs')
        stream.seek(0)
        with zipfile.ZipFile(stream) as z:
            if len(z.namelist())!=len(set(z.namelist())) or set(z.namelist())!=set(m['entries']):raise ValueError('entry inventory differs')
            for name in z.namelist():
                p=Path(name)
                if p.is_absolute() or '..' in p.parts:raise ValueError('unsafe entry')
                raw=z.read(name);entry=m['entries'][name]
                if len(raw)!=entry['size'] or hashlib.sha256(raw).hexdigest()!=entry['sha256']:raise ValueError('entry SHA differs: '+name)
                if extract:
                    target=extract/p;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
    return dict(parts=len(m['parts']),entries=len(m['entries']),all_sha_verified=True,archive_bytes=size)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--verify',type=Path,required=True);p.add_argument('--extract',type=Path);a=p.parse_args();print(json.dumps(verify(a.verify,a.extract)))
