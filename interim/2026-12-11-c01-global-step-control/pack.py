"""Pack newly solved raw fields with per-entry and per-part SHA-256 (never executables)."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def pack(source,destination,baseline=None,compat=None):
    index=destination/'archive-index.json'
    if index.exists(): raise FileExistsError(index)
    archive=source.parent/'c01-evidence.zip'
    files=[]
    for folder,prefix in [(source,'paths'),(baseline,'legacy-baseline'),(compat,'legacy-current')]:
        if folder is None: continue
        for path in sorted(folder.rglob('*')):
            if not path.is_file() or path.name in ('probe','probe.exe'):continue
            files.append((path,prefix+'/'+path.relative_to(folder).as_posix()))
    entries=[]
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
        for path,name in files:
            entries.append(dict(path=name,bytes=path.stat().st_size,sha256=sha(path)))
            z.write(path,name)
    parts=[]
    with archive.open('rb') as f:
        for i,chunk in enumerate(iter(lambda:f.read(48*1024*1024),b'')):
            path=destination/f'evidence.zip.part{i:03}'
            with path.open('xb') as out:out.write(chunk)
            parts.append(dict(file=path.name,bytes=len(chunk),sha256=hashlib.sha256(chunk).hexdigest()))
    index.write_text(json.dumps(dict(schema='c01-control-archive/1',bytes=archive.stat().st_size,sha256=sha(archive),parts=parts,entries=entries),ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(len(entries),'entries',len(parts),'parts',archive.stat().st_size,flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--destination',type=Path,required=True);p.add_argument('--baseline',type=Path);p.add_argument('--compat',type=Path);a=p.parse_args()
    pack(a.source.resolve(),a.destination.resolve(),a.baseline,a.compat)
if __name__=='__main__': main()
