"""Seal only this task's bounded cloud observations; never edit old evidence."""
import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import zipfile
from pcfail_evidence import verify


def seal(source,destination):
    if (destination/'archive-index.json').exists():raise FileExistsError('existing seal retained')
    entries={}
    with tempfile.TemporaryFile() as stream:
        with zipfile.ZipFile(stream,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
            for path in sorted(source.rglob('*')):
                if not path.is_file():continue
                name=path.relative_to(source).as_posix();raw=path.read_bytes()
                entries[name]=dict(size=len(raw),sha256=hashlib.sha256(raw).hexdigest());z.writestr(name,raw)
        stream.seek(0);parts=[];h=hashlib.sha256();total=0
        while raw:=stream.read(16*1024*1024):
            name='evidence.zip.part'+str(len(parts)).zfill(3)
            with (destination/name).open('xb') as f:f.write(raw)
            h.update(raw);total+=len(raw);parts.append(dict(name=name,size=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
    index=dict(schema='pcfail-evidence/1',archive_bytes=total,archive_sha256=h.hexdigest(),parts=parts,entries=entries)
    with (destination/'archive-index.json').open('x',encoding='utf-8') as f:json.dump(index,f,indent=2);f.write('\n')
    return verify(destination,index)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--destination',type=Path,required=True);a=p.parse_args()
    print(json.dumps(seal(a.source,a.destination)))
