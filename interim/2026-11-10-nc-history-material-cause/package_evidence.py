"""Package exact small-case evidence and sealed sources, without solver binaries."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def run(data):
    archive=data/'new-evidence.zip';entries={}
    def add(z,name,raw):
        if name in entries:raise ValueError('duplicate package member')
        z.writestr(name,raw);entries[name]=dict(size=len(raw),sha256=hashlib.sha256(raw).hexdigest())
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for folder in ('input128','input256','trace128','trace256','direct128','direct256','restarts','prior'):
            for p in sorted((data/folder).rglob('*')):
                if p.is_file():add(z,'data/'+p.relative_to(data).as_posix(),p.read_bytes())
        for p in sorted(data.glob('*.log')):add(z,'process-logs/'+p.name,p.read_bytes())
        for name in ('direct-processes.json',):add(z,'data/'+name,(data/name).read_bytes())
        add(z,'expected-material-history.json',(data/'measurements-causal.json').read_bytes())
        for name in ('source-snapshot.zip','source-snapshot.json'):
            add(z,name,(data/'validation'/name).read_bytes())
        for p in sorted(HERE.glob('*.py')):
            add(z,'source/interim/'+HERE.name+'/'+p.name,p.read_bytes())
    digest=hashlib.sha256();size=0;parts=[]
    with archive.open('rb') as stream:
        for i,raw in enumerate(iter(lambda:stream.read(2*1024*1024),b'')):
            name=f'evidence.zip.part{i:03d}'
            with (HERE/name).open('xb') as out:out.write(raw)
            digest.update(raw);size+=len(raw);parts.append(dict(name=name,size=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
    result=dict(archive_sha256=digest.hexdigest(),archive_size=size,parts=parts,entries=entries,
                solver_binaries_included=False,numerical_acceptance=False)
    with (HERE/'archive-sha256.json').open('x',encoding='utf-8') as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(dict(bytes=size,parts=len(parts),entries=len(entries)),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--data',type=Path,required=True)
    run(p.parse_args().data.resolve())
