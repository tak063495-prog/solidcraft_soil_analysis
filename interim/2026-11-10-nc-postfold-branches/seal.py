"""Seal all small-case raw fields/linear systems and the exact audit source."""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile

HERE=Path(__file__).resolve().parent


class Parts:
    def __init__(self,out):
        self.out=out;self.buffer=bytearray();self.parts=[];self.sha=hashlib.sha256();self.size=0
    def write(self,data):
        self.sha.update(data);self.size+=len(data);self.buffer.extend(data)
        while len(self.buffer)>=16*1024*1024:self.emit(16*1024*1024)
        return len(data)
    def emit(self,count):
        data=bytes(self.buffer[:count]);del self.buffer[:count]
        name=f'evidence.tar.gz.part{len(self.parts):03d}'
        (self.out/name).write_bytes(data)
        self.parts.append(dict(name=name,size=len(data),sha256=hashlib.sha256(data).hexdigest()))
    def flush(self):pass
    def finish(self):
        if self.buffer:self.emit(len(self.buffer))


def run(source,prepared,out):
    if (out/'archive-sha256.json').exists():raise ValueError('archive already exists')
    files={}
    for path in source.rglob('*'):
        if path.is_file():files['study/'+path.relative_to(source).as_posix()]=path
    for name in ['preparation.json','observer/preparation.json','prior/extraction.json']:
        files['prepared/'+name]=prepared/name
    for path in HERE.glob('*.py'):files['drivers/'+path.name]=path
    entries={name:dict(size=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest()) for name,path in sorted(files.items())}
    parts=Parts(out)
    with gzip.GzipFile(fileobj=parts,mode='wb',mtime=0,compresslevel=6) as compressed:
        with tarfile.open(fileobj=compressed,mode='w|',format=tarfile.PAX_FORMAT) as archive:
            for name,path in sorted(files.items()):
                info=tarfile.TarInfo(name);info.size=entries[name]['size'];info.mode=0o644;info.mtime=0
                with path.open('rb') as stream:archive.addfile(info,stream)
    parts.finish()
    result=dict(archive_sha256=parts.sha.hexdigest(),archive_size=parts.size,parts=parts.parts,entries=entries,
                numerical_acceptance=False)
    (out/'archive-sha256.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(archive_size=parts.size,parts=len(parts.parts),entries=len(entries))),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['source','prepared','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();run(a.source,a.prepared,a.output)
