"""Losslessly recompress completed JSONL streams; verify every raw byte."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def compact(evidence):
    records=[]
    for folder in sorted(evidence.iterdir()):
        if not (folder/'report.json').exists():
            continue
        for path in sorted(folder.glob('*.jsonl.gz')):
            old_size=path.stat().st_size
            with path.open('rb') as stream:old_sha=hashlib.file_digest(stream,'sha256').hexdigest()
            temporary=Path(str(path)+'.compact');checksum=hashlib.sha256();size=0
            with gzip.open(path,'rb') as source,temporary.open('xb') as raw,gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0,compresslevel=6) as destination:
                for chunk in iter(lambda:source.read(1024*1024),b''):
                    checksum.update(chunk);size+=len(chunk);destination.write(chunk)
            with gzip.open(temporary,'rb') as stream:
                verified=hashlib.file_digest(stream,'sha256').hexdigest()
            if verified!=checksum.hexdigest():raise ValueError('recompression changed raw bytes')
            temporary.replace(path)
            with path.open('rb') as stream:new_sha=hashlib.file_digest(stream,'sha256').hexdigest()
            records.append(dict(path=path.relative_to(evidence).as_posix(),raw_sha256=verified,raw_size=size,
                old_compressed_size=old_size,new_compressed_size=path.stat().st_size,
                old_compressed_sha256=old_sha,new_compressed_sha256=new_sha,raw_bytes_exact=True))
            print(path.relative_to(evidence),old_size,path.stat().st_size,flush=True)
    result=dict(all_raw_bytes_exact=True,streams=records,numerical_values_modified=False)
    (evidence/'lossless-recompression.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--evidence',type=Path,required=True);a=p.parse_args()
    compact(a.evidence.resolve())
