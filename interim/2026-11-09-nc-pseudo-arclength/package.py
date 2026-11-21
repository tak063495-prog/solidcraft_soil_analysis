"""Seal raw small-case evidence, without binaries; split below GitHub file limits."""
import argparse,hashlib,json,pathlib,tarfile
HERE=pathlib.Path(__file__).resolve().parent

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def package(w):
    archive=w.parent/(w.name+'-evidence.tar.xz');assert not archive.exists()
    files=[p for p in sorted(w.rglob('*')) if p.is_file() and not any(n in ('target','__pycache__','source-reproduction') for n in p.relative_to(w).parts)]
    entries={p.relative_to(w).as_posix():sha(p) for p in files}
    with tarfile.open(archive,'w:xz',preset=3) as z:
        for p in files:z.add(p,arcname=p.relative_to(w).as_posix(),recursive=False)
    parts=[]
    with archive.open('rb') as src:
        index=1
        while True:
            b=src.read(40*1024*1024)
            if not b:break
            dest=HERE/f'evidence.tar.xz.part{index:02d}'
            with dest.open('xb') as f:f.write(b)
            parts.append(dict(name=dest.name,bytes=len(b),sha256=hashlib.sha256(b).hexdigest()));index+=1
    result=dict(archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,parts=parts,entries=entries,entry_count=len(entries),uncompressed_bytes=sum(p.stat().st_size for p in files),binaries_included=False)
    with (HERE/'archive-sha256.json').open('x') as f:json.dump(result,f,indent=2)
    print(json.dumps({k:v for k,v in result.items() if k!='entries'},indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--evidence',type=pathlib.Path,required=True);a=p.parse_args();package(a.evidence.resolve())
