"""Package this diagnosis only; production and old evidence remain untouched."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile
HERE=Path(__file__).resolve().parent

def package(work):
    paths=[]
    for source,prefix in [(work/'prior','prior'),(work/'diagnostic','diagnostic'),(work/'final','measurements')]:
        for p in sorted(source.rglob('*')):
            if p.is_file():paths.append((p,str(Path(prefix)/p.relative_to(source)).replace('\\','/')))
    for n in ('src/stages.rs','examples/nc_stop_probe.rs','examples/nc_midpoint_map_probe.rs','original-source-sha256.json'):
        paths.append((work/'source'/n,'captured-source/'+n))
    paths.append((work/'measurements/material-process.json','original-material-process.json'))
    manifest={}
    archive=HERE/'evidence.zip'
    with archive.open('xb') as stream,zipfile.ZipFile(stream,'w',compression=zipfile.ZIP_LZMA) as z:
        for p,n in paths:
            b=p.read_bytes();z.writestr(n,b);manifest[n]=hashlib.sha256(b).hexdigest()
    record=dict(archive_sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),entries=manifest)
    with (HERE/'archive-sha256.json').open('x') as f:json.dump(record,f,indent=2)
    return len(manifest)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--work',type=Path,required=True);a=p.parse_args();print(package(a.work.resolve()))
