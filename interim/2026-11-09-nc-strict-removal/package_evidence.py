"""Package raw evidence without binaries; retain initial postprocessing failures."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import zipfile

def package(root,out):
    if out.exists():raise FileExistsError(out)
    entries={}
    with zipfile.ZipFile(out,'x',compression=zipfile.ZIP_LZMA) as archive:
        for folder,dirs,files in os.walk(root):
            dirs[:]=sorted(d for d in dirs if d not in ('target','.git','__pycache__','solver'))
            for name in sorted(files):
                p=Path(folder)/name; rel=p.relative_to(root).as_posix()
                if name.startswith('nc_stop_probe-') or name.endswith(('.exe','.so','.rlib','.rmeta')):continue
                data=p.read_bytes();entries[rel]=hashlib.sha256(data).hexdigest();archive.writestr(rel,data)
    manifest=dict(schema='soil-fem-nc-evidence-archive/1',archive_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),
                  entry_count=len(entries),entries=entries,excludes=['executables','target/','job/solver/'],passed=False,history_converged=False)
    with (out.parent/'archive-sha256.json').open('x') as f:json.dump(manifest,f,indent=2);f.write('\n')
    return manifest

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--evidence',type=Path,required=True);p.add_argument('--archive',type=Path,required=True)
    a=p.parse_args();r=package(a.evidence.resolve(),a.archive.resolve());print(r['entry_count'],r['archive_sha256'])
