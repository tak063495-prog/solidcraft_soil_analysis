"""Package new bounded evidence and source; never alter earlier evidence."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def package(work):
    paths=[]
    for name in ('evidence','baseline-trial','candidate-trial','candidate-accepted-control','measurements'):
        paths.extend((p,str(p.relative_to(work)).replace('\\','/')) for p in sorted((work/name).rglob('*')) if p.is_file())
    for name in ('matrices.json','linear-results.json','linear-audit.json','candidate-matrices.json','candidate-linear-results.json','candidate-linear-audit.json',
                 'linear-process.json','candidate-linear-process.json','candidate-trial-process.json','control-processes.json'):
        paths.append((work/name,name))
    for variant in ('candidate','baseline'):
        for n in ('src/stages.rs','src/lib.rs','src/nc_linear_refinement.rs','examples/nc_stop_probe.rs','examples/linear_probe.rs','original-source-sha256.json'):
            paths.append((work/variant/n,'captured-source/'+variant+'/'+n))
    for name in ('start-pr.html','start-commit.html','start-remote-refs.txt'):
        paths.append((ROOT/'work'/name,'start/'+name))
    manifest={}
    with (HERE/'evidence.zip').open('xb') as stream,zipfile.ZipFile(stream,'w',compression=zipfile.ZIP_LZMA) as z:
        for path,name in paths:
            data=path.read_bytes();z.writestr(name,data);manifest[name]=hashlib.sha256(data).hexdigest()
    result=dict(archive_sha256=hashlib.sha256((HERE/'evidence.zip').read_bytes()).hexdigest(),entries=manifest)
    with (HERE/'archive-sha256.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    return len(manifest)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work',type=Path,required=True);a=p.parse_args();print(package(a.work.resolve()))
