"""Verify both raw archives and recompute the saved measurements without solving."""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile
import analyze
import audit

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def remeasure(out):
    out.mkdir(parents=True,exist_ok=False)
    tools=read(HERE/'audit-tool-sha256.json')
    for name,digest in tools.items():
        if sha(ROOT/name)!=digest:raise ValueError('audit dependency changed: '+name)
    manifests=read(HERE/'archive-sha256.json')['archives']
    verified={}
    for name,manifest in manifests.items():
        archive=HERE/name
        if sha(archive)!=manifest['archive_sha256']:raise ValueError('archive SHA256 mismatch: '+name)
        with tarfile.open(archive,'r:xz') as z:
            if set(z.getnames())!=set(manifest['entries']) or len(z.getmembers())!=len(manifest['entries']):raise ValueError('archive inventory mismatch: '+name)
            if not all(m.isfile() for m in z.getmembers()):raise ValueError('non-file archive member')
            for n,digest in manifest['entries'].items():
                p=Path(n)
                if p.is_absolute() or '..' in p.parts:raise ValueError('unsafe evidence path')
                b=z.extractfile(n).read()
                if hashlib.sha256(b).hexdigest()!=digest:raise ValueError('entry SHA256 mismatch: '+n)
                if name=='evidence.tar.xz':
                    dest=out/'evidence'/p;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(b)
        verified[name]=len(manifest['entries'])
    evidence=out/'evidence'
    source=read(evidence/'final-source-sha256.json')
    for n,digest in source.items():
        if sha(evidence/'isolated'/n)!=digest:raise ValueError('saved source differs: '+n)
    prior=read(evidence/'prior-archive-sha256.json')
    for p in (evidence/'prior').rglob('*'):
        if p.is_file() and sha(p)!=prior['entries'][p.relative_to(evidence/'prior').as_posix()]:
            raise ValueError('prior entry provenance differs')
    analyze.measure(evidence)
    audit.measure(evidence,out/'recomputed')
    matches={
        'derivative-measurements.json':read(evidence/'derivative-measurements.json')==read(HERE/'derivative-measurements.json'),
        'measurements.json':read(out/'recomputed/measurements.json')==read(HERE/'measurements.json')}
    result=dict(passed=False,history_converged=False,solver_executed=False,material_map_executed=False,
        archives_verified=verified,source_files_verified=len(source),exact_saved_json_matches=matches,
        all_matches=all(matches.values()))
    (out/'remeasurement.json').write_text(json.dumps(result,indent=2)+'\n')
    if not result['all_matches']:raise ValueError('measurement differs; retain actual output')
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    print(json.dumps(remeasure(a.output.resolve()),indent=2))
