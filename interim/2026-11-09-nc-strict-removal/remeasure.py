"""Verify every archive byte, then recompute all declared measurements without solving."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile
from audit_evidence import measure

HERE=Path(__file__).resolve().parent

def remeasure(archive,out):
    manifest=json.loads((HERE/'archive-sha256.json').read_text())
    if hashlib.sha256(archive.read_bytes()).hexdigest()!=manifest['archive_sha256']:raise ValueError('archive checksum mismatch')
    out.mkdir(parents=True,exist_ok=False);evidence=out/'evidence';evidence.mkdir()
    with zipfile.ZipFile(archive) as z:
        if len(z.namelist())!=len(set(z.namelist())) or set(z.namelist())!=set(manifest['entries']):raise ValueError('entry inventory differs')
        for name,digest in manifest['entries'].items():
            p=Path(name)
            if p.is_absolute() or '..' in p.parts:raise ValueError('unsafe archive path')
            data=z.read(name)
            if hashlib.sha256(data).hexdigest()!=digest:raise ValueError('entry checksum mismatch: '+name)
            dest=evidence/p;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(data)
    measure(evidence,out/'audit')
    matches={}
    for p in sorted((out/'audit').glob('*.json')):
        matches[p.name]=json.loads(p.read_text())==json.loads((evidence/'audit-final'/p.name).read_text())
    result=dict(passed=False,history_converged=False,solver_executed=False,material_map_executed=False,
                archive_entries_verified=len(manifest['entries']),exact_saved_json_matches=matches,all_matches=all(matches.values()))
    (out/'remeasurement.json').write_text(json.dumps(result,indent=2)+'\n')
    if not result['all_matches']:raise ValueError('measurement differs; retain actual output')
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--archive',type=Path,default=HERE/'evidence.zip');p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();print(json.dumps(remeasure(a.archive.resolve(),a.output.resolve()),indent=2))
