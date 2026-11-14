"""Verify every archive byte and regenerate JSON, without running any solver."""
import argparse,hashlib,json,pathlib,shutil,tarfile
from analyze import measure
HERE=pathlib.Path(__file__).resolve().parent

def remeasure(out):
    out.mkdir(parents=True,exist_ok=False);manifest=json.loads((HERE/'archive-sha256.json').read_text());archive=out/'joined.tar.xz';h=hashlib.sha256()
    with archive.open('xb') as target:
        for part in manifest['parts']:
            data=(HERE/part['name']).read_bytes()
            if len(data)!=part['bytes'] or hashlib.sha256(data).hexdigest()!=part['sha256']:raise ValueError('archive part differs')
            h.update(data);target.write(data)
    if h.hexdigest()!=manifest['archive_sha256']:raise ValueError('joined archive differs')
    e=out/'evidence';e.mkdir();count=0
    with tarfile.open(archive,'r:xz') as z:
        members=z.getmembers()
        if len(members)!=len(manifest['entries']) or set(m.name for m in members)!=set(manifest['entries']):raise ValueError('archive inventory differs')
        for m in members:
            p=pathlib.Path(m.name)
            if not m.isfile() or p.is_absolute() or '..' in p.parts:raise ValueError('unsafe archive entry')
            data=z.extractfile(m).read()
            if hashlib.sha256(data).hexdigest()!=manifest['entries'][m.name]:raise ValueError('entry differs: '+m.name)
            dest=e/p;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(data);count+=1
    expected=json.loads((e/'diagnostic-source-sha256.json').read_text())
    for name,digest in expected.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('diagnostic source differs: '+name)
    result=measure(e,out/'recomputed');saved=json.loads((HERE/'measurements.json').read_text());exact=result==saved
    record=dict(passed=False,history_converged=False,solver_executed=False,material_map_executed=False,archive_entries_verified=count,exact_saved_json_match=exact)
    (out/'remeasurement.json').write_text(json.dumps(record,indent=2)+'\n')
    if not exact:raise ValueError('remeasured JSON differs; retain both versions')
    return record
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=pathlib.Path,required=True);a=p.parse_args();print(json.dumps(remeasure(a.output.resolve()),indent=2))
