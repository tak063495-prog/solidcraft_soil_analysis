"""Package only this task's new raw evidence; never edit protected evidence."""
import argparse,hashlib,json,zipfile
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def package(work):
    archive=HERE/'evidence.zip';assert not archive.exists()
    entries={};files=[]
    for p in sorted(work.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:
            files.append((p,str(Path('raw')/p.relative_to(work)).replace('\\','/')))
    files.append((ROOT/'target/release/soil-fem-3d','producer/soil-fem-3d-linux'))
    for n in ('src/lib.rs','src/stages.rs','src/hydrostatic_volume.rs','tools/verify_high_order.py','tools/verify_geostatic.py','tools/verify_joint_water.py','tools/hydrostatic_volume_audit.py','examples/hydrostatic_volume_probe.rs','tests/hydrostatic_volume.rs','tests/test_hydrostatic_volume_audit.py'):
        files.append((ROOT/n,'captured-source/'+n))
    for p in sorted(HERE.glob('*.py')):
        files.append((p,'captured-study/'+p.name))
    sources=[p for base in ('src','examples','tools','tests') for p in (ROOT/base).rglob('*') if p.is_file() and p.suffix in ('.rs','.py')]
    sources+=list(HERE.glob('*.py'))
    provenance={str(p.relative_to(ROOT)):digest(p) for p in sorted(sources)}
    (HERE/'source-sha256.json').write_text(json.dumps(provenance,indent=2)+'\n')
    with archive.open('xb') as stream,zipfile.ZipFile(stream,'w',compression=zipfile.ZIP_LZMA) as z:
        for p,n in files:
            z.write(p,n);entries[n]=digest(p)
    report=dict(archive_sha256=digest(archive),bytes=archive.stat().st_size,entries=entries)
    (HERE/'archive-sha256.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(entries=len(entries),bytes=report['bytes'],archive_sha256=report['archive_sha256'])))
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work',type=Path,required=True);a=p.parse_args();package(a.work.resolve())
