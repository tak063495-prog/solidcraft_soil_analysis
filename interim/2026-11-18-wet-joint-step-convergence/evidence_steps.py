"""Lossless small-case evidence: SHA-bound extraction and full remeasurement."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import zipfile
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('c01_measure_steps',HERE/'measure_steps.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
read=m.read;write=m.write;sha=m.st.sha;ROOT=m.st.ROOT
REPORTS=('measurements.json','switch-events.json','all-entity-differences.json','chord-defect.json')


def validate_report(r):
    for k in ('passed','physical_model_changed','tolerances_changed','production_solver_changed','c02_c03_c04_accepted','p6_complete'):
        if r[k] is not False:raise ValueError('unsupported acceptance or model change: '+k)
    if r['c01_decision']!='未達' or r['global_steps']!=[16,32,64,128,256]:raise ValueError('changed scope/decision')
    for p in r['pairs'].values():
        if p['full_history_acceptance'] is not False:raise ValueError('unsupported full-history acceptance')
    return r


def pack(paths,measurements,archive):
    files={}
    for p in sorted(paths.rglob('*')):
        if p.is_file() and p.name not in ('probe','probe.exe'):
            files['paths/'+p.relative_to(paths).as_posix()]=p
    for n in REPORTS:files[n]=measurements/n
    source_paths=[]
    for folder in ('src','tools','examples','tests','interim/2026-11-15-wet-joint-convergence','interim/2026-11-18-wet-joint-step-convergence','interim/2026-11-08-wet-joint-finer-full-history'):
        source_paths += [p for p in (ROOT/folder).rglob('*') if p.is_file() and p.suffix in ('.rs','.py')]
    source_paths += [ROOT/n for n in ('Cargo.toml','Cargo.lock','rust-toolchain.toml','requirements.lock.txt','requirements-reference.lock.txt')]
    for p in sorted(source_paths):files['sources/'+p.relative_to(ROOT).as_posix()]=p
    entries={}
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for name,p in files.items():
            entries[name]=sha(p)
            z.write(p,name,compress_type=zipfile.ZIP_STORED if p.suffix=='.npz' else zipfile.ZIP_DEFLATED)
    parts=[]
    with archive.open('rb') as f:
        while data:=f.read(48*1024*1024):
            p=HERE/f'evidence.zip.part{len(parts)+1:02d}';p.write_bytes(data)
            parts.append(dict(name=p.name,bytes=len(data),sha256=sha(p)))
    write(HERE/'archive-sha256.json',dict(archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,parts=parts,entries=entries,solver_binary_archived=False))


def extract(destination):
    manifest=read(HERE/'archive-sha256.json');destination.mkdir(parents=True,exist_ok=False)
    archive=destination/'evidence.zip'
    with archive.open('xb') as f:
        for part in manifest['parts']:
            p=HERE/part['name'];assert p.stat().st_size==part['bytes'] and sha(p)==part['sha256'],'changed archive part'
            with p.open('rb') as src:shutil.copyfileobj(src,f)
    assert sha(archive)==manifest['archive_sha256'],'changed archive'
    with zipfile.ZipFile(archive) as z:
        assert len(z.namelist())==len(set(z.namelist())) and set(z.namelist())==set(manifest['entries'])
        for name,digest in manifest['entries'].items():
            path=Path(name);assert not path.is_absolute() and '..' not in path.parts
            data=z.read(name);assert hashlib.sha256(data).hexdigest()==digest,('changed entry',name)
            p=destination/path;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)
    archive.unlink();validate_report(read(destination/'measurements.json'));return manifest


def replay(destination):
    manifest=extract(destination);actual=m.measure(destination/'paths')
    for name,r in zip(REPORTS,actual):assert r==read(destination/name),('remeasured mismatch',name)
    return dict(all_entries_sha256_verified=len(manifest['entries']),all_four_measurement_objects_exact=True,numerical_acceptance=False,solver_rerun=False,software_replay_passed=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path);p.add_argument('--pack-paths',type=Path);p.add_argument('--measurements',type=Path);p.add_argument('--archive',type=Path);a=p.parse_args()
    if a.pack_paths:pack(a.pack_paths.resolve(),a.measurements.resolve(),a.archive.resolve())
    else:
        result=replay(a.output.resolve());write(a.output/'replay-verification.json',result);print(json.dumps(result))
