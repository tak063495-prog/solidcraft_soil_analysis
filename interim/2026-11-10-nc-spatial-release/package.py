"""Archive accepted fields, real linear systems and source snapshots."""
import argparse
import hashlib
import json
from pathlib import Path
import tarfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def package(evidence):
    paths = [p for p in evidence.iterdir() if p.is_dir() and (p/'controls.json').exists()]
    if not paths or any(not (p/'report.json').exists() for p in paths):
        raise ValueError('a numerical process is missing its final report')
    sources = [p for p in ROOT.glob('tools/*.py') if p.name in {
        'spatial_release_control.py', 'mcc_rate_event_floor.py', 'mcc_rate_rk4.py', 'mcc_rate_fem.py',
        'low_preload_nc_reference.py', 'probe_low_preload_nc.py', 'check_low_preload_nc_restarts.py',
        'joint_reference.py', 'probe_nc_midpoint.py', 'make_free_contact.py'}]
    sources += list(HERE.glob('*.py'))
    sources += list((HERE.parent/'2026-11-09-nc-pseudo-arclength').glob('*.py'))
    sources += [ROOT/'tests/low_preload_cut.rs', ROOT/'tests/test_spatial_release_control.py', ROOT/'tests/test_mcc_event_floor.py', ROOT/'tests/test_nc_spatial_evidence.py']
    snapshots = evidence/'sources'
    for p in sources:
        q = snapshots/p.relative_to(ROOT); q.parent.mkdir(parents=True, exist_ok=True); q.write_bytes(p.read_bytes())
    files = []
    for p in evidence.rglob('*'):
        if not p.is_file():
            continue
        name = p.relative_to(evidence).as_posix()
        if '__pycache__' in p.parts or 'observer-target' in p.parts or 'solver' in p.parts or 'cli' in p.parts or p.suffix == '.tar':
            continue
        if 'material' in p.parts and p.name != 'io-sha256.jsonl':
            continue
        if p.name == 'evidence.tar.gz':
            continue
        if p.parent == evidence and p.name.startswith('partial-'):
            continue
        if p.parent == evidence and p.suffix == '.log' and not (evidence/p.stem/'controls.json').exists() and p.name != 'prepare.log':
            continue
        files.append(p)
    files.sort()
    entries = {p.relative_to(evidence).as_posix(): dict(sha256=sha(p), size=p.stat().st_size) for p in files}
    archive = evidence/'evidence.tar.gz'
    with tarfile.open(archive, 'w:gz', compresslevel=6) as z:
        for p in files:
            z.add(p, arcname=p.relative_to(evidence).as_posix(), recursive=False)
    parts = []
    with archive.open('rb') as f:
        index = 0
        while b := f.read(16*1024*1024):
            p = HERE/f'evidence.tar.gz.part{index:03d}'; p.write_bytes(b)
            parts.append(dict(name=p.name, sha256=sha(p), size=len(b))); index += 1
    manifest = dict(archive_sha256=sha(archive), archive_size=archive.stat().st_size, parts=parts, entries=entries,
        transient_native_material_io='request/result bytes removed after primal/replay checks; SHA/size manifests retained; accepted full fields and actual solved systems included',
        binaries_included=False)
    (HERE/'archive-sha256.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return dict(parts=len(parts), entries=len(entries), archive_size=archive.stat().st_size)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--evidence', type=Path, required=True)
    print(json.dumps(package(p.parse_args().evidence.resolve()), indent=2))
