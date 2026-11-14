"""Seal C-04 raw evidence; isolated remeasurement is never a fresh solve."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import zipfile

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('c04_archive_study', HERE/'study.py')
s = importlib.util.module_from_spec(spec); spec.loader.exec_module(s)


def pack(runs, source, output, retained):
    files = []; entries = {}
    for p in sorted(runs.rglob('*')):
        if p.is_file():
            relative = p.relative_to(runs).as_posix()
            name = 'reports/'+p.name if relative in ('measurements.json', 'completion.json') else 'data/final/'+relative
            files.append((p, name))
    for root in retained:
        for p in sorted(root.rglob('*')):
            if p.is_file() and '__pycache__' not in p.parts:
                files.append((p, 'retained/'+root.name+'/'+p.relative_to(root).as_posix()))
    for directory in ('src', 'tools', 'tests', 'examples'):
        for p in sorted((source/directory).rglob('*')):
            if p.is_file() and '__pycache__' not in p.parts and p.suffix in ('.rs', '.py', '.json', '.csv', '.txt'):
                files.append((p, 'source/'+p.relative_to(source).as_posix()))
    for p in sorted((source/'interim').glob('**/*.py')):
        if '__pycache__' not in p.parts: files.append((p, 'source/'+p.relative_to(source).as_posix()))
    for name in ('Cargo.toml', 'Cargo.lock', 'rust-toolchain.toml', 'requirements-reference.lock.txt', 'requirements.lock.txt'):
        files.append((source/name, 'source/'+name))
    assert len(files) == len({name for _, name in files})
    output.mkdir(exist_ok=True, parents=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as tmp:
        archive = Path(tmp)/'evidence.zip'
        with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=6, allowZip64=True) as z:
            for p, name in files:
                before = s.sha(p); z.write(p, name); assert s.sha(p) == before, p
                entries[name] = dict(sha256=before, size=p.stat().st_size)
        parts = []
        with archive.open('rb') as stream:
            for index in range((archive.stat().st_size+19_999_999)//20_000_000):
                target = output/f'evidence.zip.part{index:03}'
                with target.open('xb') as f: f.write(stream.read(20_000_000))
                parts.append(dict(name=target.name, sha256=s.sha(target), size=target.stat().st_size))
        s.write(output/'archive-sha256.json', dict(zip_sha256=s.sha(archive), zip_size=archive.stat().st_size, parts=parts, entries=entries,
                raw_fresh_solves_and_real_checkpoints_retained=True, replay_is_not_new_solve=True))


def replay(archive_folder, output):
    output.mkdir(parents=True, exist_ok=False)
    manifest = s.read(archive_folder/'archive-sha256.json'); archive = output/'evidence.zip'
    with archive.open('xb') as joined:
        for item in manifest['parts']:
            p = archive_folder/item['name']
            assert p.stat().st_size == item['size'] and s.sha(p) == item['sha256']
            with p.open('rb') as f: shutil.copyfileobj(f, joined)
    assert archive.stat().st_size == manifest['zip_size'] and s.sha(archive) == manifest['zip_sha256']
    root = output/'unpacked'; root.mkdir()
    with zipfile.ZipFile(archive) as z:
        assert len(z.namelist()) == len(manifest['entries']) and set(z.namelist()) == set(manifest['entries'])
        for name, info in manifest['entries'].items():
            relative = Path(name); assert not relative.is_absolute() and '..' not in relative.parts
            path = root/relative; path.parent.mkdir(parents=True, exist_ok=True)
            with z.open(name) as src, path.open('xb') as dst: shutil.copyfileobj(src, dst)
            assert path.stat().st_size == info['size'] and s.sha(path) == info['sha256']
    script = root/'source/interim/2026-11-22-wet-joint-transient-water/measure.py'
    result = output/'remeasured.json'; original = s.read(root/'reports/measurements.json')
    bootstrap = "import pathlib,runpy,sys;p=sys.argv[1];sys.path.insert(0,str(pathlib.Path(p).parent));sys.argv=sys.argv[1:];runpy.run_path(p,run_name='__main__')"
    command = [sys.executable, '-I', '-X', 'utf8', '-c', bootstrap, str(script), '--output', str(root/'data/final'), '--report', str(result),
               '--water-steps', *map(str, original['hydrostatic_water_steps']), '--hydraulic-steps', *map(str, original['Richards_physical_steps_per_second'])]
    proc = s.invoke(command, output/'replay.log', 900)
    assert proc['exit_code'] == 0 and s.read(result) == original
    diagnostic = script.with_name('diagnose.py'); diagnostic_result = output/'cause-diagnostics.json'
    proc2 = s.invoke([sys.executable, '-I', '-X', 'utf8', '-c', bootstrap, diagnostic, '--output', root/'data/final', '--report', diagnostic_result], output/'diagnostic.log', 900)
    assert proc2['exit_code'] == 0 and s.read(diagnostic_result) == s.read(archive_folder/'cause-diagnostics.json')
    s.write(output/'verification.json', dict(complete=True, archive_parts_verified=len(manifest['parts']), entries_verified=len(manifest['entries']),
            measurements_object_exact=True, cause_diagnostics_object_exact=True, isolated_sealed_source=True,
            new_solve=False, real_checkpoint_resume_executed_by_replay=False, processes=[proc, proc2]))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--mode', choices=['pack', 'replay'], required=True)
    p.add_argument('--output', type=Path, required=True); p.add_argument('--runs', type=Path); p.add_argument('--source', type=Path)
    p.add_argument('--retained', nargs='*', type=Path, default=[]); p.add_argument('--archive', type=Path, default=HERE)
    a = p.parse_args()
    if a.mode == 'pack': pack(a.runs.resolve(), a.source.resolve(), a.output.resolve(), [v.resolve() for v in a.retained])
    else: replay(a.archive.resolve(), a.output.resolve())
