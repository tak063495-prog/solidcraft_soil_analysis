"""Verify the small archive, then remeasure with its sealed UTF-8 sources."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import zipfile

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def extract(output):
    output.mkdir(parents=True, exist_ok=False)
    manifest = read(HERE/'archive-sha256.json')
    archive = HERE/'evidence.zip'
    if hashlib.sha256(archive.read_bytes()).hexdigest() != manifest['archive_sha256']:
        raise ValueError('archive changed')
    with zipfile.ZipFile(archive) as z:
        if len(z.namelist()) != len(manifest['entries']) or set(z.namelist()) != set(manifest['entries']):
            raise ValueError('archive inventory changed')
        for name, expected in manifest['entries'].items():
            path = Path(name)
            if path.is_absolute() or '..' in path.parts or '\\' in name:
                raise ValueError('unsafe archive path')
            raw = z.read(name)
            if len(raw) != expected['bytes'] or hashlib.sha256(raw).hexdigest() != expected['sha256']:
                raise ValueError('archive entry changed: '+name)
            target = output/path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(raw)
    return manifest


def run(output):
    manifest = extract(output)
    sources = output/'sources'
    with zipfile.ZipFile(output/'source-snapshot.zip') as z:
        expected = read(output/'source-snapshot.json')
        if set(z.namelist()) != set(expected) or len(z.namelist()) != len(expected):
            raise ValueError('source inventory changed')
        for name, digest in expected.items():
            path = Path(name)
            if path.is_absolute() or '..' in path.parts or '\\' in name:
                raise ValueError('unsafe source path')
            raw = z.read(name)
            if hashlib.sha256(raw).hexdigest() != digest:
                raise ValueError('source changed: '+name)
            target = sources/path; target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(raw)
    script = sources/'interim/2026-11-10-wet-joint-small-refinement/study.py'
    command = [sys.executable, '-I', '-X', 'utf8', str(script), '--replay', str(output/'paths'),
               '--output', str(output/'sealed-remeasurement.json')]
    subprocess.run(command, check=True)
    report = dict(archive_entries_verified=len(manifest['entries']), sealed_sources_verified=len(expected),
                  isolated_utf8_command=command, **read(output/'sealed-remeasurement.json'))
    (output/'replay-verification.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    print(json.dumps(run(p.parse_args().output.resolve())))
