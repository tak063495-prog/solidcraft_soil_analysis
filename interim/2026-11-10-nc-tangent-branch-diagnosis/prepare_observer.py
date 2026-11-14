"""Build an isolated, read-only material observer from the starting checkout."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ADDITIONS = ROOT/'interim/2026-11-09-nc-turning-path-diagnosis'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run(out):
    out.mkdir(parents=True, exist_ok=False)
    source = out/'source'
    source.mkdir()
    shutil.copytree(ROOT/'src', source/'src')
    originals = {p.relative_to(ROOT).as_posix(): sha(p) for p in (ROOT/'src').rglob('*.rs')}
    for name in ('Cargo.toml', 'Cargo.lock', 'rust-toolchain.toml'):
        shutil.copyfile(ROOT/name, source/name)
    for name, addition in [('mcc.rs', 'mcc-diagnostic-addition.rs'), ('mcc_midpoint.rs', 'midpoint-diagnostic-addition.rs')]:
        path = source/'src'/name
        path.write_bytes(path.read_bytes().replace(b'\r\n', b'\n')+(ADDITIONS/addition).read_bytes())
    (source/'examples').mkdir()
    shutil.copyfile(ADDITIONS/'material_diagnostic.rs', source/'examples/material_diagnostic.rs')
    subprocess.run(['cargo', 'fmt', '--manifest-path', str(source/'Cargo.toml')], check=True)
    subprocess.run(['cargo', 'build', '--locked', '--offline', '--release', '--manifest-path', str(source/'Cargo.toml'),
                    '--target-dir', str(out/'target'), '--example', 'material_diagnostic'], check=True)
    binary = out/'target/release/examples'/('material_diagnostic.exe' if __import__('sys').platform == 'win32' else 'material_diagnostic')
    record = dict(starting_commit='113d83b8a37b3b68ce46eb00b6664f042334eae0', source_original_sha256=originals,
                  observer_source_sha256={p.relative_to(source).as_posix(): sha(p) for p in source.rglob('*') if p.is_file()},
                  binary_sha256=sha(binary), production_modified=False, numerical_acceptance=False)
    (out/'preparation.json').write_text(json.dumps(record, indent=2)+'\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    run(p.parse_args().output.resolve())
