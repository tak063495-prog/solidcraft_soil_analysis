"""Build the hash-verified historical observer and fresh input-only NC prefix."""
import argparse
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
REF = 'f19ce642647207a87b9dda6cc31b006b85b5936c'


def run(out):
    path = ROOT/'interim/2026-11-10-nc-mean-stress-release/prepare.py'
    spec = importlib.util.spec_from_file_location('historical_nc_prepare', path)
    prior = importlib.util.module_from_spec(spec); spec.loader.exec_module(prior)
    original = shutil.copytree
    def frozen_copy(source, destination, *args, **kwargs):
        if Path(source).resolve() != ROOT/'src':
            return original(source, destination, *args, **kwargs)
        destination = Path(destination); destination.mkdir()
        listing = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', REF, '--', 'src'], cwd=ROOT).decode().splitlines()
        if not listing:
            raise ValueError('historical observer source unavailable; do not substitute current source')
        for name in listing:
            target = destination/Path(name).relative_to('src')
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(subprocess.check_output(['git', 'show', REF+':'+name], cwd=ROOT))
        return destination
    prior.shutil.copytree = frozen_copy
    try:
        prior.run(out)
    finally:
        prior.shutil.copytree = original
    record = json.loads((out/'preparation.json').read_text())
    record.update(observer_source_ref=REF, observer_source_verified_against_existing_69_source_manifest=True,
        checkout_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT).decode().strip(),
        historical_observer_is_not_current_production_binary=True)
    (out/'preparation.json').write_text(json.dumps(record, indent=2)+'\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    run(p.parse_args().output.resolve())
