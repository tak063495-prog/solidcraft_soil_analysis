"""Prepare six solver settings for one original 96-Tet10 input."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

METHODS = ('repeated_refinement', 'extended_residual', 'ruiz', 'gmres_polish')
LEGACY = Path(__file__).resolve().parents[1] / '2026-11-22-impl-linear-solver' / 'prepare_inputs.py'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--source', type=Path)
    choice.add_argument('--cloud-original', action='store_true')
    args = parser.parse_args()
    if args.source:
        original = json.loads((args.source / 'stages.json').read_bytes())
        if 'dense_residual_method' in original['settings']:
            parser.error('source must be the original baseline stage case')
    source = ['--source', str(args.source)] if args.source else ['--cloud-original']
    subprocess.run([sys.executable, LEGACY, args.output, *source], check=True)
    case = json.loads((args.output / 'stages.json').read_bytes())
    for method in METHODS:
        case['settings']['dense_residual_method'] = method
        path = args.output / f'stages-{method}.json'
        path.write_bytes((json.dumps(case, ensure_ascii=False, indent=2, allow_nan=False) + '\n').encode())


if __name__ == '__main__':
    main()
