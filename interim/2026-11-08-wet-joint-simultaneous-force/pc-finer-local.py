"""User PC only: existing small specimen at 2048/64; never certify acceptance."""
import argparse
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools'))
from diagnose_cut_volume_history import inputs, validate_path, summary, write
from diagnose_wet_joint_simultaneous import recorded
from wet_joint_cycle_reference import cut_volume_reference, differences


def run(output):
    output.mkdir(parents=True, exist_ok=False)
    values = inputs()
    for name, value in zip(('mesh', 'ground', 'stages'), values):
        write(output/(name+'.json'), value)
    start = time.monotonic()
    result = cut_volume_reference(*values, 2048, 64, 7)
    elapsed = time.monotonic()-start
    validate_path(result, *values, 2048, 64)
    write(output/'g2048-l64.json', result)
    # The old mechanical solution is read only after the new input-only solve.
    old, provenance = recorded(Path(__file__).parent, ['paths/g2048-l32.json'])
    baseline = old['paths/g2048-l32.json']
    validate_path(baseline, *values, 2048, 32)
    rows = differences(baseline, result)
    report = dict(schema='soil-fem-wet-joint-pc-local/1', passed=False,
                  history_converged=False, production_adopted=False,
                  elapsed_seconds=elapsed, baseline=provenance, fields=rows, maxima=summary(rows),
                  limitation='Seven full fields and accepted schedule metadata only; no full-history acceptance certificate.')
    write(output/'measurements.json', report)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = run(args.output.resolve())
    print(json.dumps(dict(passed=result['passed'], history_converged=result['history_converged'])))
