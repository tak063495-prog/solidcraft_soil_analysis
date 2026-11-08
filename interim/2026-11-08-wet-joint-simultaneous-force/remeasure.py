"""Verify every saved entry, then remeasure without any solver invocation."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools'))
from diagnose_cut_volume_history import write
from diagnose_wet_joint_simultaneous import measure
from wet_joint_cycle_reference import serializable


def remeasure(evidence):
    manifest = json.loads((Path(__file__).parent/'archive-sha256.json').read_text())
    for name, digest in manifest['entries'].items():
        assert hashlib.sha256((evidence/name).read_bytes()).hexdigest() == digest, ('changed entry', name)
    actual = serializable(measure(evidence/'paths'))
    expected = json.loads((evidence/'paths/measurements.json').read_text())
    assert actual == expected, 'saved measurements differ'
    assert not actual['passed'] and not actual['history_converged'] and not actual['production_adopted']
    return dict(schema='soil-fem-wet-joint-simultaneous-remeasurement/1', passed=False,
                history_converged=False, production_adopted=False, solver_rerun=False,
                independent_path_rerun=False, all_entries_sha256_verified=True,
                all_measurements_exactly_reproduced=True,
                exact_volume_equilibrium_accepted=actual['force_alignment']['exact_volume_equilibrium_accepted'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    write(args.output, remeasure(args.evidence.resolve()))
    print('Saved evidence reproduced; passed=false')
