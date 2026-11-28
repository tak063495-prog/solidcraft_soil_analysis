"""Small existing wet JOINT consistency evidence; never a convergence pass."""
import argparse
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import audit_wet_joint_accepted_history as replay
from wet_joint_cycle_reference import serializable


def run(out):
    out.mkdir(parents=True, exist_ok=False)
    values = replay.history.inputs()
    for name, value in zip(('mesh', 'ground', 'stages'), values):
        (out/(name+'.json')).write_text(json.dumps(value, indent=2)+'\n')
    checks = {}
    for count in (4, 8):
        folder = out/f'g{count}'
        writer = replay.history.HistoryWriter(folder/'accepted', chunk_size=3)
        started = time.monotonic()
        result = replay.history.cut_volume_reference(*values, count, 8, 7, history_sink=writer)
        writer.close()
        (folder/'result.json').write_text(json.dumps(serializable(result), indent=2)+'\n')
        checks[str(count)] = dict(audit=replay.audit(folder/'accepted', result, values), elapsed_seconds=time.monotonic()-started)
    report = dict(passed=False, history_converged=False, production_adopted=False,
                  checks=checks, common_accepted_comparison=replay.history.compare_accepted(out/'g4/accepted', out/'g8/accepted', 2),
                  native_rerun=False, independent_native_comparison=False,
                  p3b_p4_complete=False, structure_water_joint_complete=False, p6_complete=False)
    (out/'measurements.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    run(p.parse_args().output.resolve())
