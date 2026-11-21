"""New input-only runs of the existing small cycle; all flags remain incomplete."""
import argparse
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import audit_wet_joint_bulk_history as bulk
from wet_joint_cycle_reference import serializable


def run(out):
    out.mkdir(parents=True, exist_ok=False)
    values = bulk.history.inputs()
    for name, value in zip(('mesh', 'ground', 'stages'), values):
        (out/(name+'.json')).write_text(json.dumps(value, indent=2)+'\n', encoding='utf-8')
    reports = {}
    for steps in (16, 32):
        folder = out/f'g{steps}-l8'; folder.mkdir()
        sink = bulk.history.HistoryWriter(folder/'accepted', chunk_size=16)
        started = time.monotonic()
        result = bulk.history.cut_volume_reference(*values, steps, 8, 7, history_sink=sink)
        sink.close()
        result = serializable(result)
        (folder/'result.json').write_text(json.dumps(result)+'\n', encoding='utf-8')
        report = bulk.audit(folder/'accepted', result, values)
        (folder/'audit.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
        reports[f'g{steps}-l8'] = dict(elapsed_seconds=time.monotonic()-started,
                                    accepted_states=report['accepted_states'],
                                    maxima=report['bulk_water_maximum_absolute_differences'])
    comparison = bulk.history.compare_accepted(out/'g16-l8/accepted', out/'g32-l8/accepted', ratio=2)
    report = dict(passed=False, input_only_reference=True, native_rerun=False,
                  full_history_convergence_accepted=False, independent_bulk_plastic_solution=False,
                  p3b_p4_complete=False, structure_water_joint_complete=False, p6_complete=False,
                  paths=reports, comparison=comparison)
    (out/'measurements.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    run(p.parse_args().output.resolve())
