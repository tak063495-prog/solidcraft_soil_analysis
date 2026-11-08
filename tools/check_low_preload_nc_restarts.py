"""Compare separate-process partial restarts in all four existing NC phases."""
import argparse
from pathlib import Path
import subprocess

import probe_low_preload_nc as low


def check(solver, output):
    output.mkdir(parents=True, exist_ok=False)
    rows = []
    for phase, count, stage_index in (('geostatic', 4, 0), ('load', 9, 1),
                                      ('unload', 13, 2), ('excavate_upper', 15, 3)):
        folder = output/phase; folder.mkdir()
        low.probe(solver, folder/'first', 60., 5., stop_after_accepted=count)
        first = low.read(folder/'first/observed-run.json')
        progress = first['geostatic']['progress'] if stage_index == 0 else first['last_committed']['stage_progress']
        actual_index = first['last_committed']['stage_index']
        if actual_index != stage_index or not 0 < progress < 1:
            raise ValueError('observation did not stop inside phase: '+phase)
        low.probe(solver, folder/'resumed', 60., 5., resume_from=folder/'first/job', stop_after_accepted=1)
        resumed = low.read(folder/'resumed/observed-run.json')
        final_count = low.accepted_count(resumed)
        if final_count <= low.accepted_count(first):
            raise ValueError('no additional accepted increment: '+phase)
        inputs = [folder/'first/input'/(n+'.json') for n in low.NAMES]
        with (folder/'direct.log').open('x', encoding='utf-8') as log:
            result = subprocess.run([str(solver), 'checkpoint-stages', *map(str, inputs),
                                     str(final_count), str(folder/'direct.json')],
                                    stdout=log, stderr=subprocess.STDOUT, timeout=90, check=False)
        if result.returncode:
            raise ValueError('direct prefix failed: '+phase)
        direct = low.read(folder/'direct.json')['run']
        low.compare_prefix(resumed, direct)
        if low.deterministic(resumed) != low.deterministic(direct):
            raise ValueError('full run (including continuation/input/layout) differs: '+phase)
        rows.append(dict(phase=phase, stage_index=stage_index, starting_progress=progress,
                         accepted_before=low.accepted_count(first), accepted_after=final_count,
                         exact_entire_run_except_elapsed_seconds=True))
    report = dict(schema='soil-fem-low-preload-all-phase-partial-restarts/1', passed=False,
                  partial_restart_comparisons_exact=True, comparisons=rows,
                  solver_sha256=low.sha(solver), excluded_recursive_fields=['elapsed_seconds'],
                  limitations=['各工程途中から次の受理点までの一致。全撤去完了への再開一致は未確認'])
    low.write(output/'verification.json', report)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--solver', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(check(args.solver.resolve(), args.output.resolve()))
