"""Continue each unfinished new path from its actual sealed terminal checkpoint."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent


def run(study, evidence, binary, budget):
    jobs = []
    for folder in sorted((study/'paths').glob('*')):
        if folder.name.endswith(('-resumed', '-extended')):
            continue
        report = json.loads((folder/'report.json').read_text(encoding='utf-8'))
        if report['full_removal_reached']:
            continue
        controls = json.loads((folder/'controls.json').read_text(encoding='utf-8'))
        name = folder.name+'-extended'
        cmd = [sys.executable, str(HERE/'oblique.py'), '--evidence', str(evidence), '--binary', str(binary),
               '--output', str(study/'paths'/name), '--map', report['map'],
               '--waypoint', *map(str, controls['withdrawal_waypoint_at_half']),
               '--step-size', str(controls['binding']['initial_step']), '--resume', str(folder/'restart.json'),
               '--count', '2048', '--budget-seconds', str(budget), '--invocation-commit', controls['invocation_commit']]
        jobs.append((name, cmd))
    def execute(item):
        name, cmd = item; started = time.monotonic()
        with (study/(name+'.log')).open('x', encoding='utf-8') as stream:
            child = subprocess.run(cmd, stdout=stream, stderr=subprocess.STDOUT)
        return dict(name=name, command=cmd, exit_code=child.returncode, elapsed_seconds=time.monotonic()-started)
    with ThreadPoolExecutor(max_workers=2) as pool:
        result = list(pool.map(execute, jobs))
    (study/'extension-processes.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    if any(p['exit_code'] for p in result):
        raise ValueError('extension process failed; do not call it a numerical verdict')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('study', 'evidence', 'binary'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--budget-seconds', type=float, default=120.)
    a = p.parse_args()
    run(a.study, a.evidence, a.binary, a.budget_seconds)
