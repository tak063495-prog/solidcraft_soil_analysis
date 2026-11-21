"""Recompute the unchanged native initial/load/unload prefix at 5 kPa."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools'))
from probe_nc_midpoint import inputs


def run(solver, evidence, output):
    output.mkdir(parents=True, exist_ok=False)
    mesh, ground, case = inputs('midpoint')
    original = json.loads(json.dumps(case))
    case['stages'] = case['stages'][:2]
    names = {s['name'] for s in case['stages']}
    case['verification_reject_once_stages'] = [n for n in case['verification_reject_once_stages'] if n in names]
    for name, value in zip(('mesh', 'ground', 'stages'), (mesh, ground, case)):
        (output/(name+'.json')).write_text(json.dumps(value, indent=2)+'\n')
    command = [str(solver), 'solve-stages', *[str(output/(n+'.json')) for n in ('mesh', 'ground', 'stages')], str(output/'result.json')]
    with (output/'solve.log').open('x') as log:
        child = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    if child.returncode:
        raise subprocess.CalledProcessError(child.returncode, command)
    current = json.loads((output/'result.json').read_text())
    before = json.loads((evidence/'prior/fresh-head-midpoint/job/recovery-checkpoint-0.json').read_text())['run']
    a = next(f['state'] for f in current['frames'] if f['name'] == 'unload')
    b = next(f['state'] for f in before['frames'] if f['name'] == 'unload')
    result = dict(passed=False, exit_code=child.returncode, prefix_complete=current['complete'], fresh_unload_entire_state_exact=a == b,
                  original_case_sha256=hashlib.sha256(json.dumps(original, sort_keys=True).encode()).hexdigest(),
                  solver_sha256=hashlib.sha256(solver.read_bytes()).hexdigest(), command=command,
                  scope='Same original first two stages and their verification injections. Excavation and its named verification injection omitted only from this prefix run. No new loading history.')
    (output/'verification.json').write_text(json.dumps(result, indent=2)+'\n')
    if not current['complete'] or a != b:
        raise ValueError('fresh native prefix differs from archived starting history')
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('solver', 'evidence', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    print(json.dumps(run(a.solver.resolve(), a.evidence.resolve(), a.output.resolve()), indent=2))
