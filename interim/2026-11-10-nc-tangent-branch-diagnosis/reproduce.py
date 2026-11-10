"""Explicit small-case diagnostics; not part of default numerical acceptance."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def run(evidence, observer, out, solver):
    out.mkdir(parents=True, exist_ok=False)
    binary = observer/'target/release/examples'/('material_diagnostic.exe' if sys.platform == 'win32' else 'material_diagnostic')
    records = []
    def invoke(label, command):
        started = time.monotonic()
        with (out/(label+'.log')).open('x') as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        record = dict(label=label, command=command, exit_code=result.returncode, elapsed_seconds=time.monotonic()-started)
        print(json.dumps(record), flush=True)
        return record
    for label, command in [
        ('lambda-zero', [sys.executable, str(HERE/'lambda_zero.py'), '--evidence', str(evidence),
                         '--observer', str(observer), '--output', str(out/'lambda-zero')]),
        ('fresh-prefix', [sys.executable, str(ROOT/'interim/2026-11-10-nc-mean-stress-release/fresh_prefix.py'),
                          '--solver', str(solver), '--evidence', str(evidence/'prepared/origin'), '--output', str(out/'fresh-prefix')]),
        ('tangent', [sys.executable, str(HERE/'diagnose.py'), '--evidence', str(evidence),
                     '--binary', str(binary), '--output', str(out/'tangent')]),
    ]:
        record = invoke(label, command)
        records.append(record)
        (out/'processes.json').write_text(json.dumps(records, indent=2)+'\n')
        if record['exit_code']:
            raise ValueError('diagnostic preparation failed; preserve log')
    specs = []
    for kind, size, count, steps in [('native', 1e-4, 128, 4), ('native', 5e-5, 256, 4),
                                      ('independent', 1e-4, 128, 4), ('independent', 5e-5, 256, 4),
                                      ('independent', 1e-4, 128, 8)]:
        label = f'arc-{kind}-{size}-m{steps}'
        command = [sys.executable, str(HERE/'arc.py'), '--evidence', str(evidence), '--binary', str(binary),
                   '--output', str(out/label), '--map', kind, '--step-size', str(size), '--count', str(count),
                   '--material-steps', str(steps), '--budget-seconds', '300']
        specs.append((label, command))
    with ThreadPoolExecutor(max_workers=2) as pool:
        for record in pool.map(lambda spec: invoke(*spec), specs):
            records.append(record)
            (out/'processes.json').write_text(json.dumps(records, indent=2)+'\n')
    if any(x['exit_code'] for x in records):
        raise ValueError('diagnostic child failed; no numerical acceptance')
    for label, command in specs:
        if '5e-05' not in label or not (out/label/'middle.json').exists():
            continue
        command = list(command)
        command[command.index('--output')+1] = str(out/(label+'-resumed'))
        command += ['--resume', str(out/label/'middle.json')]
        records.append(invoke(label+'-resumed', command))
        (out/'processes.json').write_text(json.dumps(records, indent=2)+'\n')
    if any(x['exit_code'] for x in records):
        raise ValueError('diagnostic restart child failed; no completion claim')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('evidence', 'observer', 'output', 'solver'):
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    run(a.evidence.resolve(), a.observer.resolve(), a.output.resolve(), a.solver.resolve())
