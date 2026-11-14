"""Explicit manual NC diagnostics; never invoked by default CI tests."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PREVIOUS = ROOT/'interim/2026-11-10-nc-mean-stress-release'


def commands(out):
    prepared = out/'prepared'
    binary = prepared/'observer-target/release/examples'/('material_diagnostic.exe' if sys.platform == 'win32' else 'material_diagnostic')
    common = ['--evidence', str(prepared/'origin'), '--binary', str(binary)]
    rows = []
    for kind, policy, steps, material, reference, order, suffix in [
        ('native', 'plain', 128, 4, 400, 'complement', ''),
        ('native', 'refined', 128, 4, 400, 'complement', ''),
        ('native', 'refined', 256, 4, 400, 'complement', ''),
        ('independent', 'plain', 128, 4, 400, 'complement', ''),
        ('independent', 'refined', 128, 4, 400, 'complement', ''),
        ('independent', 'refined', 256, 4, 400, 'complement', ''),
        ('independent', 'refined', 128, 8, 400, 'complement', '-m8'),
        ('independent', 'refined', 128, 4, 800, 'complement', '-pre800'),
        ('native', 'refined', 128, 4, 400, 'mean', ''),
        ('independent', 'refined', 128, 4, 400, 'mean', ''),
    ]:
        label = f'{kind}-{order}-{policy}{steps}{suffix}'
        command = [sys.executable, str(HERE/'continuation.py'), '--map', kind, '--linear-policy', policy,
            '--order', order+'-first', '--step-size', str(1/steps), '--material-steps', str(material),
            '--minimum-step', str(2**-28), '--budget-seconds', '900', '--count', '4096',
            '--reference', str(prepared/f'reference{reference}'), '--output', str(out/label), *common]
        rows.append((label, command))
    return rows


def run(out, jobs, reuse_prepared=False):
    if reuse_prepared:
        if not (out/'prepared/preparation.json').is_file() or set(p.name for p in out.iterdir()) != {'prepared'}:
            raise ValueError('reuse requires only a completed preparation in output')
    else:
        out.mkdir(parents=True, exist_ok=False)
    records = []
    def invoke(label, command, timeout=None):
        start = time.monotonic(); expired = False
        with (out/(label+'.log')).open('x') as log:
            child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
            try:
                code = child.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                expired = True; child.terminate(); child.wait(); code = child.returncode
        record = dict(label=label, command=command, exit_code=code, external_timeout=expired,
            elapsed_seconds=time.monotonic()-start)
        print(json.dumps(record), flush=True)
        return record
    preparation_commands = [
        ('prepare', [sys.executable, str(HERE/'prepare.py'), '--output', str(out/'prepared')]),
        ('reference800', [sys.executable, str(ROOT/'interim/2026-11-09-nc-pseudo-arclength/capture_reference_history.py'),
            '--output', str(out/'prepared/reference800'), '--steps', '800']),
        ('fresh-native-prefix', [sys.executable, str(PREVIOUS/'fresh_prefix.py'),
            '--solver', str(ROOT/'target/release'/('soil-fem-3d.exe' if sys.platform == 'win32' else 'soil-fem-3d')),
            '--evidence', str(out/'prepared/origin'), '--output', str(out/'fresh-native-prefix')]),
    ]
    if reuse_prepared:
        records.append(dict(label='prepare', reused_completed_preparation=json.loads((out/'prepared/preparation.json').read_text()), exit_code=0))
        preparation_commands = preparation_commands[1:]
    for label, command in preparation_commands:
        record = invoke(label, command); records.append(record)
        (out/'processes.json').write_text(json.dumps(records, indent=2)+'\n')
        if record['exit_code']:
            raise ValueError('preparation failed; never substitute native history into independent path')
    specs = commands(out)
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = [pool.submit(invoke, label, command, 930) for label, command in specs]
        for future in futures:
            records.append(future.result())
            (out/'processes.json').write_text(json.dumps(records, indent=2)+'\n')
    if any(row['exit_code'] for row in records):
        raise ValueError('diagnostic process incomplete; inspect logs and preserve timeout as inconclusive')
    restart_specs = []
    for label, command in specs:
        if '-refined' not in label or '-m8' in label or '-pre800' in label:
            continue
        checkpoint = out/label/'tail.json'
        if not checkpoint.exists():
            continue
        resumed = label+'-resumed'
        command = list(command)
        command[command.index('--output')+1] = str(out/resumed)
        command += ['--resume', str(checkpoint)]
        restart_specs.append((resumed, command))
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = [pool.submit(invoke, label, command, 930) for label, command in restart_specs]
        for future in futures:
            records.append(future.result())
            (out/'processes.json').write_text(json.dumps(records, indent=2)+'\n')
    if any(row['exit_code'] for row in records):
        raise ValueError('restart process incomplete; do not claim completion')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--jobs', type=int, choices=[1, 2, 3, 4], default=4)
    p.add_argument('--list-only', action='store_true')
    p.add_argument('--reuse-prepared', action='store_true')
    args = p.parse_args()
    if args.list_only:
        print(json.dumps(commands(args.output.resolve()), indent=2))
    else:
        run(args.output.resolve(), args.jobs, args.reuse_prepared)
