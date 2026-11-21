"""Run the bounded spatial-order diagnostics on the existing small NC case."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent


def commands(w):
    prepared = w/'prepared'
    binary = prepared/'observer-target/release/examples'/('material_diagnostic.exe' if sys.platform == 'win32' else 'material_diagnostic')
    common = ['--evidence', str(prepared/'origin'), '--binary', str(binary), '--reference', str(prepared/'reference400')]
    result = []
    for kind in ['native', 'independent']:
        for order in ['left', 'center']:
            for steps in [128, 256]:
                label = f'{kind}-{order}{steps}'
                result.append((label, [sys.executable, str(HERE/'study.py'), '--map', kind, '--order', order+'-first',
                    *common, '--output', str(w/label), '--step-size', str(1/steps), '--budget-seconds', str(600 if steps == 128 else 900)]))
    for label, script, extras in [
        ('independent-left128-m8', 'study.py', ['--material-steps', '8']),
        ('independent-event-left128', 'event_study.py', [])]:
        result.append((label, [sys.executable, str(HERE/script), '--map', 'independent', '--order', 'left-first',
            *common, '--output', str(w/label), '--step-size', str(1/128), '--budget-seconds', '600', *extras]))
    for kind, order, seed in [('native', 'left', w/'native-left128/midpoint.json'),
                              ('independent', 'center', w/'restart-seeds/independent-center128.json')]:
        label = f'{kind}-{order}-resumed'
        result.append((label, [sys.executable, str(HERE/'study.py'), '--map', kind, '--order', order+'-first',
            *common, '--output', str(w/label), '--step-size', str(1/128), '--resume', str(seed),
            '--count', '12', '--budget-seconds', '300']))
    for kind, order in [('native', 'left'), ('independent', 'center')]:
        label=f'{kind}-{order}-tail'
        result.append((label, [sys.executable, str(HERE/'study.py'), '--map', kind, '--order', order+'-first',
            *common, '--output', str(w/label), '--step-size', str(1/128), '--resume', str(w/f'{kind}-{order}128/restart.json'),
            '--count', '1024', '--budget-seconds', '600']))
    return result


def capture(w):
    """Copy an actual live checkpoint; never reconstruct/reseal its payload."""
    p = w/'independent-center128/restart.json'; q = w/'restart-seeds/independent-center128.json'
    if q.exists() or not p.exists():
        return
    try:
        raw = p.read_bytes(); saved = json.loads(raw)
        checksum = hashlib.sha256(json.dumps(saved['payload'], sort_keys=True, allow_nan=False).encode()).hexdigest()
        if saved['payload_sha256'] == checksum and .25 <= saved['payload']['coordinate'] < .5:
            q.parent.mkdir(exist_ok=True); q.write_bytes(raw)
    except (ValueError, KeyError):
        pass  # Writer may currently be replacing restart.json; retry next poll.


def run(w):
    w.mkdir(parents=True, exist_ok=False)
    subprocess.run([sys.executable, str(HERE.parent/'2026-11-10-nc-component-release/prepare.py'), '--output', str(w/'prepared')], check=True)
    records = []
    for label, cmd in commands(w):
        start = time.monotonic()
        with (w/(label+'.log')).open('x') as log:
            child = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)
            while child.poll() is None:
                if label == 'independent-center128':
                    capture(w)
                time.sleep(.1)
        records.append(dict(label=label, command=cmd, exit_code=child.returncode, elapsed_seconds=time.monotonic()-start))
        (w/'numerical-processes.json').write_text(json.dumps(records, indent=2)+'\n')
        if child.returncode:
            raise subprocess.CalledProcessError(child.returncode, cmd)
    subprocess.run([sys.executable, str(HERE/'measure.py'), '--evidence', str(w), '--output', str(w/'measurements.json')], check=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True); p.add_argument('--list-only', action='store_true')
    a = p.parse_args()
    if a.list_only:
        print(json.dumps(commands(a.output.resolve()), indent=2))
    else:
        run(a.output.resolve())
