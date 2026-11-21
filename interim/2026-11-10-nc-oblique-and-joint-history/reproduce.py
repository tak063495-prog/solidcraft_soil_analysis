"""Run bounded, unchanged-map NC controls and actual partial-checkpoint restarts."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent


def run(args):
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    paths = args.output/'paths'; paths.mkdir()
    def command(kind, label, waypoint, steps, output, extra=()):
        return [sys.executable, str(HERE/'oblique.py'), '--evidence', str(args.evidence),
                '--binary', str(args.binary), '--output', str(output), '--map', kind,
                '--waypoint', *map(str, waypoint), '--step-size', str(1/steps),
                '--budget-seconds', str(args.budget_seconds), '--invocation-commit', args.invocation_commit,
                *extra]
    def execute(item):
        name, cmd = item; start = time.monotonic()
        with (args.output/(name+'.log')).open('x', encoding='utf-8') as stream:
            result = subprocess.run(cmd, stdout=stream, stderr=subprocess.STDOUT)
        return dict(name=name, command=cmd, exit_code=result.returncode, elapsed_seconds=time.monotonic()-start)
    jobs = []
    for kind in ('native', 'independent'):
        for label, waypoint in (('a80-b20', (.8, .2)), ('a20-b80', (.2, .8))):
            for steps in (128, 256):
                name = f'{kind}-{label}-{steps}'
                jobs.append((name, command(kind, label, waypoint, steps, paths/name)))
    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(execute, jobs))
    (args.output/'processes.json').write_text(json.dumps(results, indent=2)+'\n', encoding='utf-8')
    if any(x['exit_code'] for x in results):
        raise ValueError('diagnostic child failed; do not count as numerical result')
    restarts = []
    import gzip
    for kind in ('native', 'independent'):
        name = f'{kind}-a80-b20-128'; original = paths/name
        checkpoint = original/'midpoint.json'
        if not checkpoint.exists():
            checkpoint = original/'early.json'
        cp = json.loads(checkpoint.read_text(encoding='utf-8'))
        with gzip.open(original/'samples.jsonl.gz', 'rt') as stream:
            samples = [json.loads(line) for line in stream]
        remaining = samples[-1]['step']-cp['payload']['step']
        if remaining < 1:
            raise ValueError('actual checkpoint must precede observed terminal state')
        cmd = command(kind, 'a80-b20', (.8, .2), 128, paths/(name+'-resumed'),
                      ['--resume', str(checkpoint), '--count', str(remaining)])
        process = execute((name+'-resumed', cmd)); results.append(process)
        if process['exit_code']:
            raise ValueError('partial restart process failed')
        with gzip.open(paths/(name+'-resumed')/'samples.jsonl.gz', 'rt') as stream:
            resumed = [json.loads(line) for line in stream]
        expected = [s for s in samples if s['step'] >= cp['payload']['step']]
        # The first saved row is the checkpoint itself. The direct row also
        # records its preceding du; the restart's initial row has no preceding du.
        first = dict(expected[0]); first.pop('reduced_increment_m', None)
        assert resumed == [first, *expected[1:]], 'partial restart fields differ'
        assert (original/'restart.json').read_bytes() == (paths/(name+'-resumed')/'restart.json').read_bytes()
        restarts.append(dict(path=name, checkpoint_step=cp['payload']['step'],
                             checkpoint_coordinate=cp['payload']['coordinate'], compared_fields=len(resumed),
                             all_saved_fields_exact=True, terminal_checkpoint_bytes_exact=True,
                             completion_verified=False))
    (args.output/'processes.json').write_text(json.dumps(results, indent=2)+'\n', encoding='utf-8')
    (args.output/'restarts.json').write_text(json.dumps(restarts, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('evidence', 'binary', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--invocation-commit', required=True)
    p.add_argument('--budget-seconds', type=float, default=90.)
    run(p.parse_args())
