"""Partition Cargo's discovered tests without dropping or ignoring a target.

Each integration target belongs to exactly one sorted round-robin shard. Core
runs lib/bin tests and doctests, and retains the default example build check.
No numerical settings or assertions are changed by this runner.
"""
import argparse
import json
from pathlib import Path
import re
import subprocess
import sys
import time

GROUPS = ('rust-core', 'rust-integration-0', 'rust-integration-1')


def partition(targets):
    names = sorted(t['name'] for t in targets if 'test' in t['kind'])
    if not names or len(names) != len(set(names)):
        raise ValueError('expected a nonempty unique integration-test inventory')
    return {f'rust-integration-{i}': names[i::2] for i in range(2)}


def commands(group, shards):
    if group == 'rust-core':
        return [['cargo', 'test', '--locked', '--release', '--lib', '--bins'],
                ['cargo', 'test', '--locked', '--release', '--doc'],
                ['cargo', 'build', '--locked', '--release', '--examples']]
    names = shards[group]
    if not names:
        raise ValueError('empty integration shard')
    return [['cargo', 'test', '--locked', '--release'] +
            [arg for name in names for arg in ('--test', name)]]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--group', choices=GROUPS, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--inventory-only', action='store_true')
    args = p.parse_args()
    root = Path(__file__).resolve().parents[1]
    meta = json.loads(subprocess.check_output(
        ['cargo', 'metadata', '--locked', '--no-deps', '--format-version', '1'], cwd=root))
    package, = [pkg for pkg in meta['packages'] if Path(pkg['manifest_path']).resolve() == root/'Cargo.toml']
    shards = partition(package['targets'])
    report = dict(group=args.group, integration_target_shards=shards,
                  integration_targets_total=sum(map(len, shards.values())),
                  commands=commands(args.group, shards), executed=not args.inventory_only,
                  exit_code=None, test_results=[], elapsed_seconds=None)
    args.output.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    status = 0
    try:
        if not args.inventory_only:
            with (args.output/'rust-tests.log').open('w', encoding='utf-8') as log:
                for command in report['commands']:
                    print(' '.join(command), flush=True)
                    child = subprocess.Popen(command, cwd=root, stdout=subprocess.PIPE,
                                             stderr=subprocess.STDOUT, text=True, encoding='utf-8', errors='replace')
                    for line in child.stdout:
                        print(line, end='', flush=True)
                        log.write(line)
                        match = re.search(r'test result: (\w+)\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out', line)
                        if match:
                            report['test_results'].append(dict(zip(
                                ('result', 'passed', 'failed', 'ignored', 'measured', 'filtered_out'),
                                [match[1], *map(int, match.groups()[1:])])))
                    status = child.wait()
                    if status:
                        break
    except BaseException:
        status = 1
        raise
    finally:
        report['exit_code'] = status
        report['elapsed_seconds'] = time.monotonic()-started
        (args.output/'rust-test-inventory.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    return status


if __name__ == '__main__':
    sys.exit(main())
