"""Explicit, bounded observation only. Mesh sizes 12/96, load steps 4/8/16."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from b02cause_binding import ROOT, OLD, legacy

HERE = Path(__file__).resolve().parent


def write(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def command(cmd, log, expected=0):
    start = time.monotonic()
    with log.open('x', encoding='utf-8') as f:
        p = subprocess.Popen(list(map(str, cmd)), cwd=ROOT, stdout=f, stderr=subprocess.STDOUT)
        try:
            code = p.wait(timeout=7200)
            timed_out = False
        except subprocess.TimeoutExpired:
            p.terminate()
            try:
                code = p.wait(timeout=20)
            except subprocess.TimeoutExpired:
                p.kill()
                code = p.wait()
            timed_out = True
    result = dict(command=list(map(str, cmd)), pid=p.pid, exit_code=code,
                  elapsed_seconds=time.monotonic()-start, timed_out=timed_out, log_sha256=sha(log))
    write(log.with_suffix('.process.json'), result)
    if timed_out or (expected is not None and code != expected):
        raise RuntimeError('retained failed/incomplete process: '+str(log))
    return result


def reverse_keys(value):
    if isinstance(value, dict):
        return {k: reverse_keys(value[k]) for k in reversed(value)}
    if isinstance(value, list):
        return [reverse_keys(v) for v in value]
    return value


def sources():
    names = subprocess.check_output(['git', 'ls-files'], cwd=ROOT, text=True).splitlines()
    paths = {ROOT/n for n in names if Path(n).suffix in ('.rs', '.py', '.ps1', '.toml') or Path(n).name == 'Cargo.lock'}
    paths.update(p for p in HERE.iterdir() if p.suffix in ('.py', '.md'))
    return {p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(paths)}


def run(output, observer, cli, steps=(4, 8, 16)):
    output.mkdir(parents=True, exist_ok=False)
    (output/'bin').mkdir()
    for name, source in [('observer', observer), ('cli', cli)]:
        target = output/'bin'/(name+source.suffix)
        shutil.copy2(source, target)
    binary, cli = [output/'bin'/(name+source.suffix) for name, source in [('observer', observer), ('cli', cli)]]
    before = sources()
    write(output/'source-binding-before.json', before)
    write(output/'executable-binding.json', {p.name: sha(p) for p in (binary, cli)})
    records, restarts = [], []
    for level in (0, 1):
        for n in steps:
            name = f'tet{12*8**level}-n{n}'
            folder = output/name
            folder.mkdir()
            m, g, c, refinement = legacy['cases'].inputs(level, n, .05, True)
            for key, value in [('mesh', m), ('ground', g), ('stages', c), ('refinement', refinement)]:
                write(folder/(key+'.json'), value)
            base = [folder/(s+'.json') for s in ('mesh', 'ground', 'stages')]
            stops = [n+n//2, 2*n+n//2, 3*n+n//2, 4*n+n//2] if (level, n) == (0, 4) else [2*n+n//2]
            native = command([binary, 'solve', *base, folder/'native', ','.join(map(str, stops))], folder/'native.log')
            independent = command([sys.executable, HERE/'b02cause_reference.py', folder, folder/'independent.jsonl'], folder/'independent.log')
            comparison = legacy['compare'].compare(folder)
            plain = command([binary, 'plain', *base, folder/'plain'], folder/'plain.log')
            original = legacy['cases'].read(folder/'native/run.json')
            assert original['complete'] and legacy['study'].exact(original, legacy['cases'].read(folder/'plain/run.json'))
            observations = [json.loads(s) for s in (folder/'native/accepted.jsonl').read_text().splitlines()]
            for stop in stops:
                target = folder/f'resumed-{stop}'
                process = command([binary, 'resume', *base, target, folder/f'native/checkpoint-{stop}.json'], folder/f'resumed-{stop}.log')
                later = [json.loads(s) for s in (target/'accepted.jsonl').read_text().splitlines()]
                expected = [s for s in observations if s['accepted_index'] > stop]
                assert process['pid'] != native['pid']
                assert legacy['study'].exact(legacy['cases'].read(target/'run.json'), original)
                assert legacy['study'].exact(later, expected)
                restarts.append(dict(case=name, stop=stop, later_states=len(later), exact_except_elapsed=True, process=process))
            records.append(dict(case=name, level=level, load_steps=n, nodes=len(m['nodes']), tet10=len(m['tetrahedra']),
                                accepted_states=len(observations), native=native, independent=independent, plain=plain, comparison=comparison))
            write(output/'plan-results.json', dict(passed=False, b02_complete=False, diagnostic_only=True, cases=records, restarts=restarts))
            print(name, 'states', len(observations), 'old gates', comparison['existing_field_gates_passed'], flush=True)
    if 4 in steps:
        for size in (12, 96):
            folder = output/f'tet{size}-n4'
            for label, options in [('gl7', ['--order', '7']), ('reverse-q', ['--reverse'])]:
                command([sys.executable, HERE/'b02cause_reference.py', folder, folder/(label+'.jsonl'), *options], folder/(label+'.log'))
        original_folder = output/'tet12-n4'
        repeat = output/'repeat-tet12-n4'
        repeat.mkdir()
        keys = output/'keys-reversed-tet12-n4'
        keys.mkdir()
        for folder, reverse in [(repeat, False), (keys, True)]:
            for key in ('mesh', 'ground', 'stages'):
                obj = legacy['cases'].read(original_folder/(key+'.json'))
                write(folder/(key+'.json'), reverse_keys(obj) if reverse else obj)
            base = [folder/(s+'.json') for s in ('mesh', 'ground', 'stages')]
            command([binary, 'plain', *base, folder/'plain'], folder/'plain.log')
            command([sys.executable, HERE/'b02cause_reference.py', folder, folder/'independent.jsonl'], folder/'independent.log')
            assert legacy['study'].exact(legacy['cases'].read(folder/'plain/run.json'), legacy['cases'].read(original_folder/'native/run.json'))
            assert (folder/'independent.jsonl').read_bytes() == (original_folder/'independent.jsonl').read_bytes()
        # Reuse unchanged scope and authenticated replay probes in child PIDs.
        command([sys.executable, OLD/'scope.py', output/'scope', '--cli', cli], output/'scope.log')
        mutation_root = output/'mutation-source'
        mutation_root.mkdir()
        (mutation_root/'solver').mkdir()
        shutil.copy2(binary, mutation_root/'solver'/binary.name)
        shutil.copytree(original_folder, mutation_root/'beam-positive')
        command([sys.executable, OLD/'mutations.py', mutation_root, output/'mutations'], output/'mutations.log')
    after = sources()
    write(output/'source-binding-after.json', after)
    assert before == after
    return records


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--observer', type=Path, required=True)
    p.add_argument('--cli', type=Path, required=True)
    p.add_argument('--steps', type=int, choices=(4, 8, 16), nargs='+', default=[4, 8, 16])
    a = p.parse_args()
    if len(set(a.steps)) != len(a.steps):
        p.error('duplicate step count')
    run(a.output.resolve(), a.observer.resolve(), a.cli.resolve(), tuple(sorted(a.steps)))
