"""Fresh production/public-API and input-only independent solves, bounded cloud plan."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import json
import shutil
import sys
import zipfile
from pathlib import Path
from b01_cases import HERE, B02, ROOT, make, read, write
from study import command, sha, exact
from compare_b01 import compare


def run(out, binary, attached, roof, steps):
    out.mkdir(parents=True, exist_ok=False)
    frozen = out/'solver'/binary.name
    frozen.parent.mkdir()
    shutil.copy2(binary, frozen)
    sources = [p for base in ('src', 'tools', 'tests', 'examples', str(HERE.relative_to(ROOT)), str(B02.relative_to(ROOT)))
               for p in (ROOT/base).rglob('*') if p.suffix in ('.rs', '.py', '.md') or (base == 'examples' and p.suffix == '.json')]
    sources += [ROOT/p for p in ('Cargo.toml', 'Cargo.lock', 'rust-toolchain.toml')]
    before = {str(p.relative_to(ROOT)): sha(p) for p in sources}
    write(out/'source-binding-before.json', before)
    with zipfile.ZipFile(out/'fixed-source-at-start.zip', 'x', compression=zipfile.ZIP_DEFLATED) as z:
        for relative in sorted(before):
            z.write(ROOT/relative, relative)
    records, restarts = [], []
    for kind, n in [('attached', n) for n in attached]+[('roof', n) for n in roof]:
        name = f'{kind}-n{n}-s{steps}'
        folder = out/name
        folder.mkdir()
        for key, value in zip(('mesh', 'ground', 'stages', 'geometry'), make(kind, n, steps)):
            write(folder/(key+'.json'), value)
        base = [folder/(k+'.json') for k in ('mesh', 'ground', 'stages')]
        stops = [2*steps+steps//2]
        if kind == 'attached' and n == min(attached):
            stops = [i*steps+steps//2 for i in range(5)]
        native = command([frozen, 'solve', *base, folder/'native', ','.join(map(str, stops))], folder/'native.log', expect=None)
        independent = command([sys.executable, B02/'reference.py', folder, folder/'independent.jsonl'], folder/'independent.log', expect=None)
        rec = dict(case=name, kind=kind, n=n, load_steps=steps,
                   native_process=native, independent_process=independent,
                   solver_sha256=sha(frozen), input_sha256={p.name: sha(p) for p in base})
        if native['exit_code'] == independent['exit_code'] == 0:
            rec['comparison'] = compare(folder)
            plain = command([frozen, 'plain', *base, folder/'plain'], folder/'plain.log')
            original = read(folder/'native/run.json')
            rec['plain_process'] = plain
            rec['all_run_equal_to_unobserved'] = original['complete'] and exact(original, read(folder/'plain/run.json'))
            assert rec['all_run_equal_to_unobserved']
            obs = [json.loads(s) for s in (folder/'native/accepted.jsonl').read_text().splitlines()]
            for stop in stops:
                checkpoint = folder/f'native/checkpoint-{stop}.json'
                target = folder/f'resumed-{stop}'
                resumed = command([frozen, 'resume', *base, target, checkpoint], folder/f'resumed-{stop}.log')
                later = [json.loads(s) for s in (target/'accepted.jsonl').read_text().splitlines()]
                row = dict(case=name, checkpoint_accepted_index=stop, producer_pid=native['pid'], resume_pid=resumed['pid'],
                           checkpoint_sha256=sha(checkpoint), later_accepted_states=len(later),
                           run_exact_except_elapsed=exact(original, read(target/'run.json')),
                           all_later_states_exact=exact(later, [s for s in obs if s['accepted_index'] > stop]))
                row['passed'] = row['run_exact_except_elapsed'] and row['all_later_states_exact'] and native['pid'] != resumed['pid']
                restarts.append(row)
                write(out/'restarts.json', restarts)
                assert row['passed']
        records.append(rec)
        write(out/'measurements.json', dict(passed=False, b01_complete=False, cases=records))
        print(name, native['exit_code'], independent['exit_code'], rec.get('comparison', {}).get('existing_field_gates_passed'), flush=True)
    after = {str(p.relative_to(ROOT)): sha(p) for p in sources}
    write(out/'source-binding-after.json', after)
    assert before == after


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('output', type=Path)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--attached-levels', nargs='*', type=int, default=[2, 4, 8])
    p.add_argument('--roof-levels', nargs='*', type=int, default=[2, 4])
    p.add_argument('--steps', type=int, default=4)
    a = p.parse_args()
    run(a.output.resolve(), a.binary.resolve(), a.attached_levels, a.roof_levels, a.steps)
