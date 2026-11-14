"""Explicit bounded diagnostic orchestration; never called by production."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
from a04cause_binding import ROOT, OLD, bind
HERE = Path(__file__).resolve().parent
BASE = '1c921b1df0a3da6784476a81c68e6353e39b3678'


def read(p):
    return json.loads(Path(p).read_text(encoding='utf-8-sig'))


def write(p, value):
    with Path(p).open('x', encoding='utf-8', newline='\n') as f:
        json.dump(value, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write('\n')


def sha(p):
    with Path(p).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def reverse_keys(value):
    if isinstance(value, dict):
        return {k: reverse_keys(value[k]) for k in reversed(value)}
    if isinstance(value, list):
        return [reverse_keys(v) for v in value]
    return value


def exact(value):
    if isinstance(value, dict):
        return {k: exact(v) for k, v in value.items() if k != 'elapsed_seconds'}
    if isinstance(value, list):
        return [exact(v) for v in value]
    return value


def invoke(command, log):
    start = time.monotonic()
    command = list(map(str, command))
    with log.open('x', encoding='utf-8') as f:
        p = subprocess.Popen(command, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT)
        timed = False
        try:
            code = p.wait(timeout=int(os.environ.get('A04CAUSE_PROCESS_BUDGET_SECONDS', '1800')))
        except subprocess.TimeoutExpired:
            timed = True
            p.kill()
            code = p.wait()
    result = dict(command=command, pid=p.pid, exit_code=code, timed_out=timed, elapsed_seconds=time.monotonic()-start)
    write(log.with_suffix('.process.json'), result)
    return result


def independent(folder, order):
    legacy = bind()
    mesh, ground, case = [read(folder/(k+'.json')) for k in ('mesh', 'ground', 'case')]
    if order is not None:
        mesh['quadrature_order'] = order
    model = legacy['soil_reference'].Model(mesh, ground, case)
    stream, sink = legacy['study'].sink_file(folder/'independent-fields.jsonl')
    try:
        result = model.solve(sink)
    finally:
        stream.close()
    result['diagnostic_hydraulic_order'] = order or mesh.get('quadrature_order', 3)
    result['input_only_no_native_solution_read'] = True
    write(folder/'independent-result.json', result)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--independent', type=Path)
    p.add_argument('--hydraulic-order', type=int, choices=(3, 5, 7))
    p.add_argument('--steps', nargs=3, type=int, default=(4, 8, 16))
    a = p.parse_args()
    if a.independent:
        independent(a.independent, a.hydraulic_order)
        return
    if a.steps[1] != 2*a.steps[0] or a.steps[2] != 2*a.steps[1] or any(n < 4 or n > 128 or n & (n-1) for n in a.steps):
        raise ValueError('three bounded dyadic levels, 4..128 required')
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    sources = sorted(set(subprocess.check_output(['git', 'ls-files'], cwd=ROOT, text=True).splitlines()) |
                     {str(f.relative_to(ROOT)) for f in HERE.glob('*') if f.suffix == '.py' or f.name == 'SPEC.ja.md'})
    sources = [s for s in sources if Path(s).suffix in ('.py', '.rs', '.toml', '.ps1') or s.endswith('SPEC.ja.md')]
    binding = {s: sha(ROOT/s) for s in sources}
    write(out/'source-binding.json', dict(base=BASE, sources=binding))
    snapshot = out/'source-snapshot'
    for s in binding:
        target = snapshot/s
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT/s, target)
    records = []
    def run(command, name, expected=0):
        result = invoke(command, out/(name+'.log'))
        records.append(dict(name=name, **result))
        write(out/(name+'-checked.json'), dict(expected_exit_code=expected, actual=result))
        if result['exit_code'] != expected or result['timed_out']:
            raise RuntimeError((name, result))
        return result
    try:
        run([sys.executable, OLD/'study.py', '--output', out/'soil', '--steps', *a.steps], 'soil')
        run([sys.executable, OLD/'audit.py', '--run', out/'soil'], 'soil-audit')
        run([sys.executable, OLD/'convergence.py', '--run', out/'soil', '--steps', *a.steps], 'soil-convergence')
        if a.steps == [4, 8, 16] or a.steps == (4, 8, 16):
            for kind in ('dry', 'hydrostatic'):
                run([sys.executable, OLD/'joint_controls.py', '--kind', kind, '--steps', '4', '--output', out/kind], kind)
            suffix = '.exe' if os.name == 'nt' else ''
            solver = out/'soil/bin'/('solver'+suffix)
            run([sys.executable, OLD/'refusals.py', '--solver', solver, '--output', out/'refusals'], 'refusals')
            for label in ('VG', 'saturated'):
                base = out/'soil'/(label+'-n4')
                for variant in ('repeat', 'reverse-keys'):
                    folder = out/(label+'-'+variant)
                    folder.mkdir()
                    for k in ('mesh', 'ground', 'case'):
                        value = read(base/(k+'.json'))
                        write(folder/(k+'.json'), reverse_keys(value) if variant == 'reverse-keys' else value)
                    run([solver, 'solve-one-way-construction', folder/'mesh.json', folder/'ground.json', folder/'case.json', folder/'native.json'], folder.name)
                    if variant == 'reverse-keys':
                        run([sys.executable, Path(__file__), '--output', out, '--independent', folder], folder.name+'-independent')
                for order in (5, 7):
                    folder = out/(label+'-quadrature'+str(order))
                    folder.mkdir()
                    for k in ('mesh', 'ground', 'case'):
                        shutil.copy2(out/'soil'/(label+'-n8')/(k+'.json'), folder/(k+'.json'))
                    run([sys.executable, Path(__file__), '--output', out, '--independent', folder, '--hydraulic-order', str(order)], folder.name)
            folder = out/'boundary-gap'
            folder.mkdir()
            for k in ('mesh', 'ground'):
                shutil.copy2(out/'soil/saturated-n4'/(k+'.json'), folder/(k+'.json'))
            value = read(out/'soil/saturated-n4/case.json')
            value['one_way']['flow']['analysis']['periods'][2]['boundaries'][1]['condition']['value']['head_m'] = 1.8
            value['one_way']['exterior_periods'][3]['reservoirs'][0]['level_m'] = 2.2
            write(folder/'case.json', value)
            run([solver, 'solve-one-way-construction', folder/'mesh.json', folder/'ground.json', folder/'case.json', folder/'native.json'], 'boundary-gap-native')
            run([sys.executable, Path(__file__), '--output', out, '--independent', folder], 'boundary-gap-independent', expected=1)
        changed = [s for s in binding if sha(ROOT/s) != binding[s]]
        if changed:
            raise RuntimeError('source changed during diagnosis: '+repr(changed))
        write(out/'driver-result.json', dict(base=BASE, diagnostic_only=True, passed=False, a04_complete=False,
            execution_complete=True, sources_unchanged=True, steps=a.steps, processes=records,
            contact_mapping=None, joint_physical_clock=None, joint_mass_balance=None,
            excluded_boundary_gap_not_a_refusal_success=(out/'boundary-gap/native.json').exists()))
    except BaseException as e:
        write(out/'failure.json', dict(execution_complete=False, error=repr(e), processes=records))
        raise


if __name__ == '__main__':
    main()
