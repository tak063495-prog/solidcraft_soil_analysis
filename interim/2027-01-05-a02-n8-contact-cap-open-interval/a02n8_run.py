"""Explicit new equilibrium solves, bounded to the two requested A-02 probes."""
import argparse
import itertools
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_key] = '1'

from a02n8_binding import ROOT, HERE, FOLLOW, study, fields, Predictor, serial
from a02n8_common import write, sha
import a02n8_reference as controlled

BASE = '612f0f64c704f0ef0dee950a81a10f7e312a745c'


def schedule(n, divisions):
    """Exact binary endpoints; later open intervals and every other stage stay fixed."""
    if n not in (4, 8) or divisions not in (1, 2, 4, 8):
        raise ValueError('bounded N/divisions only')
    return {'open': [i/(n*divisions) for i in range(1, divisions+1)]
                    + [i/n for i in range(2, n+1)]}


def native(folder, binary, n, grid):
    values = [study.read(p) for p in study.input_paths(folder)]
    process = study.invoke([binary, *study.input_paths(folder), folder/'native'], folder/'native.log')
    write(folder/'native-process.json', process)
    if process['exit_code'] or process['timed_out']:
        raise RuntimeError('native numerically incomplete; retain accepted prefix')
    result = study.read(folder/'native/run.json')
    assert result['complete']
    names = ['geostatic'] + [s['name'] for s in values[2]['stages']]
    expected = [(name, p) for name in names for p in grid.get(name, [i/n for i in range(1, n+1)])]
    groups = [result['geostatic']['attempts']] + [
        [a for a in result['attempts'] if a['stage'] == s] for s in names[1:]]
    for name, attempts in zip(names, groups):
        assert all(a['accepted'] for a in attempts), 'unexpected cutback retained as incomplete'
        assert [a['progress_to'] for a in attempts] == [p for s, p in expected if s == name]
        assert all(a['relative_force_residual'] <= values[2]['settings']['relative_force_tolerance'] for a in attempts)
        assert all(x['true_relative_residual'] <= 1e-10 for a in attempts for x in a['linear_solves'])
    model = controlled.Model(*values)
    with (folder/'native-fields.jsonl').open('x', encoding='utf-8', newline='\n') as stream:
        sink = fields.Sink(stream)
        with (folder/'native/accepted.jsonl').open() as states, (folder/'native/constitutive.jsonl').open() as diagnostics, (folder/'native/structure-reports.jsonl').open() as reports:
            actual = []
            for state, diagnostic, report in itertools.zip_longest(states, diagnostics, reports):
                assert all(x is not None for x in (state, diagnostic, report))
                state, diagnostic, report = map(json.loads, (state, diagnostic, report))
                phase = state['stage_index']
                coord = (names[phase], diagnostic['progress'])
                assert diagnostic['stage'] == phase and all(state['active_elements'])
                actual.append(coord)
                sink(dict(stage=coord[0], progress=coord[1]), study.native_fields(model, state, diagnostic, report))
            assert actual == expected
    direct = study.invoke([binary, *study.input_paths(folder), folder/'whole', '--whole'], folder/'whole.log')
    write(folder/'whole-process.json', direct)
    assert direct['exit_code'] == 0
    assert study.exact(study.read(folder/'whole/run.json')) == study.exact(result)
    for audit in result['transition_audits']:
        if audit['stage'] in ('birth', 'rebirth'):
            for point in audit['prepared_state']['joints']['points']:
                for key in ('plastic_slip_m', 'accumulated_plastic_slip_m', 'plastic_dissipation_kpa_m', 'numerical_dissipation_kpa_m', 'endpoint_work_kpa_m'):
                    assert not any(fields.np.asarray(point[key]).ravel())
    death = next(i for i, s in enumerate(values[2]['stages']) if s['kind'] == 'joint_death')
    archive = result['last_committed']['joints']['lifecycle'][0]['archived'][0]['points']
    assert study.exact(archive) == study.exact(list(enumerate(result['frames'][death]['state']['joints']['points'])))
    return len(actual)


def worker(out, build, n, cap, divisions, engine, predictor, repartition, parts):
    values = study.fixture('stress_free', n, True)
    assert len(values[0]['tetrahedra']) == 24 and len(values[0]['nodes']) == 81
    assert values[2]['joint_system'].get('global_step_control') is None
    values[2]['joint_system']['maximum_increment_m'] = cap
    grid = schedule(n, divisions) if divisions > 1 else {}
    if divisions > 1:
        values[2]['settings']['minimum_increment'] = 1/(n*divisions)
    study.make_inputs(out, values)
    write(out/'global-input-endpoints.json', grid)
    paths = list(HERE.glob('*.py')) + list(FOLLOW.glob('*.py'))
    paths += list((ROOT/'interim/2026-12-21-a02-joint-birth-cause-diagnosis').glob('*.py'))
    paths += list((ROOT/'interim/2026-11-25-joint-birth-projection-a02').glob('*.py'))
    binding = dict(base=BASE, n=n, maximum_contact_increment_m=cap, divisions=divisions,
                   engine=engine, predictor=predictor, repartition=repartition, minimum_parts=parts,
                   source_sha256={str(p.relative_to(ROOT)): sha(p) for p in sorted(paths)},
                   build_binding=study.read(build/'binding.json'),
                   input_sha256={p.name: sha(p) for p in study.input_paths(out) + [out/'global-input-endpoints.json']})
    write(out/'binding.json', binding)
    exe = '.exe' if os.name == 'nt' else ''
    if engine == 'native':
        os.environ.update(A02N8_P=predictor, A02N8_R=repartition, A02N8_MINIMUM=str(parts), A02N8_TRACE='1')
        os.environ.pop('A02N8_SCHEDULE', None)
        if grid:
            os.environ['A02N8_SCHEDULE'] = str(out/'global-input-endpoints.json')
        count = native(out, build/'target/release'/('a02n8_probe'+exe), n, grid)
    else:
        service = Predictor(build/'target/release'/('a02n8_predictor'+exe)) if predictor == 'native' else None
        try:
            with (out/'independent-fields.jsonl').open('x', encoding='utf-8', newline='\n') as stream, (out/'events.jsonl').open('x', encoding='utf-8', newline='\n') as events:
                def observe(v):
                    keep = ('kind', 'stage', 'progress', 'iteration', 'relative', 'counts', 'target', 'traction', 'gap', 'reset', 'cut', 'used', 'trust', 'accepted', 'projected_load', 'original_linear_relative')
                    events.write(json.dumps({k: x for k, x in v.items() if k in keep}, default=serial, allow_nan=False)+'\n')
                    events.flush()
                result = controlled.run(controlled.Model(*values), n, fields.Sink(stream), minimum_parts=parts,
                                        predictor=predictor, repartition=repartition, service=service,
                                        observer=observe, stage_progress=grid)
            write(out/'independent-result.json', result)
            count = result['accepted_states']
        finally:
            if service:
                service.close()
    assert count == 7*n + divisions-1
    assert all(sha(ROOT/p) == digest for p, digest in binding['source_sha256'].items())
    write(out/'completion.json', dict(complete=True, accepted_states=count,
                                     observer_free_compatibility_states=count if engine == 'native' else 0,
                                     source_binding_unchanged=True, a02_passed=False))


def jobs():
    result = []
    for cap, label in ((.0005, 'cap'), (.05, 'normal')):
        for e, p, r in itertools.product(('native', 'independent'), repeat=3):
            result.append(dict(path=f'exchange/{label}/{e}-p{p}-r{r}', n=8, cap=cap, divisions=1,
                               engine=e, predictor=p, repartition=r, parts=1))
    for n, b, e in itertools.product((4, 8), (2, 4, 8), ('native', 'independent')):
        result.append(dict(path=f'open/n{n}-b{b}/{e}', n=n, cap=.05, divisions=b,
                           engine=e, predictor=e, repartition=e, parts=4))
    return result


def batch(out, build):
    out.mkdir(parents=True, exist_ok=False)
    rows = []
    for job in jobs():
        target = out/job['path']
        target.parent.mkdir(parents=True, exist_ok=True)
        command = [sys.executable, str(Path(__file__).resolve()), 'worker', '--output', str(target), '--build', str(build)]
        for key in ('n', 'cap', 'divisions', 'engine', 'predictor', 'repartition', 'parts'):
            command += ['--'+key, str(job[key])]
        start = time.monotonic()
        with target.with_suffix('.log').open('x', encoding='utf-8') as log:
            process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
            timed_out = False
            try:
                code = process.wait(timeout=960)
            except subprocess.TimeoutExpired:
                timed_out = True
                process.kill()
                code = process.wait()
        row = dict(**job, pid=process.pid, exit_code=code, timed_out=timed_out, elapsed_seconds=time.monotonic()-start)
        if (target/'completion.json').exists():
            row.update(study.read(target/'completion.json'))
        else:
            row.update(complete=False, classification='numerically_incomplete')
            file = target/('native/accepted.jsonl' if job['engine'] == 'native' else 'independent-fields.jsonl')
            row['accepted_prefix_lines'] = sum(1 for _ in file.open()) if file.exists() else 0
        write(target.with_suffix('.process.json'), row)
        rows.append(row)
        print(job['path'], code, row['complete'], round(row['elapsed_seconds'], 3), flush=True)
    write(out/'paths.json', rows)
    if not all(r['complete'] and r['exit_code'] == 0 for r in rows):
        raise RuntimeError('incomplete numerical paths retained; not refusal, failure of physics, or acceptance')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=('batch', 'worker'))
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--build', type=Path, required=True)
    p.add_argument('--n', type=int, choices=(4, 8))
    p.add_argument('--cap', type=float, choices=(.0005, .05))
    p.add_argument('--divisions', type=int, choices=(1, 2, 4, 8))
    p.add_argument('--parts', type=int, choices=(1, 4))
    for key in ('engine', 'predictor', 'repartition'):
        p.add_argument('--'+key, choices=('native', 'independent'))
    a = p.parse_args()
    if a.mode == 'batch':
        batch(a.output.resolve(), a.build.resolve())
    else:
        try:
            worker(a.output.resolve(), a.build.resolve(), a.n, a.cap, a.divisions, a.engine, a.predictor, a.repartition, a.parts)
        except BaseException as error:
            if a.output.exists():
                write(a.output/'failure.json', dict(error=repr(error), traceback=traceback.format_exc(), classification='numerically_incomplete', a02_passed=False))
            raise
