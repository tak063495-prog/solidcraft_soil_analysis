"""Explicit bounded fresh equilibrium matrix, with raw states and process logs."""
import argparse
import copy
import hashlib
import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from a03follow_binding import ROOT, HERE, OLD, study, fields, reference
from a03follow_controls import fixture, controls, Sink, layout, canonical, compare


def write(path, value):
    with Path(path).open('x', encoding='utf-8', newline='\n') as f:
        json.dump(value, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write('\n')


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def command(args, folder, name, minimum=1, reverse=False):
    env = dict(os.environ)
    env.pop('A03FOLLOW_REVERSE_FACE_POINTS', None)
    env['A03FOLLOW_MINIMUM'] = str(minimum)
    if reverse:
        env['A03FOLLOW_REVERSE_FACE_POINTS'] = '1'
    start = time.monotonic()
    with (folder/(name+'.log')).open('x', encoding='utf-8') as log:
        p = subprocess.Popen(list(map(str, args)), cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
        timed = False
        try:
            code = p.wait(timeout=int(os.environ.get('A03FOLLOW_PROCESS_SECONDS', '3600')))
        except subprocess.TimeoutExpired:
            timed = True
            p.kill()
            code = p.wait()
    value = dict(command=list(map(str, args)), pid=p.pid, exit_code=code, timed_out=timed,
                 elapsed_seconds=time.monotonic()-start, minimum_parts=minimum, reverse_face_points=reverse)
    write(folder/(name+'-process.json'), value)
    if code or timed:
        raise RuntimeError(('numerically incomplete; retain all prefix files', value))
    return value


def make_plan():
    specs = {}
    def add(role, n, shape='wide', stress='k05', minimum=1, element=False, face=False, reverse=False, volume=5):
        name = f'{shape}-{stress}-n{n}-m{minimum}-e{int(element)}f{int(face)}q{int(reverse)}-gl{volume}'
        if name not in specs:
            specs[name] = dict(name=name, n=n, shape=shape, stress=stress, minimum=minimum,
                               element=element, face=face, reverse=reverse, volume=volume,
                               engines=['independent'] if volume == 7 else ['native', 'independent'], roles=[])
        specs[name]['roles'].append(role)
    for shape, stress, n in itertools.product(('wide', 'tall'), ('k05', 'k08'), (4, 8, 16)):
        add('geometry', n, shape, stress)
    for shape, stress in itertools.product(('wide', 'tall'), ('k05', 'k08')):
        add('volume', 8, shape, stress, volume=7)
    for n, minimum in itertools.product((4, 8, 16, 32), (1, 2, 4)):
        add('cross', n, minimum=minimum)
    for e, f, q in itertools.product((False, True), repeat=3):
        add('order', 8, element=e, face=f, reverse=q)
    return dict(a03_complete=False, diagnostic_only=True, planned_before_solves=True,
                geometry='54 nodes / 12 Tet10 / 324 material GP / 32 contact GP',
                common_cross_coordinates_per_stage=4, finest_is_exact=False, specs=list(specs.values()))


def snapshot(out):
    paths = subprocess.check_output(['git', 'ls-files', 'src', 'tools', 'examples', 'tests',
                                    'Cargo.toml', 'Cargo.lock', 'rust-toolchain.toml', '.github'], cwd=ROOT, text=True).splitlines()
    paths = [p for p in paths if Path(p).suffix in ('.rs', '.py', '.toml', '.yml') or Path(p).name == 'Cargo.lock']
    paths += [p.relative_to(ROOT).as_posix() for d in (OLD, HERE) for p in d.glob('*.py')]
    paths += [p.relative_to(ROOT).as_posix() for p in HERE.glob('SPEC.ja.md')]
    paths = sorted(set(paths))
    binding = {p: sha(ROOT/p) for p in paths}
    write(out/'source-binding.json', binding)
    for p in paths:
        target = out/'source-snapshot'/p
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT/p, target)
    return binding


def independent(folder, spec):
    with controls(spec['volume'], spec['reverse'], spec['minimum']):
        model = reference.Model(*[study.read(p) for p in study.paths(folder)])
        mapping = layout(model)
        write(folder/'independent-layout.json', mapping)
        with (folder/'independent-fields.jsonl').open('x', encoding='utf-8') as f, (folder/'independent-prepared.jsonl').open('x', encoding='utf-8') as p:
            result = reference.solve(model, spec['n'], Sink(f), Sink(p, work=False))
    result.update(volume_order=spec['volume'], minimum_parts=spec['minimum'], reverse_faces=spec['reverse'])
    write(folder/'independent-result.json', result)


def native(folder, spec, probe, solver):
    n = spec['n']
    proc = command([probe, *study.paths(folder), folder/'native'], folder, 'native', spec['minimum'], spec['reverse'])
    run = study.read(folder/'native/run.json')
    assert run['complete'] and run['provenance']['solver_binary_sha256'] == sha(probe)
    groups = [run['geostatic']['attempts']] + [[a for a in run['attempts'] if a['stage'] == s['name']] for s in run['case']['stages']]
    for attempts in groups:
        assert len(attempts) == n and all(a['accepted'] for a in attempts)
        assert [a['progress_to'] for a in attempts] == [i/n for i in range(1, n+1)]
        assert all(a['relative_force_residual'] <= run['case']['settings']['relative_force_tolerance'] for a in attempts)
        assert all(l['true_relative_residual'] <= 1e-10 for a in attempts for l in a['linear_solves'])
    audits = []
    with controls(5, spec['reverse'], spec['minimum']):
        model = reference.Model(*[study.read(p) for p in study.paths(folder)])
        # Legacy layout audit writes an independent geometry file before any comparison.
        geometry = study.geometry(model, run, folder)
        (folder/'independent-layout.json').rename(folder/'native-public-layout.json')
        mapping = layout(model)
        write(folder/'native-layout.json', mapping)
        prep = run['transition_audits'][0]
        selected = set(prep['selected_elements'])
        fresh = run['case']['stages'][0]['states']['UpperClay']
        assert all(study.exact(s) == study.exact(fresh) for q, s in zip(run['integration_points'], prep['prepared_state']['material_states']) if q['element_index'] in selected)
        names = ['geostatic'] + [s['name'] for s in run['case']['stages']]
        with (folder/'native-fields.jsonl').open('x', encoding='utf-8') as f:
            sink = Sink(f)
            with (folder/'native/accepted.jsonl').open(encoding='utf-8') as a, (folder/'native/constitutive.jsonl').open(encoding='utf-8') as b:
                count = 0
                for state, diagnostic in itertools.zip_longest(a, b):
                    assert state is not None and diagnostic is not None
                    state, diagnostic = json.loads(state), json.loads(diagnostic)
                    if count == n:
                        sink.placement(prep['prepared_state']['displacement_m'])
                    values, audit = fields.native_fields(model, state, diagnostic, run)
                    assert not np.any(values['material_return']) and not np.any(values['plastic_work'])
                    if count >= 2*n:
                        assert np.all(values['contact_status'] == 1) and not np.any(values['plastic_dissipation'])
                    sink(fields.row(names[count//n], (count % n+1)/n), values)
                    audits.append(dict(row=fields.row(names[count//n], (count % n+1)/n), **audit))
                    count += 1
                assert count == 4*n
        with (folder/'native-prepared.jsonl').open('x', encoding='utf-8') as f:
            sink = Sink(f, work=False)
            for entry in study.read(folder/'native/preparations.json'):
                audit = next(a for a in run['transition_audits'] if a['stage'] == entry['stage'])
                values, check = fields.native_fields(model, audit['prepared_state'], entry['diagnostic'], run, True)
                sink(fields.row(entry['stage'], 0.), values)
    write(folder/'native-audits.json', dict(layout=geometry, rows=audits, fresh_state_exact=True))
    whole = command([probe, *study.paths(folder), folder/'whole', '--whole'], folder, 'whole', spec['minimum'], spec['reverse'])
    assert study.exact(run) == study.exact(study.read(folder/'whole/run.json'))
    standard = None
    if spec['n'] == 8 and spec['minimum'] == 1 and not any(spec[k] for k in ('element', 'face', 'reverse')):
        standard = command([solver, 'solve-stages', *study.paths(folder), folder/'production-run.json'], folder, 'production')
        left, right = copy.deepcopy(run), study.read(folder/'production-run.json')
        for v in (left, right):
            v['provenance'].pop('solver_binary_sha256')
        assert study.exact(left) == study.exact(right)
    return dict(process=proc, whole_process=whole, production_process=standard, states=4*n,
                observed_unobserved_exact_except_elapsed=True, production_default_exact=standard is not None)


def run(out, build):
    out = out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    plan = make_plan()
    write(out/'plan.json', plan)
    binding = snapshot(out)
    build = build.resolve()
    suffix = '.exe' if os.name == 'nt' else ''
    probe = build/('target/release/a03follow_probe'+suffix)
    solver = ROOT/('target/release/soil-fem-3d'+suffix)
    write(out/'binaries.json', {name: dict(path=str(path), sha256=sha(path)) for name, path in [('probe', probe), ('production', solver)]})
    shutil.copy2(build/'binding.json', out/'isolated-build-binding.json')
    shutil.copy2(build/'stages.patch', out/'stages.patch')
    shutil.copy2(build/'elements.patch', out/'elements.patch')
    paths = []
    for spec in plan['specs']:
        folder = out/spec['name']
        study.inputs(folder, fixture(spec['n'], spec['shape'], spec['stress'], spec['element'], spec['face']))
        write(folder/'spec.json', spec)
        if 'native' in spec['engines']:
            value = native(folder, spec, probe, solver)
            paths.append(dict(path=spec['name'], engine='native', **spec, **value))
        proc = command([sys.executable, HERE/'a03follow_run.py', '--independent', folder], folder, 'independent', spec['minimum'], spec['reverse'])
        paths.append(dict(path=spec['name'], engine='independent', **spec, process=proc, states=4*spec['n']))
        if 'native' in spec['engines']:
            a, b = (list(fields.load_fields(folder/(engine+'-fields.jsonl'))) for engine in ('native', 'independent'))
            # Same-layout checks retain the original acceptance gates; no new all-field gate.
            compared = compare(a, b)
            for key in ('u', 'stress', 'traction'):
                assert compared['fields'][key]['relative_l2'] < 1e-7
            for key in ('internal', 'applied_load', 'effective_internal', 'pore_internal', 'contact_internal', 'gravity_load', 'surcharge_load', 'exterior_water_load', 'joint_water_load', 'holding_load'):
                assert compared['fields'][key]['maximum_absolute'] < 1e-10
            assert all(v['mismatches'] == 0 for v in compared['fields'].values() if 'mismatches' in v)
            write(folder/'implementation-comparison.json', compared)
        print(spec['name'], 'complete', flush=True)
    assert binding == {p: sha(ROOT/p) for p in binding}
    result = dict(passed=False, a03_complete=False, diagnostic_only=True, paths=paths,
                  formal_paths=len(paths), formal_states=sum(p['states'] for p in paths),
                  compatibility_paths=sum(p['engine'] == 'native' for p in paths),
                  compatibility_states=sum(p['states'] for p in paths if p['engine'] == 'native'),
                  standard_cli_paths=sum(p['engine'] == 'native' and p['production_process'] is not None for p in paths),
                  standard_cli_states=sum(p['states'] for p in paths if p['engine'] == 'native' and p['production_process'] is not None),
                  sources_unchanged_after_solves=True, refusals_unchanged_not_resolved=True)
    write(out/'driver-result.json', result)
    print(json.dumps({k: v for k, v in result.items() if k != 'paths'}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path)
    p.add_argument('--build', type=Path)
    p.add_argument('--independent', type=Path)
    a = p.parse_args()
    if a.independent:
        independent(a.independent, study.read(a.independent/'spec.json'))
    else:
        if a.output is None or a.build is None:
            p.error('--output and --build required')
        run(a.output, a.build)
