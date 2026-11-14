"""Small fresh C-04 solves, independent observations, and real PID restarts."""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
import copy
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module
mechanics = load_module('c04_mechanics_reference', HERE/'reference.py')
hydraulics = load_module('c04_hydraulic_reference', HERE/'hydraulic_reference.py')
import fem_reference as fem

BASE = '7ae0a2a669e0e1df93f3dd9791bbeda056334c3e'
LEVELS = (1.25, .25, .9, .4, 1.1, .75)
WATER_NAMES = ('rise_after_shear', 'fall_while_open', 'rise_after_recontact',
               'fall_after_reverse', 'rise_after_reload', 'return_level')
HYDRAULIC_LEVELS = (.3, .5, .1, .4, .2, .45, .3)


def read(p): return json.loads(Path(p).read_text(encoding='utf-8'))
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def convert(v):
    if isinstance(v, np.ndarray): return np.asarray(v, dtype=float).tolist() if v.dtype.kind == 'f' else v.tolist()
    if isinstance(v, np.floating): return float(v)
    if isinstance(v, np.generic): return v.item()
    raise TypeError(type(v).__name__)
def write(p, value):
    with Path(p).open('x', encoding='utf-8') as f:
        json.dump(value, f, default=convert, ensure_ascii=False, allow_nan=False, indent=2)
        f.write('\n')
def strip_clock(v):
    if isinstance(v, dict): return {k: strip_clock(x) for k, x in v.items() if k != 'elapsed_seconds'}
    if isinstance(v, list): return [strip_clock(x) for x in v]
    return v
def invoke(cmd, log, timeout=900):
    cmd = list(map(str, cmd)); started = time.monotonic()
    with Path(log).open('x', encoding='utf-8') as f:
        f.write(repr(cmd)+'\n'); f.flush()
        child = subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT)
        try: code = child.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            child.kill(); child.wait()
            write(Path(log).with_suffix('.failure.json'), dict(complete=False, passed=False, reason='timeout', pid=child.pid, command=cmd))
            raise
    return dict(command=cmd, pid=child.pid, exit_code=code, elapsed_seconds=time.monotonic()-started)


def mechanical_case(mechanical_steps, water_steps):
    case = fem.normalize_stage_case(read(ROOT/'examples/stages-joint-water-cycle.json'))
    case['settings']['initial_increment'] = 1/mechanical_steps
    case['settings']['linear_solver'] = dict(method='dense_lu')
    case['settings']['maximum_free_dofs'] = 1000
    case['verification_reject_once_stages'] = []
    phases = []
    for s, name, level in zip(case['stages'], WATER_NAMES, LEVELS):
        s['initial_increment'] = 1/mechanical_steps; phases.append(s)
        phases.append(dict(kind='hydrostatic_water', name=name, target_level_m=level, initial_increment=1/water_steps))
    case['stages'] = phases
    return case


def hydraulic_case(steps):
    case = copy.deepcopy(read(ROOT/'examples/one-way-mc.json')['flow'])
    for mat in case['materials']:
        mat['saturated_conductivity_m_s'] = (np.eye(3)*.02).tolist()
        mat['retention']['saturated_specific_storage_per_m'] = .02
    case['initial_head'] = dict(head_m=HYDRAULIC_LEVELS[0])
    periods = []
    for phase, level in enumerate(HYDRAULIC_LEVELS):
        for i in range(1, steps+1):
            periods.append(dict(name=f'p{phase}-s{i}', end_time_s=phase+i/steps,
                                boundaries=[dict(face_set=face, condition=dict(kind='head', value=dict(head_m=level))) for face in ('bottom', 'top')]))
    case['analysis'] = dict(kind='transient', periods=periods)
    case['settings'].update(initial_time_step_s=1/steps, maximum_time_step_s=1/steps,
                            minimum_time_step_s=1/(1024*steps), maximum_steps=100000,
                            store_integration_points=True, maximum_free_faces=1000, linear_solver=dict(method='dense_lu'))
    return case


def schedule(case):
    origin = case['settings']['initial_increment']
    return [dict(name='geostatic', kind='initial', n=round(1/origin))]+[
        dict(name=s['name'], kind=s['kind'], n=round(1/s['initial_increment'])) for s in case['stages']]


def initialize(out, probe, cli, element, mechanical_steps):
    out.mkdir(parents=True, exist_ok=False)
    mesh = read(ROOT/('examples/joint-water-column-mesh.json' if element == 'tet4' else 'examples/joint-water-column-tet10-mesh.json'))
    ground = fem.normalize_ground(read(ROOT/'examples/ground-joint-water-mc.json'))
    ground['initialization']['geostatic']['hydrostatic_face_sets'] = ['exterior']
    write(out/'mesh.json', mesh); write(out/'ground.json', ground)
    write(out/'hydraulic-mesh.json', read(ROOT/'examples/seepage-column-tet4-mesh.json'))
    for label, path in (('probe', probe), ('solver', cli)):
        shutil.copy2(path, out/(label+('.exe' if path.suffix == '.exe' else '')))
    write(out/'identity.json', dict(starting_commit=BASE, element=element, mechanical_steps=mechanical_steps,
          binaries={p.name: sha(p) for p in out.iterdir() if p.name in ('probe', 'probe.exe', 'solver', 'solver.exe')},
          inputs={p.name: sha(p) for p in out.glob('*.json')}, production_changed=False,
          hydrostatic_clock='Each reporting stage has duration 1 s by case definition. Production is quasistatic and stores progress/level, not physical time or water storage.',
          hydraulic_clock='Richards physical seconds; seven fixed boundary periods of 1 s each. Subperiods only expose every accepted step.',
          hydraulic_joint_coupling=False))


def binary(out, name):
    p = out/(name+('.exe' if (out/(name+'.exe')).exists() else ''))
    assert sha(p) == read(out/'identity.json')['binaries'][p.name]
    return p


def native(out, folder, mechanical_steps, water_steps, timeout, whole=False):
    folder.mkdir(parents=True, exist_ok=False)
    case = mechanical_case(mechanical_steps, water_steps); write(folder/'case.json', case)
    proc = invoke([binary(out, 'probe'), out/'mesh.json', out/'ground.json', folder/'case.json', folder/'raw']+(['--whole'] if whole else []), folder/'solve.log', timeout)
    write(folder/'process.json', proc)
    if proc['exit_code'] != 0: raise RuntimeError(('native incomplete', folder, proc['exit_code']))
    result = read(folder/'raw/run.json'); assert result['complete']
    phases = schedule(case)
    groups = [result['geostatic']['attempts']]+[[a for a in result['attempts'] if a['stage'] == s['name']] for s in case['stages']]
    audit = []
    for phase, attempts in zip(phases, groups):
        assert all(a['accepted'] for a in attempts), ('cutback', phase['name'])
        assert [a['progress_to'] for a in attempts] == [i/phase['n'] for i in range(1, phase['n']+1)], ('changed schedule', phase['name'])
        force = max(a['relative_force_residual'] for a in attempts)
        linear = max((s['true_relative_residual'] for a in attempts for s in a['linear_solves']), default=0.)
        assert force <= case['settings']['relative_force_tolerance'] and linear <= case['settings']['relative_linear_tolerance']
        audit.append(dict(**phase, maximum_force=force, maximum_true_linear=linear))
    write(folder/'schedule.json', audit)
    if whole: return
    model = mechanics.Model(read(out/'mesh.json'), read(out/'ground.json'), case)
    coords = [dict(phase=phase, index=i, stage=p['name'], progress=i/p['n']) for phase, p in enumerate(phases) for i in range(1, p['n']+1)]
    with (folder/'raw/accepted.jsonl').open() as sf, (folder/'raw/constitutive.jsonl').open() as df, (folder/'fields.jsonl').open('x') as ff:
        for count, (coord, sline, dline) in enumerate(itertools.zip_longest(coords, sf, df), 1):
            assert coord is not None and sline is not None and dline is not None, 'missing accepted state or diagnostic'
            state = json.loads(sline); diagnostic = json.loads(dline)
            assert state.get('pore_pressure_kpa') is None, 'hydrostatic input must not be relabelled transient GP state'
            model.set_level(state['hydrostatic_water_level_m'])
            j = state['joints']['points']; v = lambda key: np.array([p[key] for p in j])
            traction = v('traction_kpa'); load = np.asarray(state['applied_nodal_load_kn'])
            f = model.fields(np.asarray(state['displacement_m']).ravel(), state['material_states'], v('jump_m'), traction, v('plastic_slip_m'),
                             v('accumulated_plastic_slip_m'), v('plastic_dissipation_kpa_m'), v('numerical_dissipation_kpa_m'), v('endpoint_work_kpa_m'), load,
                             [mechanics.old.CONTACT.index(p['status']) for p in j])
            internal = np.asarray(diagnostic['internal'])
            actual_bulk = internal-(model.internal(state['material_states'], traction)-model.pore_force)
            f['bulk_water_force'] = actual_bulk.reshape(-1, 3)
            f['reaction'] = (internal-load).reshape(-1, 3)
            f['physical_pressure'] = np.asarray(state['joint_water']['physical_pressure_kpa'])
            f['water_load'] = np.asarray(state['joint_water']['nodal_load_kn']).reshape(-1, 3)
            holding = (model.seed_force-(model.dry_base+model.water.hydrostatic(model.active, .75)['load']).ravel()-model.exterior_water.ravel())*(1-coord['progress']) if coord['phase'] == 0 else np.zeros(model.nd)
            f['exterior_water_load'] = (load-holding-state['boundary_factor']*model.boundary).reshape(-1, 3)-model.dry_base-f['water_load']
            f['applied_load'] = load.reshape(-1, 3)
            f['total_normal_pressure'] = f['effective_contact_pressure']+f['physical_pressure']
            row = dict(stage=coord['stage'], progress=coord['progress'], factor=state['boundary_factor'], level_m=model.level, time_s=coord['phase']+coord['progress'], relative_force_residual=diagnostic['relative_force_residual'])
            ff.write(json.dumps(dict(row=row, fields=f), default=convert, allow_nan=False)+'\n')
    assert count == sum(p['n'] for p in phases)
    print('native', folder.name, count, flush=True)


def independent(out, folder, mechanical_steps, water_steps):
    folder.mkdir(parents=True, exist_ok=False)
    case = mechanical_case(mechanical_steps, water_steps); write(folder/'case.json', case)
    model = mechanics.Model(read(out/'mesh.json'), read(out/'ground.json'), case)
    start = time.monotonic()
    with (folder/'fields.jsonl').open('x') as f:
        def sink(row, fields):
            f.write(json.dumps(dict(row=row, fields=fields), default=convert, allow_nan=False)+'\n'); f.flush()
        result = mechanics.run(model, mechanical_steps, case['settings']['material_maximum_substep'], case['joint_system']['maximum_increment_m'], sink)
    write(folder/'result.json', result)
    write(folder/'process.json', dict(pid=os.getpid(), elapsed_seconds=time.monotonic()-start, reference_uses_native_results=False))
    print('independent', folder.name, len(result['rows']), flush=True)


def hydraulic(out, folder, steps, timeout):
    folder.mkdir(parents=True, exist_ok=False)
    case = hydraulic_case(steps); write(folder/'case.json', case)
    mesh_path = out/'hydraulic-mesh.json'
    proc = invoke([binary(out, 'solver'), 'solve-richards', mesh_path, folder/'case.json', folder/'native.json'], folder/'solve.log', timeout)
    write(folder/'process.json', proc)
    if proc['exit_code'] != 0: raise RuntimeError(('Richards incomplete', proc))
    native_run = read(folder/'native.json'); assert native_run['complete']
    assert all(a['accepted'] for a in native_run['attempts'])
    assert len(native_run['attempts']) == 7*steps
    assert len(native_run['frames']) == 7*steps+1
    start = time.monotonic()
    independent_run = hydraulics.Model(read(mesh_path), case).solve()
    write(folder/'independent.json', independent_run)
    write(folder/'independent-process.json', dict(pid=os.getpid(), elapsed_seconds=time.monotonic()-start, reference_uses_native_results=False))
    print('hydraulic', steps, len(independent_run['rows']), flush=True)


def restart(out, source, target, timeout):
    target.mkdir(parents=True, exist_ok=False)
    direct = strip_clock(read(source/'raw/run.json'))
    base_pid = read(source/'process.json')['pid']; cps = sorted((source/'raw').glob('checkpoint-*.json'))
    assert len(cps) == len(schedule(read(source/'case.json')))
    rows = []
    for cp in cps:
        saved = read(cp); state = saved['run']['last_committed']
        assert state['stage_progress'] == .5
        suffix = target/cp.stem
        proc = invoke([binary(out, 'probe'), out/'mesh.json', out/'ground.json', source/'case.json', suffix, '--resume', cp], target/(cp.stem+'.log'), timeout)
        assert proc['exit_code'] == 0 and proc['pid'] != base_pid
        assert strip_clock(read(suffix/'run.json')) == direct
        stop = int(cp.stem.split('-')[-1]); count = 0
        with (source/'raw/accepted.jsonl').open() as a, (suffix/'accepted.jsonl').open() as b:
            for _ in range(stop): next(a)
            for x, y in itertools.zip_longest(a, b):
                assert x is not None and y is not None and json.loads(x) == json.loads(y)
                count += 1
        rows.append(dict(checkpoint=cp.name, checkpoint_sha256=sha(cp), saved_stage_index=state['stage_index'], saved_progress=.5,
                         baseline_pid=base_pid, resume_pid=proc['pid'], suffix_states_exact=count, entire_run_exact_except_elapsed=True, process=proc))
    write(target/'verification.json', dict(complete=True, real_file_separate_pid=True, comparisons=rows))


def hydraulic_restart(out, source, target, timeout):
    target.mkdir(parents=True, exist_ok=False)
    direct = read(source/'native.json'); n = len(direct['attempts']); rows = []
    steps = round(1/direct['case']['settings']['maximum_time_step_s'])
    for count in (steps//2, 2*steps+steps//2, 5*steps+steps//2):
        cp = target/f'checkpoint-{count:04}.json'; resumed = target/f'resumed-{count:04}.json'
        save = invoke([binary(out, 'solver'), 'checkpoint-richards', out/'hydraulic-mesh.json', source/'case.json', count, cp], target/f'save-{count}.log', timeout)
        load = invoke([binary(out, 'solver'), 'resume-richards', out/'hydraulic-mesh.json', source/'case.json', cp, resumed], target/f'resume-{count}.log', timeout)
        assert save['exit_code'] == load['exit_code'] == 0 and save['pid'] != load['pid']
        assert strip_clock(read(resumed)) == strip_clock(direct)
        rows.append(dict(accepted_count=count, saved_time_s=read(cp)['run']['last_committed']['time_s'], checkpoint_sha256=sha(cp),
                         checkpoint_pid=save['pid'], resume_pid=load['pid'], baseline_pid=read(source/'process.json')['pid'],
                         entire_run_exact_except_elapsed=True, suffix_frames_exact=n-count, save_process=save, resume_process=load))
    write(target/'verification.json', dict(complete=True, real_file_separate_pid=True, comparisons=rows))


def refusals(out, folder, timeout):
    folder.mkdir(parents=True, exist_ok=False)
    case = read(ROOT/'examples/one-way-mc.json')
    case['flow'] = hydraulic_case(4)
    case['exterior_periods'] = [dict(name=p['name'], reservoirs=[]) for p in case['flow']['analysis']['periods']]
    wet = mechanical_case(32, 8)
    case['mechanics']['stages'] = []
    case['mechanics']['pressure_on_original_surface_only'] = True
    case['mechanics']['settings']['kinematics'] = 'explicit'
    case['mechanics']['boundary'] = wet['boundary']
    case['mechanics']['joint_system'] = wet['joint_system']
    ground = read(ROOT/'examples/ground-one-way-mc.json')
    write(folder/'ground.json', ground)
    rows = []
    for policy in ('continuous', 'soil', 'prescribed'):
        c = copy.deepcopy(case); config = c['mechanics']['joint_system']
        config['joints'][0]['water'] = dict(kind=policy)
        if policy == 'soil':
            config['structure_bodies'] = ['UpperClay']; config['joints'][0]['water']['soil_side'] = 'a'
        if policy == 'prescribed': config['joints'][0]['water']['level_m'] = .75
        path = folder/(policy+'.json'); write(path, c)
        log = folder/(policy+'.log'); output = folder/(policy+'-result.json')
        proc = invoke([binary(out, 'solver'), 'solve-one-way', out/'mesh.json', folder/'ground.json', path, output], log, timeout)
        text = log.read_text()
        assert proc['exit_code'] != 0 and 'ONE_WAY requires dry-reference weights' in text and not output.exists(), text
        rows.append(dict(policy=policy, supported=False, refusal_verified=True, process=proc, diagnostic=text.splitlines()[-1]))
    write(folder/'verification.json', dict(Richards_joint_ONE_WAY_accepted=False, cases=rows,
          bidirectional=dict(supported=False, executed=False, evidence='No bidirectional controller/schema/CLI in src/main.rs, src/one_way.rs or src/one_way_construction.rs; no new feedback model was added.'),
          birth_death_generations=dict(supported=False, executed=False, reason='ExternalMechanics rejects JOINT before any hydraulic-to-interface mapping or contact generation exists.'),
          pr19_exponential_Richards=dict(supported=False, executed=False, reason='Production retention input is van Genuchten/Mualem, not PR #19 exponential conductivity/storage; no constitutive substitution.')))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mode', choices=['init', 'native', 'independent', 'hydraulic', 'restart', 'hydraulic-restart', 'refusals', 'whole'], required=True)
    p.add_argument('--output', type=Path, required=True); p.add_argument('--folder', default='path')
    p.add_argument('--source', default='native-w32'); p.add_argument('--mechanical-steps', type=int, default=32)
    p.add_argument('--water-steps', type=int, default=8); p.add_argument('--hydraulic-steps', type=int, default=4)
    p.add_argument('--element', choices=['tet4', 'tet10'], default='tet4'); p.add_argument('--timeout', type=int, default=900)
    p.add_argument('--probe', type=Path, default=ROOT/'target/release/examples/wet_joint_transient_probe')
    p.add_argument('--solver', type=Path, default=ROOT/'target/release/soil-fem-3d')
    a = p.parse_args(); out = a.output.resolve(); folder = out/a.folder
    try:
        if a.mode == 'init': initialize(out, a.probe.resolve(), a.solver.resolve(), a.element, a.mechanical_steps)
        elif a.mode == 'independent': independent(out, folder, a.mechanical_steps, a.water_steps)
        elif a.mode == 'native' or a.mode == 'whole': native(out, folder, a.mechanical_steps, a.water_steps, a.timeout, a.mode == 'whole')
        elif a.mode == 'hydraulic': hydraulic(out, folder, a.hydraulic_steps, a.timeout)
        elif a.mode == 'restart': restart(out, out/a.source, folder, a.timeout)
        elif a.mode == 'hydraulic-restart': hydraulic_restart(out, out/a.source, folder, a.timeout)
        elif a.mode == 'refusals': refusals(out, folder, a.timeout)
    except Exception as e:
        if out.exists(): write(out/f'failure-{a.mode}-{time.time_ns()}.json', dict(complete=False, passed=False, error=repr(e)))
        raise


if __name__ == '__main__': main()
