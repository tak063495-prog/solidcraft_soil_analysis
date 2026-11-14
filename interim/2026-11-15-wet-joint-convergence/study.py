"""Separate space, global steps, contact partitions and volume rules; never certify P6.

All solves start from inputs. Historical fields are used only for comparison.
Use a new output folder; raw accepted fields and failures are retained.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import shutil
import sys
import time
from types import SimpleNamespace
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import fem_reference as fem
from joint_water_reference import Water
from verify_high_order import physical
from verify_joint_coupling import contact
from hydrostatic_cut_volume_reference import force as cut_force
from refine_tet10 import refine
from diagnose_wet_joint_convergence import configuration, schedule
from diagnose_cut_volume_history import FIELDS
from wet_joint_cycle_reference import serializable
from audit_wet_joint_accepted_history import history as codec
spec = importlib.util.spec_from_file_location('convergence_reference_probe', HERE/'reference_probe.py')
probe = importlib.util.module_from_spec(spec); spec.loader.exec_module(probe)
BASE = 'bc6eb1c39749559db91e20d2a4cba2d20679649c'
STAGES = ('geostatic', 'shear', 'open', 'recontact', 'reverse', 'reload', 'return')


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    with path.open('x', encoding='utf-8') as f:
        json.dump(serializable(value), f, indent=2, ensure_ascii=False, allow_nan=False); f.write('\n')


class Writer:
    def __init__(self, folder):
        folder.mkdir(parents=True, exist_ok=False)
        self.folder = folder; self.rows = []; self.values = []; self.chunks = []; self.count = 0

    def __call__(self, row, fields):
        self.rows.append(row); self.values.append({k: np.asarray(fields[k]).copy() for k in FIELDS})
        if len(self.rows) == 32:
            self.flush()

    def flush(self):
        if not self.rows:
            return
        path = self.folder/f'{len(self.chunks):04d}.npz'
        arrays = {k: np.stack([x[k] for x in self.values]) for k in FIELDS}
        with path.open('xb') as f:
            np.savez_compressed(f, **{k:codec.encode(v) for k,v in arrays.items()})
        self.chunks.append(dict(file=path.name, sha256=sha(path), rows=self.rows, start=self.count, count=len(self.rows)))
        self.count += len(self.rows); self.rows = []; self.values = []

    def close(self):
        self.flush(); write(self.folder/'manifest.json', dict(count=self.count, chunks=self.chunks, encoding=codec.ENCODING))


def states(folder):
    manifest = read(folder/'manifest.json'); count = 0
    for chunk in manifest['chunks']:
        path = folder/chunk['file']
        assert sha(path) == chunk['sha256'], 'changed accepted field'
        assert chunk['start'] == count and chunk['count'] == len(chunk['rows'])
        with np.load(path, allow_pickle=False) as data:
            assert set(data.files) == set(FIELDS)
            decoded = {k:codec.decode(data[k],manifest['encoding']) if manifest.get('encoding') else data[k] for k in FIELDS}
            for i, row in enumerate(chunk['rows']):
                value = {k:decoded[k][i].copy() for k in FIELDS}
                assert all(np.isfinite(v).all() for v in value.values())
                yield row, value
                count += 1
    assert count == manifest['count'], 'incomplete accepted history'


def operators(mesh, ground, case, rule):
    water = Water(mesh, ground, case); s = water.contact
    xyz = np.asarray(mesh['nodes']); nd = xyz.size; k = np.zeros((nd, nd)); seed_force = np.zeros(nd)
    ops = []; regions = {r['id']: r['body'] for r in mesh['regions']}
    mats = {m['id']: m['mohr_coulomb'] for m in ground['materials']}
    cfg = ground['initialization']; level = cfg['water_level_m']
    for ei, mid in enumerate(mesh['material_ids']):
        q = fem.volume(mesh, ei, independent=rule == 'cut', independent_order=7)
        ids = (3*q['nodes'][:, None]+np.arange(3)).ravel(); d = fem.elastic_tensor(mats[mid]['young_kpa'], mats[mid]['poisson'])
        initial = cfg['geostatic']['seed']['states'][regions[mesh['region_ids'][ei]]]
        s0 = np.asarray(initial['effective_stress_compression_positive_kpa'])
        k[np.ix_(ids, ids)] += np.einsum('qai,ab,qbj,q->ij', q['b'], d, q['b'], q['weight'], optimize=True)
        seed_force[ids] -= np.einsum('qai,a,q->i', q['b'], s0, q['weight'])
        pore = water.gamma*np.maximum(level-q['position'][:, 2], 0)
        seed_force[ids] += (cut_force(mesh, ei, level, water.gamma)['force'].ravel() if rule == 'cut'
                           else -np.einsum('qai,q,q->i', q['b'][:, :3], pore, q['weight']))
        ops.append((fem.volume(mesh, ei), ids))
    active = np.ones(len(mesh['tetrahedra']), dtype=bool)
    interface = water.hydrostatic(active, level)
    base, _ = physical(mesh, ground, active, cfg['surface_surcharge_kpa'], case['pressure_on_original_surface_only'], joint_system=s)
    base = (base+interface['load']).ravel(); boundary = fem.conditions_load(mesh, case['boundary'], active).ravel()
    basis, err = fem.staged_constraint_space(mesh, case, dict(active_elements=active, boundary_factor=0., displacement_m=np.zeros_like(xyz)))
    assert err < 1e-10
    def reduce(v):
        return np.concatenate([b.T@v.reshape(-1, 3)[n] for n, b in basis.items()])
    return SimpleNamespace(water=water, system=s, k=k, seed=seed_force, base=base, boundary=boundary,
                           holding=seed_force-base, ops=ops, interface=interface, reduce=reduce, level=level)


def native_fields(state, op):
    j = state['joints']['points']; u = np.asarray(state['displacement_m']); m = state['material_states']; s = op.system
    traction = np.asarray([x['traction_kpa'] for x in j]); jump = np.asarray([x['jump_m'] for x in j])
    qn = (traction*s.normal).sum(axis=1); tau = traction-qn[:, None]*s.normal
    strain = np.concatenate([np.einsum('qai,i->qa', q['b'], u.ravel()[ids]) for q, ids in op.ops])
    level = op.level
    pore = np.concatenate([op.water.gamma*np.maximum(level-q['position'][:, 2], 0) for q, _ in op.ops])
    internal = op.k@u.ravel()+op.seed+np.einsum('qag,qa,q->g', s.c, traction, s.weight)
    return dict(u=u, stress=np.asarray([x['effective_stress_compression_positive_kpa'] for x in m]),
                strain=strain, void=np.asarray([x['void_ratio'] for x in m]), pore=pore,
                plastic_strain=np.asarray([x['plastic_strain_compression_positive'] for x in m]), jump=jump, traction=traction,
                plastic_slip=np.asarray([x['plastic_slip_m'] for x in j]), accumulated=np.asarray([x['accumulated_plastic_slip_m'] for x in j]),
                elastic_energy=qn**2/(2*s.kn)+(tau*tau).sum(axis=1)/(2*s.kt),
                plastic_dissipation=np.asarray([x['plastic_dissipation_kpa_m'] for x in j]),
                numerical_dissipation=np.asarray([x['numerical_dissipation_kpa_m'] for x in j]),
                endpoint_work=np.asarray([x['endpoint_work_kpa_m'] for x in j]),
                physical_pressure=np.asarray(state['joint_water']['physical_pressure_kpa']),
                water_load=np.asarray(state['joint_water']['nodal_load_kn']).reshape(-1, 3),
                reaction=(internal-np.asarray(state['applied_nodal_load_kn'])).reshape(-1, 3), gap=(jump*s.normal).sum(axis=1))


def audit(folder, mesh, ground, case, rule, steps):
    op = operators(mesh, ground, case, rule); s = op.system
    count = 0; residual = 0.; pressure_delta = 0.; pore_delta = 0.; stress_delta = 0.; prev = None
    for row, f in states(folder):
        phase, step = divmod(count, steps)
        assert row['stage'] == STAGES[phase] and row['progress'] == (step+1)/steps, 'missing or shifted history'
        origin = 0. if phase <= 1 else case['stages'][phase-2]['target_factor']
        expected_factor = 0. if phase == 0 else origin+(case['stages'][phase-1]['target_factor']-origin)*row['progress']
        assert row['factor'] == expected_factor, 'input boundary-factor schedule differs'
        count += 1
        np.testing.assert_allclose(f['physical_pressure'], op.interface['pressure'], atol=1e-8, rtol=0)
        np.testing.assert_allclose(f['water_load'], op.interface['load'], atol=1e-8, rtol=0)
        for condition in case['boundary']['constraints']:
            nodes = mesh['node_sets'][condition['node_set']]
            xyz = np.asarray(mesh['nodes'])[nodes]
            expected = expected_factor*(condition['value_m']+(xyz-condition.get('reference_point_m',[0.,0.,0.]))@np.asarray(condition.get('gradient',[0.,0.,0.])))
            np.testing.assert_allclose(f['u'][nodes]@np.asarray(condition['direction']),expected,atol=1e-12,rtol=0)
        u = f['u'].ravel(); q = f['traction']; gap = (f['jump']*s.normal).sum(axis=1); qn = (q*s.normal).sum(axis=1)
        np.testing.assert_allclose(f['jump'], np.einsum('qag,g->qa', s.c, u), atol=1e-10, rtol=0)
        np.testing.assert_allclose(qn, s.kn*np.minimum(gap, 0), atol=1e-7, rtol=1e-9)
        tau = q-qn[:, None]*s.normal
        assert np.all(np.linalg.norm(tau, axis=1) <= s.mu*np.maximum(-qn, 0)+1e-7)
        np.testing.assert_allclose(f['endpoint_work'], f['elastic_energy']+f['plastic_dissipation']+f['numerical_dissipation'], atol=1e-9, rtol=1e-8)
        assert not f['plastic_strain'].any(), 'elastic bulk branch only'
        assert np.min(f['accumulated']) >= 0 and np.all(f['accumulated']+1e-10 >= np.linalg.norm(f['plastic_slip'], axis=1))
        if prev is not None:
            assert np.all(f['accumulated']+1e-12 >= prev['accumulated'])
            pore_delta = max(pore_delta, float(np.max(abs(f['pore']-prev['pore']))))
            pressure_delta = max(pressure_delta, float(np.max(abs(f['physical_pressure']-prev['physical_pressure']))))
        force = op.k@u+op.seed+np.einsum('qag,qa,q->g', s.c, q, s.weight)
        load = op.base+(1-row['progress'])*op.holding if row['stage'] == 'geostatic' else op.base+row['factor']*op.boundary
        residual = max(residual, float(np.linalg.norm(op.reduce(force-load))/max(np.linalg.norm(load), 1.)))
        start = 0
        for ei, ((volume, ids), mid) in enumerate(zip(op.ops, mesh['material_ids'])):
            mat = next(m['mohr_coulomb'] for m in ground['materials'] if m['id'] == mid)
            d = fem.elastic_tensor(mat['young_kpa'], mat['poisson']); n = len(volume['weight'])
            # Both bodies have the same input isotropic effective seed; no history is injected.
            body = next(r['body'] for r in mesh['regions'] if r['id'] == mesh['region_ids'][ei])
            seed = np.asarray(ground['initialization']['geostatic']['seed']['states'][body]['effective_stress_compression_positive_kpa'])
            expected_stress = seed-f['strain'][start:start+n]@d.T
            np.testing.assert_allclose(f['stress'][start:start+n],expected_stress,atol=1e-7,rtol=1e-9)
            np.testing.assert_allclose(f['strain'][start:start+n],np.einsum('qai,i->qa',volume['b'],u[ids]),atol=1e-10,rtol=0)
            initial=ground['initialization']['geostatic']['seed']['states'][body]
            np.testing.assert_allclose(f['void'][start:start+n],(1+initial['void_ratio'])*np.exp(f['strain'][start:start+n,:3].sum(axis=1))-1,atol=1e-8,rtol=0)
            np.testing.assert_allclose(f['pore'][start:start+n],op.water.gamma*np.maximum(op.level-volume['position'][:,2],0),atol=1e-8,rtol=0)
            assert np.max(probe.normalized_yield(mat,f['stress'][start:start+n]))<=1e-8
            stress_delta = max(stress_delta, float(np.max(abs(f['stress'][start:start+n]-expected_stress))))
            start += n
        prev = f
    assert count == 7*steps, 'incomplete path'
    native = folder.name == 'native'
    tolerance = 1.05e-8+1e-10 if native else 1e-10
    assert residual < tolerance, ('existing force audit', residual, tolerance)
    return dict(accepted_states=count, maximum_reassembled_force_residual=residual, existing_force_audit_tolerance=tolerance,
                maximum_pressure_change_kpa=pressure_delta, maximum_pore_change_kpa=pore_delta,
                maximum_elastic_stress_compatibility_error_kpa=stress_delta,
                all_saved_states_audited=True, numerical_acceptance=False)


def route_list():
    rows = []
    for g in (16, 32, 64, 128):
        rows.append(dict(id=f'h0-g{g}-l8-cut', h=0, g=g, local=8, rule='cut', maximum=None, native=False))
    for g, local in ((64,16),(64,32),(32,4),(128,16)):
        rows.append(dict(id=f'h0-g{g}-l{local}-cut', h=0,g=g,local=local,rule='cut',maximum=None,native=False))
    rows.append(dict(id='h0-g64-m5e-05-discrete',h=0,g=64,local=1,rule='discrete',maximum=5e-5,native=False))
    for h, g, maximum in ((0,32,5e-5),(0,64,5e-5),(0,128,5e-5),(0,64,1e-4),(0,64,2.5e-5),(1,32,5e-5),(1,64,5e-5)):
        rows.append(dict(id=f'h{h}-g{g}-m{maximum:g}-cut',h=h,g=g,local=1,rule='cut',maximum=maximum,native=True))
    return rows


def solve_route(output, route, solver):
    name = route['id']; folder = output/name; folder.mkdir()
    mesh = read(output/f'mesh-h{route["h"]}.json'); ground = read(output/'ground.json'); case = read(output/'stages.json')
    write(folder/'configuration.json', route)
    writer = Writer(folder/'independent'); started = time.monotonic()
    try:
        value = probe._reference(mesh, ground, case, route['g'], route['local'], 7 if route['rule']=='cut' else mesh['quadrature_order'],
                                 route['rule']=='discrete', route['rule']=='cut', writer,
                                 maximum_m=route['maximum'], maximum_dofs=15000 if route.get('user_pc') else 3000)
        writer.close(); write(folder/'independent-result.json', value)
        check = audit(folder/'independent', mesh, ground, case, route['rule'], route['g'])
        write(folder/'independent-audit.json', check)
    except Exception as error:
        write(folder/'independent-failure.json',dict(passed=False, error=repr(error), elapsed_seconds=time.monotonic()-started)); raise
    write(folder/'independent-timing.json',dict(elapsed_seconds=time.monotonic()-started))
    print('independent complete', name, flush=True)
    if route['native']:
        run_native(output, route, solver)


def run_native(output, route, solver):
    name=route['id']; folder=output/name
    mesh=read(output/f'mesh-h{route["h"]}.json'); ground=read(output/'ground.json'); case=read(output/'stages.json')
    if (folder/'native-raw/verification.json').exists():
        convert_native(output, route)
        write(folder/'native-reuse.json',dict(complete=True, raw_verification_sha256=sha(folder/'native-raw/verification.json')))
        print('native raw re-audited', name, flush=True)
        return
    original = configuration(case, route['g'], route['maximum'])
    write(folder/'native-stages.json', original)
    command = [str(solver), str(output/f'mesh-h{route["h"]}.json'), str(output/'ground.json'), str(folder/'native-stages.json'), str(folder/'native-raw')]
    if route.get('user_pc'):
        command.append('--user-pc')
    started = time.monotonic()
    with (folder/'native.log').open('x') as log:
        process = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=route.get('native_timeout_seconds',1800))
    assert process.returncode == 0, ('native solver incomplete', name)
    convert_native(output, route)
    write(folder/'native-timing.json',dict(elapsed_seconds=time.monotonic()-started, exit_code=0, command=command))
    print('native complete', name, flush=True)


def convert_native(output, route):
    name=route['id'];folder=output/name
    mesh=read(output/f'mesh-h{route["h"]}.json');ground=read(output/'ground.json');case=read(output/'stages.json')
    verification=read(folder/'native-raw/verification.json')
    assert verification['complete'] and verification['uninterrupted_run_exact_except_elapsed']
    assert verification['accepted_states']==7*route['g']
    run = read(folder/'native-raw/run.json'); schedule(run, route['g'])
    assert fem.normalize_ground(run['ground']) == ground, 'native ground input differs'
    assert fem.normalize_stage_case(run['case']) == fem.normalize_stage_case(configuration(case,route['g'],route['maximum'])), 'native case input differs'
    op = operators(mesh, ground, case, 'cut'); writer = Writer(folder/'native')
    with (folder/'native-raw/accepted.jsonl').open(encoding='utf-8') as f:
        for i, line in enumerate(f):
            state = json.loads(line); phase, step = divmod(i, route['g'])
            row = dict(stage=STAGES[phase], progress=(step+1)/route['g'], factor=state['boundary_factor'])
            assert state['stage_index'] == phase
            writer(row, native_fields(state, op))
    writer.close()
    write(folder/'native-audit.json', audit(folder/'native', mesh, ground, case, 'cut', route['g']))


def freeze_solver(output, solver):
    """Pin executable identity even if another build updates target/ during a study."""
    snapshot=output/('native-probe.exe' if solver.suffix.lower()=='.exe' else 'native-probe')
    metadata=output/'native-probe-sha256.json'
    if not snapshot.exists():
        shutil.copy2(solver,snapshot)
        write(metadata,dict(sha256=sha(snapshot),source_executable=str(solver),binary_archived=False))
    assert sha(snapshot)==read(metadata)['sha256'], 'pinned native executable changed'
    return snapshot


def run(output, solver, only=None, native_only=False):
    if not output.exists():
        output.mkdir(parents=True)
        values = [read(ROOT/'examples'/n) for n in ('joint-water-column-tet10-mesh.json','ground-joint-water-mc.json','stages-joint-water-cycle.json')]
        mesh, ground, case = values; refined, refinement = refine(mesh)
        write(output/'mesh-h0.json',mesh); write(output/'mesh-h1.json',refined); write(output/'refinement.json',refinement)
        write(output/'ground.json',fem.normalize_ground(ground)); write(output/'stages.json',fem.normalize_stage_case(case))
        write(output/'routes.json',route_list())
    solver=freeze_solver(output,solver)
    existing = read(output/'routes.json')
    if existing != route_list():
        raise ValueError('route inventory differs: start a new output or explicitly record an extension')
    failures = {}
    for route in route_list():
        if only is None or route['id'] in only:
            try:
                if native_only:
                    run_native(output, route, solver)
                else:
                    solve_route(output, route, solver)
            except Exception as error:
                path = output/route['id']/'failure.json'
                failures[route['id']] = dict(passed=False, error=repr(error), completed=False)
                write(path, failures[route['id']])
                print('RETAINED FAILURE', route['id'], repr(error), flush=True)
    return failures


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path,required=True)
    p.add_argument('--solver',type=Path,default=ROOT/'target/release/examples/wet_joint_history_probe')
    p.add_argument('--only',nargs='+'); p.add_argument('--native-only',action='store_true'); a=p.parse_args(); failures=run(a.output.resolve(),a.solver.resolve(),a.only,a.native_only); sys.exit(bool(failures))
