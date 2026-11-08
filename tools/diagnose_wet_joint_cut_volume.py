"""Isolate the hydrostatic-volume kink of the existing small wet JOINT.

Only the new cut-volume reference is solved. Recorded GJ3/GL7/native 256-step
histories are read with SHA checks, never rerun or rewritten. This diagnostic
does not certify contact-history convergence or adopt a production change.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
import hashlib
import json
from pathlib import Path
import time
import zipfile
import numpy as np
import fem_reference as fem
from hydrostatic_cut_volume_reference import force
from joint_water_reference import Water
from wet_joint_cycle_reference import cut_volume_reference, differences, serializable

ROOT = Path(__file__).resolve().parents[1]
BASE = '8397d04b569ce264bb939cad731ff373a0d930a3'
PRIOR = ROOT/'interim/2026-11-08-wet-joint-separated-convergence'
INPUTS = ('joint-water-column-tet10-mesh', 'ground-joint-water-mc', 'stages-joint-water-cycle')


def write(path, value):
    path.write_text(json.dumps(serializable(value), ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def error(a, b):
    a = np.asarray(a); b = np.asarray(b)
    assert a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all()
    return dict(maximum_absolute=float(np.max(abs(a-b), initial=0)),
                relative_l2=float(np.linalg.norm(a-b)/max(np.linalg.norm(b), 1e-12)))


def components(mesh, ground, rule, order=7):
    """No displacement/contact result enters this force decomposition."""
    xyz = np.asarray(mesh['nodes']); nd = xyz.size
    k = np.zeros((nd, nd)); effective = np.zeros(nd); pore = np.zeros(nd)
    rows = []; volume = 0.; pressure = 0.; pressure_first = np.zeros(3)
    materials = {m['id']: m['mohr_coulomb'] for m in ground['materials']}
    bodies = {r['id']: r['body'] for r in mesh['regions']}
    cfg = ground['initialization']; level = cfg['water_level_m']; gamma = cfg['water_unit_weight_kn_m3']
    for ei, mid in enumerate(mesh['material_ids']):
        q = fem.volume(mesh, ei, independent=rule != 'gj', independent_order=order)
        ids = (3*q['nodes'][:, None]+np.arange(3)).ravel()
        mat = materials[mid]; d = fem.elastic_tensor(mat['young_kpa'], mat['poisson'])
        s0 = cfg['geostatic']['seed']['states'][bodies[mesh['region_ids'][ei]]]['effective_stress_compression_positive_kpa']
        ke = np.einsum('qai,ab,qbj,q->ij', q['b'], d, q['b'], q['weight'], optimize=True)
        k[np.ix_(ids, ids)] += ke
        effective[ids] -= np.einsum('qai,a,q->i', q['b'], s0, q['weight'])
        exact = force(mesh, ei, level, gamma)
        if rule == 'cut':
            term = exact['force'].ravel()
        else:
            p = gamma*np.maximum(level-q['position'][:, 2], 0)
            term = -np.einsum('qi,q,q->i', q['b'][:, :3].sum(axis=1), p, q['weight'])
        pore[ids] += term
        volume += exact['wet_volume']; pressure += exact['pressure_integral']
        pressure_first += exact['pressure_first_moment']
        rows.append(dict(element=ei, partially_wet=bool(xyz[q['nodes'][:4], 2].min() < level < xyz[q['nodes'][:4], 2].max()),
                         force_error=error(term, exact['force'].ravel()),
                         force_kn=term, exact_force_kn=exact['force'], wet_volume_m3=exact['wet_volume']))
    return dict(stiffness=k, effective_seed=effective, pore_force=pore, elements=rows,
                wet_volume_m3=volume, pressure_integral_kpa_m3=pressure, pressure_first_moment=pressure_first)


def linear_initial_correction(mesh, ground, case, components_from, components_to):
    """Input-only closed-stick correction for the initial equilibrium only."""
    xyz = np.asarray(mesh['nodes']); water = Water(mesh, ground, case); s = water.contact
    basis, residual = fem.staged_constraint_space(mesh, case, dict(
        active_elements=[True]*len(mesh['tetrahedra']), boundary_factor=0., displacement_m=np.zeros_like(xyz)))
    assert residual == 0
    nf = sum(b.shape[1] for b in basis.values()); t = np.zeros((xyz.size, nf)); start = 0
    for node, b in basis.items():
        t[3*node:3*node+3, start:start+b.shape[1]] = b; start += b.shape[1]
    nn = np.einsum('qa,qb->qab', s.normal, s.normal)
    dc = s.kn[:, None, None]*nn+s.kt[:, None, None]*(np.eye(3)-nn)
    contact_k = np.einsum('qai,qab,qbj,q->ij', s.c, dc, s.c, s.weight, optimize=True)
    k = np.asarray(components_to['stiffness'])+contact_k
    delta = np.asarray(components_to['pore_force'])-components_from['pore_force']
    free = np.linalg.solve(t.T@k@t, -t.T@delta)
    return (t@free).reshape(-1, 3)


def recorded_inputs():
    manifest = json.loads((PRIOR/'archive-sha256.json').read_text())
    archive = PRIOR/'evidence.zip'
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == manifest['archive_sha256']
    paths = {'discrete': 'separated/g256-l8-discrete.json', 'gl7': 'separated/g256-l8-gl7.json'}
    values = {}; used = {}
    with zipfile.ZipFile(archive) as z:
        for key, path in paths.items():
            raw = z.read(path); digest = hashlib.sha256(raw).hexdigest()
            assert digest == manifest['entries'][path], ('baseline changed', path)
            values[key] = json.loads(raw); used[path] = digest
        for name, source in zip(('mesh', 'ground', 'stages'), INPUTS):
            path = 'separated/'+name+'.json'; raw = z.read(path)
            assert hashlib.sha256(raw).hexdigest() == manifest['entries'][path]
            assert json.loads(raw) == json.loads((ROOT/'examples'/(source+'.json')).read_text())
            used[path] = manifest['entries'][path]
    return values, dict(archive_sha256=manifest['archive_sha256'], entries=used, rerun=False)


def measure(mesh, ground, case, c, new, elapsed):
    """Recompute measurements from saved new evidence and immutable baselines."""
    originals = [json.loads((ROOT/'examples'/(name+'.json')).read_text()) for name in INPUTS]
    expected_inputs = (originals[0], fem.normalize_ground(originals[1]), fem.normalize_stage_case(originals[2]))
    assert (mesh, ground, case) == expected_inputs, 'existing input only'
    checks = []
    for ei in range(len(mesh['tetrahedra'])):
        exact = force(mesh, ei, .75, 10)
        checks.append(dict(element=ei, analytic_vs_cut_gl3=error(exact['force'], force(mesh, ei, .75, 10, 3)['force']),
                           analytic_vs_cut_gl5=error(exact['force'], force(mesh, ei, .75, 10, 5)['force'])))
    old, provenance = recorded_inputs()
    assert new['steps'] == 256 and new['contact_steps'] == 8 and len(new['frames']) == 7
    assert new['volume_rule'] == 'affine_cut_hydrostatic' and not new['reference_uses_native_results']
    assert len(new['path']) == 7*256 and new['maximum_relative_force_residual'] < 1e-10
    for name in ['geostatic']+[s['name'] for s in case['stages']]:
        rows = [p for p in new['path'] if p['stage'] == name]
        assert [p['progress'] for p in rows] == [i/256 for i in range(1, 257)], 'changed cut-reference schedule'
        assert all(p['relative_force_residual'] < 1e-10 for p in rows)
    assert all(r['steps'] == 256 and r['contact_steps'] == 8 for r in old.values())
    correction = {}
    for name, key in [('discrete', 'gj3'), ('gl7', 'gl7')]:
        # Validate the closed-stick hypothesis before using the linear ablation.
        for r in (old[name], new):
            f = r['frames'][0]; s = Water(mesh, ground, case).contact
            traction = np.asarray(f['traction']); normal = np.einsum('qa,qa->q', traction, s.normal)
            shear = traction-normal[:, None]*s.normal
            assert np.all(np.asarray(f['gap']) < 0)
            assert np.all(np.linalg.norm(shear, axis=1) < s.mu*(-normal))
        predicted = linear_initial_correction(mesh, ground, case, c[key], c['cut'])
        actual = np.asarray(new['frames'][0]['u'])-old[name]['frames'][0]['u']
        correction[name] = dict(predicted_delta_m=predicted, actual_delta_m=actual,
                                unexplained_difference=error(actual, predicted))
    report = dict(schema='soil-fem-wet-joint-cut-volume/1', starting_commit=BASE, passed=False,
                  diagnostic_completed=True, production_adopted=False, baseline=provenance,
                  new_solve=dict(global_steps=256, contact_steps=8, elapsed_seconds=elapsed,
                                 maximum_relative_force_residual=new['maximum_relative_force_residual']),
                  component_differences={key: {field: error(c[key][field], c['cut'][field])
                    for field in ('stiffness', 'effective_seed', 'pore_force')} for key in ('gj3', 'gl7', 'gl11')},
                  analytic_vs_clipped_quadrature=checks,
                  analytic_cube_check=dict(wet_volume_m3=c['cut']['wet_volume_m3'],
                    expected_wet_volume_m3=.75, pressure_integral_kpa_m3=c['cut']['pressure_integral_kpa_m3'],
                    expected_pressure_integral_kpa_m3=10*.75**2/2,
                    pressure_first_moment=c['cut']['pressure_first_moment'],
                    expected_pressure_first_moment=[10*.75**2/4, 10*.75**2/4, 10*.75**3/6]),
                  new_cut_vs_recorded={key: differences(new, value) for key, value in old.items()},
                  initial_closed_stick_pressure_correction=correction,
                  reference_uses_native_results=False, history_converged=False,
                  full_history_acceptance='unresolved: this changes volume integration only; recorded global-step history sensitivity remains',
                  nc_executed=False, windows_ci='unverified', p3b_p4_complete=False,
                  structure_water_joint_complete=False, p6_complete=False)
    return report


def run(output):
    output.mkdir(parents=True, exist_ok=False)
    for name, source in zip(('mesh', 'ground', 'stages'), INPUTS):
        (output/(name+'.json')).write_bytes((ROOT/'examples'/(source+'.json')).read_bytes())
    mesh, ground, case = [json.loads((output/(name+'.json')).read_text()) for name in ('mesh', 'ground', 'stages')]
    ground = fem.normalize_ground(ground); case = fem.normalize_stage_case(case)
    # Component integration does not run or read any solution.
    c = {key: components(mesh, ground, rule, order) for key, rule, order in
         [('gj3', 'gj', 3), ('gl7', 'gl', 7), ('gl11', 'gl', 11), ('cut', 'cut', 7)]}
    write(output/'components.json', c)
    start = time.monotonic()
    new = cut_volume_reference(mesh, ground, case, 256, 8, 7)
    elapsed = time.monotonic()-start
    write(output/'cut-g256-l8-gl7.json', new); write(output/'timing.json', dict(elapsed_seconds=elapsed))
    print('new cut-volume reference completed', flush=True)
    report = measure(mesh, ground, case, c, new, elapsed)  # Baselines are first read after this solve.
    write(output/'measurements.json', report)
    return report


def remeasure(folder):
    read = lambda name: json.loads((folder/(name+'.json')).read_text())
    mesh, ground, case = [read(name) for name in ('mesh', 'ground', 'stages')]
    return measure(mesh, fem.normalize_ground(ground), fem.normalize_stage_case(case),
                   read('components'), read('cut-g256-l8-gl7'), read('timing')['elapsed_seconds'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--remeasure', type=Path, help='read saved evidence only, write a new JSON file')
    args = parser.parse_args()
    if args.remeasure:
        assert not args.output.exists(), 'new output file required'
        result = remeasure(args.remeasure.resolve()); write(args.output, result)
    else:
        result = run(args.output.resolve())
    print(json.dumps(dict(passed=result['passed'], diagnostic_completed=result['diagnostic_completed'])))
