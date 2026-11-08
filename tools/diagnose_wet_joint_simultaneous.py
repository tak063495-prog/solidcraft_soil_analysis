"""Existing small wet JOINT: finer simultaneous paths and stored-force audit.

Input-only solves precede all historical result reads. This diagnostic never
certifies full-history convergence or installs a production volume rule.
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
from diagnose_cut_volume_history import inputs, validate_path, summary, write
from diagnose_wet_joint_cut_volume import components, error
from diagnose_wet_joint_convergence import schedule
from joint_water_reference import Water
from verify_high_order import prepare, internal, physical
from verify_joint_water import audit
from wet_joint_cycle_reference import cut_volume_reference, differences

ROOT = Path(__file__).resolve().parents[1]
BASE = '9d3374bd2573f1931a4206383964b527d26f8c62'
HISTORY = ROOT/'interim/2026-11-08-nc-midpoint-resume-and-cut-history'
NATIVE = ROOT/'interim/2026-11-08-wet-joint-separated-convergence'
PAIRS = (('g512-l8', 'g1024-l32'), ('g512-l32', 'g1024-l32'),
         ('g1024-l8', 'g1024-l32'), ('g1024-l8', 'g2048-l32'),
         ('g1024-l32', 'g2048-l32'))


def recorded(folder, names):
    manifest = json.loads((folder/'archive-sha256.json').read_text())
    archive = folder/'evidence.zip'
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    assert digest == manifest['archive_sha256'], 'archive changed'
    values = {}; entries = {}
    with zipfile.ZipFile(archive) as stream:
        for name in names:
            raw = stream.read(name); h = hashlib.sha256(raw).hexdigest()
            assert h == manifest['entries'][name], ('entry changed', name)
            values[name] = json.loads(raw); entries[name] = h
    return values, dict(archive=str(archive.relative_to(ROOT)), archive_sha256=digest,
                        entries=entries, rerun=False)


def force_alignment(mesh, ground, case, run):
    """Audit saved native fields; no independent mechanical solution claimed."""
    expected_case = fem.normalize_stage_case(run['case'])
    # Numerical controls differ, while geometry/material/load/water scope is bound.
    from diagnose_wet_joint_convergence import configuration
    assert expected_case == fem.normalize_stage_case(configuration(case, 256, 5e-6))
    assert run['ground'] == ground
    schedules = schedule(run, 256)
    original = audit(mesh, ground, run)
    sampled = components(mesh, ground, 'gj', 3)
    exact = components(mesh, ground, 'cut', 7)
    delta = (exact['pore_force']-sampled['pore_force']).reshape(-1, 3)
    q, offset = prepare(mesh, run['integration_points'])
    water = Water(mesh, ground, run['case'], run['joint_integration_points'])
    rows = []
    for frame in run['frames']:
        state = frame['state']; active = np.asarray(state['active_elements'])
        assert active.all() and state.get('pore_pressure_kpa') is None
        force, _ = internal(mesh, ground, state['material_states'], active, q, offset)
        interface = water.hydrostatic(active, ground['initialization']['water_level_m'])
        force += water.audit(state, frame['joint_results'], interface)
        load, _ = physical(mesh, ground, active, state['surface_pressure_kpa'],
                            run['case']['pressure_on_original_surface_only'], joint_system=water.contact)
        load += interface['load']
        load += state['boundary_factor']*fem.conditions_load(mesh, run['case']['boundary'], active)
        basis, constraint_error = fem.staged_constraint_space(mesh, run['case'], state)
        assert constraint_error < 1e-10
        def reduced(value):
            return np.concatenate([b.T@value[n] for n, b in basis.items()])
        scale = max(np.linalg.norm(load), 1.)
        old = float(np.linalg.norm(reduced(force-load))/scale)
        new = float(np.linalg.norm(reduced(force+delta-load))/scale)
        assert old < 1.05e-8+1e-10, 'saved production equilibrium failed'
        rows.append(dict(frame=frame['name'], sampled_relative_force_residual=old,
                         exact_volume_relative_force_residual=new,
                         exact_volume_meets_existing_force_tolerance=new < 1e-8,
                         physical_pressure=error(state['joint_water']['physical_pressure_kpa'], interface['pressure']),
                         interface_load=error(np.asarray(state['joint_water']['nodal_load_kn']).reshape(-1, 3), interface['load'])))
    return dict(stored_field_audit=original, schedules=schedules, frames=rows,
                volume_force_difference=error(sampled['pore_force'], exact['pore_force']),
                volume_force_correction_kn=delta, elements=sampled['elements'],
                exact_volume_equilibrium_accepted=all(r['exact_volume_meets_existing_force_tolerance'] for r in rows),
                native_result_modified=False, independent_plastic_solution=False)


def compare_paths(paths):
    result = {}
    for a, b in PAIRS:
        rows = differences(paths[a], paths[b])
        result[a+'__'+b] = dict(fields=rows, maxima=summary(rows))
    return result


def measure(folder):
    mesh, ground, case = inputs()
    for name, expected in zip(('mesh', 'ground', 'stages'), (mesh, ground, case)):
        assert json.loads((folder/(name+'.json')).read_text()) == expected, 'changed input'
    paths = {}; checks = {}; failures = {}
    for g, local in ((1024, 32), (2048, 32)):
        key = f'g{g}-l{local}'
        value = json.loads((folder/(key+'.json')).read_text())
        if 'failure' in value:
            failures[key] = value
        else:
            checks[key] = validate_path(value, mesh, ground, case, g, local); paths[key] = value
    old, history_provenance = recorded(HISTORY, [f'cut-history/{key}.json' for key in ('g512-l8', 'g512-l32', 'g1024-l8')])
    for name, value in old.items():
        key = Path(name).stem; g, local = [int(x[1:]) for x in key.split('-')]
        validate_path(value, mesh, ground, case, g, local); paths[key] = value
    names = ['separated/g256-m5e-06/'+n+'.json' for n in ('mesh', 'ground', 'stages', 'stage-results')]
    old_native, native_provenance = recorded(NATIVE, names)
    native_mesh, native_ground, native_case, run = [old_native[n] for n in names]
    assert native_mesh == mesh and fem.normalize_ground(native_ground) == ground
    assert fem.normalize_stage_case(native_case) == run['case']
    alignment = force_alignment(mesh, ground, case, run)
    protected = ROOT/'interim/2026-11-08-wet-joint-separated-convergence/measurements.json'
    return dict(schema='soil-fem-wet-joint-simultaneous/1', starting_commit=BASE,
                passed=False, history_converged=False, production_adopted=False,
                diagnostic_completed=not failures, reference_uses_native_results=False,
                native_rerun=False, nc_rerun=False, failures=failures, path_field_checks=checks,
                timings=json.loads((folder/'timings.json').read_text()),
                comparisons=compare_paths(paths) if not failures else {}, force_alignment=alignment,
                baselines=dict(history=history_provenance, native=native_provenance),
                protected_old_measurements=dict(path=str(protected.relative_to(ROOT)),
                    sha256=hashlib.sha256(protected.read_bytes()).hexdigest(), rewritten=False,
                    accumulated_128_to_256_relative_l2=0.01703832492167064,
                    plastic_slip_128_to_256_relative_l2=0.22431468120870923),
                full_history_acceptance=dict(accepted=False, reasons=[
                    'No existing all-history error criterion certifies measured path differences.',
                    'Native volume force still differs from analytic partially wet volume force.',
                    'New references store all accepted progress/residual rows and seven full fields; interior full fields are not stored.',
                    'Refined analytic-volume production path and restart through that path are absent.',
                    '2048/64 and subsequent simultaneous refinement remain a user-PC comparison.']),
                windows_ci='unverified', p3b_p4_complete=False,
                structure_water_joint_complete=False, p6_complete=False)


def run(output):
    output.mkdir(parents=True, exist_ok=False)
    values = inputs()
    for name, value in zip(('mesh', 'ground', 'stages'), values):
        write(output/(name+'.json'), value)
    timings = []
    for steps, local in ((1024, 32), (2048, 32)):
        key = f'g{steps}-l{local}'; start = time.monotonic()
        try:
            value = cut_volume_reference(*values, steps, local, 7)
            validate_path(value, *values, steps, local)
        except (AssertionError, np.linalg.LinAlgError) as error_value:
            value = dict(failure=str(error_value), passed=False, steps=steps, contact_steps=local)
        timings.append(dict(route=key, elapsed_seconds=time.monotonic()-start, complete='failure' not in value))
        write(output/(key+'.json'), value)
        print(key, timings[-1], flush=True)
    write(output/'timings.json', timings)
    report = measure(output); write(output/'measurements.json', report)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--remeasure', type=Path)
    args = parser.parse_args()
    if args.remeasure:
        result = measure(args.remeasure.resolve()); write(args.output, result)
    else:
        result = run(args.output.resolve())
    print(json.dumps({k: result[k] for k in ('passed', 'diagnostic_completed', 'history_converged')}))
