"""Measure path sensitivity with exact affine hydrostatic volume force.

Existing small wet JOINT only. No native solve or production change. Read the
immutable 256-step reference only AFTER solving the new input-only paths.
Always passed=false: there is no existing full-history acceptance certificate.
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
from joint_water_reference import Water
from wet_joint_cycle_reference import cut_volume_reference, differences, serializable

ROOT = Path(__file__).resolve().parents[1]
BASE = '561afe93ab4ea82e9f56535bf9444f4cb3038577'
PRIOR = ROOT/'interim/2026-11-08-wet-joint-cut-volume'
INPUTS = ('joint-water-column-tet10-mesh', 'ground-joint-water-mc', 'stages-joint-water-cycle')
FIELDS = dict(u=(54, 3), stress=(324, 6), strain=(324, 6), void=(324,), pore=(324,),
              plastic_strain=(324, 6), jump=(32, 3), traction=(32, 3), plastic_slip=(32, 3),
              accumulated=(32,), elastic_energy=(32,), plastic_dissipation=(32,),
              numerical_dissipation=(32,), endpoint_work=(32,), physical_pressure=(32,),
              water_load=(54, 3), reaction=(54, 3), gap=(32,))


def write(path, value):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(serializable(value), stream, ensure_ascii=False, indent=2, allow_nan=False)


def inputs():
    mesh, ground, case = [json.loads((ROOT/'examples'/(n+'.json')).read_text()) for n in INPUTS]
    return mesh, fem.normalize_ground(ground), fem.normalize_stage_case(case)


def validate_path(value, mesh, ground, case, steps, local):
    """Reject missing/changed paths and invalid full fields; not a history gate."""
    assert (mesh, ground, case) == inputs(), 'existing wet JOINT inputs only'
    assert value['steps'] == steps and value['contact_steps'] == local
    assert value['volume_rule'] == 'affine_cut_hydrostatic'
    assert value['volume_order'] == 7 and not value['reference_uses_native_results']
    names = ['geostatic']+[s['name'] for s in case['stages']]
    assert len(value['path']) == steps*len(names), 'missing path'
    expected_names = ['initial_equilibrated']+names[1:]
    assert [f['name'] for f in value['frames']] == expected_names, 'missing frames'
    for name in names:
        rows = [p for p in value['path'] if p['stage'] == name]
        assert [p['progress'] for p in rows] == [i/steps for i in range(1, steps+1)], 'changed schedule'
        assert all(np.isfinite(p['relative_force_residual']) and p['relative_force_residual'] < 1e-10 for p in rows)
    assert value['maximum_relative_force_residual'] < 1e-10
    s = Water(mesh, ground, case).contact
    statuses = []
    for frame in value['frames']:
        for key, shape in FIELDS.items():
            a = np.asarray(frame[key]); assert a.shape == shape and np.isfinite(a).all(), ('invalid full field', key)
        jump = np.asarray(frame['jump']); q = np.asarray(frame['traction'])
        gap = np.einsum('qa,qa->q', jump, s.normal); qn = np.einsum('qa,qa->q', q, s.normal)
        tau = q-qn[:, None]*s.normal
        np.testing.assert_allclose(np.einsum('qag,g->qa', s.c, np.asarray(frame['u']).ravel()), jump, atol=1e-10, rtol=0)
        np.testing.assert_allclose(frame['gap'], gap, atol=1e-10, rtol=0)
        np.testing.assert_allclose(qn, s.kn*np.minimum(gap, 0), atol=1e-7, rtol=1e-9)
        assert np.all(np.linalg.norm(tau, axis=1) <= s.mu*np.maximum(-qn, 0)+1e-7), 'Coulomb limit'
        energy = qn**2/(2*s.kn)+np.einsum('qa,qa->q', tau, tau)/(2*s.kt)
        np.testing.assert_allclose(frame['elastic_energy'], energy, atol=1e-9, rtol=1e-8)
        np.testing.assert_allclose(frame['endpoint_work'], energy+np.asarray(frame['plastic_dissipation'])+
                                   np.asarray(frame['numerical_dissipation']), atol=1e-9, rtol=1e-8)
        for key in ('plastic_dissipation', 'numerical_dissipation', 'accumulated'):
            assert np.min(frame[key]) >= -1e-12, ('negative history', key)
        assert np.all(np.asarray(frame['accumulated'])+1e-10 >= np.linalg.norm(frame['plastic_slip'], axis=1))
        assert not np.asarray(frame['plastic_strain']).any(), 'bulk plasticity outside scope'
        statuses.append(dict(frame=frame['name'], open_points=int((gap > 0).sum()),
                             accumulated_max_m=float(np.max(frame['accumulated']))))
    return statuses


def recorded_256():
    archive = PRIOR/'evidence.zip'; manifest = json.loads((PRIOR/'archive-sha256.json').read_text())
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == manifest['archive_sha256']
    entries = {}; values = []
    with zipfile.ZipFile(archive) as z:
        for name in ('mesh.json', 'ground.json', 'stages.json', 'cut-g256-l8-gl7.json'):
            raw = z.read(name); digest = hashlib.sha256(raw).hexdigest()
            assert digest == manifest['entries'][name], ('baseline changed', name)
            entries[name] = digest; values.append(json.loads(raw))
    mesh, ground, case = values[:3]
    assert (mesh, fem.normalize_ground(ground), fem.normalize_stage_case(case)) == inputs()
    validate_path(values[3], *inputs(), 256, 8)
    return values[3], dict(archive_sha256=manifest['archive_sha256'], entries=entries, rerun=False)


def summary(rows):
    fields = {}
    for key in FIELDS:
        fields[key] = {metric: max((dict(frame=r['frame'], value=r[key][metric]) for r in rows),
                                  key=lambda r: r['value']) for metric in ('maximum_absolute', 'relative_l2')}
    return fields


def measure(folder):
    read = lambda name: json.loads((folder/(name+'.json')).read_text())
    mesh, ground, case = [read(n) for n in ('mesh', 'ground', 'stages')]
    assert (mesh, ground, case) == inputs(), 'existing wet JOINT inputs only'
    baseline, provenance = recorded_256()
    paths = {'g256-l8': baseline}; checks = {}; failures = {}; timings = read('timings')
    for steps, local in ((128, 8), (512, 8), (1024, 8), (512, 32)):
        key = f'g{steps}-l{local}'; result = read(key)
        if 'failure' in result:
            failures[key] = result; continue
        checks[key] = validate_path(result, mesh, ground, case, steps, local); paths[key] = result
    comparisons = {}
    for a, b in (('g128-l8', 'g256-l8'), ('g256-l8', 'g512-l8'),
                 ('g512-l8', 'g1024-l8'), ('g512-l8', 'g512-l32')):
        if a in paths and b in paths:
            rows = differences(paths[a], paths[b]); comparisons[a+'__'+b] = dict(fields=rows, maxima=summary(rows))
    old_path = ROOT/'interim/2026-11-08-wet-joint-separated-convergence/measurements.json'
    return dict(schema='soil-fem-cut-volume-history-sensitivity/1', starting_commit=BASE, passed=False,
                diagnostic_completed=not failures, history_converged=False, production_adopted=False,
                reference_uses_native_results=False, native_rerun=False, baseline=provenance,
                timings=timings, failures=failures, path_field_checks=checks, comparisons=comparisons,
                protected_old_measurements=dict(path=str(old_path.relative_to(ROOT)),
                    sha256=hashlib.sha256(old_path.read_bytes()).hexdigest(), rewritten=False,
                    accumulated_128_to_256_relative_l2=0.01703832492167064,
                    plastic_slip_128_to_256_relative_l2=0.22431468120870923),
                limitation='Exact volume integration is fixed; global/contact path differences are measurements. '
                           'No existing full-history criterion certifies convergence; production still uses its original volume rule.',
                windows_ci='unverified', p3b_p4_complete=False, structure_water_joint_complete=False, p6_complete=False)


def run(output):
    output.mkdir(parents=True, exist_ok=False)
    mesh, ground, case = inputs()
    for name, value in zip(('mesh', 'ground', 'stages'), (mesh, ground, case)):
        write(output/(name+'.json'), value)
    timings = []
    for steps, local in ((128, 8), (512, 8), (1024, 8), (512, 32)):
        key = f'g{steps}-l{local}'; started = time.monotonic()
        try:
            value = cut_volume_reference(mesh, ground, case, steps, local, 7)
            validate_path(value, mesh, ground, case, steps, local)
        except (AssertionError, np.linalg.LinAlgError) as error:
            value = dict(failure=str(error), passed=False, steps=steps, contact_steps=local)
        timings.append(dict(route=key, elapsed_seconds=time.monotonic()-started, complete='failure' not in value))
        write(output/(key+'.json'), value); print(key, timings[-1], flush=True)
    write(output/'timings.json', timings)
    # Read historical fields only after all new input-only solves.
    report = measure(output); write(output/'measurements.json', report)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--remeasure', type=Path, help='saved evidence only; output must be a new JSON file')
    args = parser.parse_args()
    result = measure(args.remeasure.resolve()) if args.remeasure else run(args.output.resolve())
    if args.remeasure:
        write(args.output.resolve(), result)
    print(json.dumps({key: result[key] for key in ('passed', 'diagnostic_completed', 'history_converged')}))
