"""Record discretization differences; no production acceptance is awarded."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def read(p):
    return json.loads(p.read_text(encoding='utf-8-sig'))


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def metric(a, b):
    a, b = np.asarray(a), np.asarray(b)
    assert a.shape == b.shape
    return dict(relative_l2=float(np.linalg.norm(a-b)/max(np.linalg.norm(b), 1e-30)),
                maximum_absolute=float(abs(a-b).max()))


def fields(a, b):
    s, t = [np.asarray(x['material_state']).reshape(-1, 14) for x in (a, b)]
    return dict(displacement_m=metric(a['displacement_m'], b['displacement_m']),
        stress_kpa=metric(s[:, :6], t[:, :6]), pc_kpa=metric(np.exp(s[:, 6]), np.exp(t[:, 6])),
        void_ratio=metric(np.exp(s[:, 7])-1, np.exp(t[:, 7])-1), plastic_strain=metric(s[:, 8:], t[:, 8:]))


def assess(study, fractional, full, output):
    paths = [study/name for name in ('reference-25.json', 'reference-50.json', 'initial-25.json', 'initial-50.json', 'initial-100.json', 'research-checks.json')]
    a, b, c, d, e, checks = map(read, paths)
    for item in (a, b, c, d, e):
        assert item['schema'] == 'soil-fem-projected-volume-mcc-diagnostic/1' and item['diagnostic_only']
        for name, digest in item['input_and_source_sha256'].items():
            assert sha(Path(name)) == digest, name
    assert a['load_fraction'] == b['load_fraction'] < 1 and c['load_fraction'] == d['load_fraction'] == e['load_fraction'] == 1
    assert checks['passed'] and not checks['production_backend_accepted']
    for name, digest in checks['sha256'].items():
        assert sha(Path(name)) == digest, name
    independent_fraction, independent_full = read(fractional), read(full)
    assert independent_fraction['load_fraction'] == b['load_fraction']
    for name, digest in independent_fraction['input_and_source_sha256'].items():
        assert sha(Path(name)) == digest, name
    expected = {Path(name).name: digest for name, digest in independent_full['input_sha256'].items()}
    for name in ('mesh.json', 'ground.json'):
        assert sha(study/'inputs'/name) == expected[name], name
    # Construction settings differ, but geostatic geometry, source state and
    # zero initial construction-boundary factor are identical raw input bytes.
    assert independent_full['frames'][0]['name'] == 'initial_equilibrated'
    assert all(independent_full['frames'][0]['active_elements'])
    full_frame = independent_full['frames'][0]
    s = np.asarray(full_frame['material_state']).reshape(-1, 14)
    report = dict(schema='soil-fem-projected-volume-discretization-assessment/1',
        diagnostic_only=True, production_backend_accepted=False, goal_complete=False, p6_complete=False,
        source_and_input_bindings_verified=True,
        fractional_reference_steps=[25, 50], full_initialization_reference_steps=[25, 50, 100],
        operators_and_material_theory_passed=True,
        fractional=dict(load_fraction=b['load_fraction'], refinement=fields(a, b),
            versus_compatible_independent=fields(b, independent_fraction),
            projected_plastic_points=b['plastic_points'],
            independent_plastic_points=int(np.count_nonzero(np.linalg.norm(np.asarray(independent_fraction['material_state']).reshape(-1, 14)[:, 8:], axis=1)>1e-10))),
        full_initialization=dict(refinement_25_to_50=fields(c, d), refinement_50_to_100=fields(d, e),
            coarse_25_to_50_displacement_gate_passed=fields(c, d)['displacement_m']['maximum_absolute'] < 1e-7,
            versus_compatible_independent_100=fields(e, full_frame),
            projected_plastic_points=e['plastic_points'],
            independent_plastic_points=int(np.count_nonzero(np.linalg.norm(s[:, 8:], axis=1)>1e-10)),
            work_conjugate_force_residual=e['relative_work_conjugate_force_residual'],
            compatible_stress_force_residual=e['relative_compatible_stress_force_residual'],
            maximum_gp_mean_pressure_projection_difference_kpa=e['maximum_gp_pressure_projection_difference_kpa']),
        interpretation='P1 volume removes the observed early free-surface GP mode on this mesh, but changes the discrete force/stress definition; integration convergence is not mesh convergence or acceptance of production stages',
        remaining=['Explicit production formulation/schema and work-conjugate reaction/release definition',
                   'Spatial convergence against independently refined compatible FEM',
                   'Frozen compatibility offsets across local DEATH, without rewriting pc/e/plastic history',
                   'Plastic load/unload/excavation, rejection, restart, output and native integration'],
        sha256={str(p.resolve()): sha(p) for p in paths+[fractional, full, Path(__file__)]})
    assert report['fractional']['refinement']['displacement_m']['maximum_absolute'] < 1e-8
    assert report['full_initialization']['refinement_50_to_100']['displacement_m']['maximum_absolute'] < 1e-7
    with output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps({k: report[k] for k in ('fractional', 'full_initialization', 'production_backend_accepted')}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ('study', 'fractional', 'full', 'output'):
        p.add_argument(name, type=Path)
    a = p.parse_args()
    assess(*(getattr(a, k).resolve() for k in ('study', 'fractional', 'full', 'output')))
