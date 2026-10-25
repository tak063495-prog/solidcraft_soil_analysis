"""Full-field nonlinear sparse/dense reference comparison on identical inputs."""
import argparse
import json
from pathlib import Path
import mcc_projected_volume_diagnostic as projected

np = projected.np


def check_bindings(result):
    for name, digest in result['input_and_source_sha256'].items():
        assert projected.sha(Path(name)) == digest, name


def compare(dense_study, sparse_study, output):
    left_initial = projected.read(dense_study/'initial-100.json')
    left_cycle = projected.read(dense_study/'reference-cycle-100.json')
    right_initial = projected.read(sparse_study/'initial.json')
    right_cycle = projected.read(sparse_study/'cycle.json')
    for result in (left_initial, left_cycle, right_initial, right_cycle):
        check_bindings(result)
    assert left_cycle['reference_steps'] == right_cycle['reference_steps'] == 100
    assert left_initial['reference_steps'] == right_initial['reference_steps'] == 100
    assert right_cycle['reference_backend'] == 'mcc_projected_sparse_reference.SparseProjectedSystem'
    for name in ('mesh.json', 'ground.json', 'case.json', 'mesh-identity.json'):
        assert projected.sha(dense_study/'inputs'/name) == projected.sha(sparse_study/'inputs'/name), name
    checks = []

    def field(name, a, b, tolerance):
        a, b = np.asarray(a), np.asarray(b)
        assert a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all()
        maximum = float(abs(a-b).max()) if a.size else 0.
        checks.append(dict(field=name, maximum_absolute_difference=maximum, tolerance=tolerance, passed=maximum < tolerance))

    def state(name, a, b):
        a, b = np.asarray(a), np.asarray(b)
        field(name+'/stress_kpa', a[..., :6], b[..., :6], 1e-6)
        field(name+'/preconsolidation_kpa', np.exp(a[..., 6]), np.exp(b[..., 6]), 1e-6)
        field(name+'/void_ratio', np.exp(a[..., 7])-1, np.exp(b[..., 7])-1, 1e-8)
        field(name+'/plastic_strain', a[..., 8:], b[..., 8:], 1e-8)

    field('initial/displacement_m', left_initial['displacement_m'], right_initial['displacement_m'], 1e-9)
    state('initial', left_initial['material_state'], right_initial['material_state'])
    assert len(left_cycle['frames']) == len(right_cycle['frames']) == 4
    for a, b in zip(left_cycle['frames'], right_cycle['frames']):
        name = a['name']
        assert name == b['name'] and a['active_elements'] == b['active_elements']
        field(name+'/displacement_m', a['displacement_m'], b['displacement_m'], 1e-9)
        state(name, a['material_state'], b['material_state'])
        for key, value, tolerance in (('pressure_dofs', 'compression_positive_kpa', 1e-6),
                                      ('volume_dofs', 'compression_positive', 1e-8),
                                      ('compatibility_offsets', 'integral_m3', 1e-10)):
            assert [(x['node_index'], x['material_id']) for x in a[key]] == [(x['node_index'], x['material_id']) for x in b[key]]
            field(name+'/'+key, [x[value] for x in a[key]], [x[value] for x in b[key]], tolerance)
        assert a['plastic_points'] == b['plastic_points']
        assert b['relative_work_conjugate_force_residual'] < 1e-8
        assert b['relative_frozen_volume_constraint_residual'] < 1e-8
        assert b['sparse_solver_diagnostics']['maximum_true_linear_residual'] < 1e-10
    assert len(left_cycle['release_audits']) == len(right_cycle['release_audits']) == 1
    for a, b in zip(left_cycle['release_audits'], right_cycle['release_audits']):
        assert a['removed_elements'] == b['removed_elements']
        assert a['material_history_unchanged_at_domain_restriction'] and b['material_history_unchanged_at_domain_restriction']
        field('excavation/holding_nodal_load_kn', a['holding_nodal_load_kn'], b['holding_nodal_load_kn'], 1e-6)
    passed = all(c['passed'] for c in checks)
    report = dict(schema='soil-fem-projected-nonlinear-reference-equivalence/1', passed=passed,
                  goal_complete=False, p6_complete=False, spatial_convergence_accepted=False,
                  full_field_checks=checks, sparse_solver=[f['sparse_solver_diagnostics'] for f in right_cycle['frames']],
                  sha256={str(p): projected.sha(p) for p in (
                      dense_study/'initial-100.json', dense_study/'reference-cycle-100.json',
                      sparse_study/'initial.json', sparse_study/'cycle.json', Path(__file__))},
                  limitations='Same-mesh reference backend equivalence; general terrain and spatial convergence remain unaccepted')
    with output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps(dict(passed=passed, checks=len(checks), failed=[c for c in checks if not c['passed']])), flush=True)
    assert passed


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ('dense', 'sparse', 'output'):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    compare(args.dense.resolve(), args.sparse.resolve(), args.output.resolve())
