"""Compare the full sparse compatible cycle to the sealed accepted PCG cycle."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def physical_case(case):
    result = {k: v for k, v in case.items() if k not in ('schema', 'settings', 'verification_reject_once_stages')}
    result['kinematics'] = case.get('settings', {}).get('kinematics', 'side_rollers')
    result['stages'] = [{k: v for k, v in s.items() if k != 'initial_increment'} for s in case['stages']]
    return result


def compare(root, accepted, study, output):
    manifest_path = accepted/'final-validation.json'
    manifest = read(manifest_path)
    assert manifest['passed']
    reference_path = accepted/'data/slope-cycle-reference-100.json'
    reference = read(reference_path)
    reference_key = reference_path.relative_to(root).as_posix()
    assert sha(reference_path) == manifest['sha256'][reference_key]
    actual_path = study/'cycle.json'
    actual = read(actual_path)
    assert actual['reference_backend'] == 'mcc_compatible_sparse_reference.SparseCompatibleSystem'
    assert actual['reference_steps'] == reference['reference_steps'] == 100
    for name, digest in actual['input_sha256'].items():
        assert sha(Path(name)) == digest, name
    for name in ('mesh.json', 'ground.json', 'mesh-identity.json'):
        old = accepted/'inputs'/name
        assert sha(old) == manifest['sha256'][old.relative_to(root).as_posix()]
        assert sha(old) == sha(study/'inputs'/name)
    old_case = accepted/'inputs/slope-cycle-case.json'
    assert sha(old_case) == manifest['sha256'][old_case.relative_to(root).as_posix()]
    assert physical_case(read(old_case)) == physical_case(read(study/'inputs/case.json'))
    checks = []

    def field(name, a, b, tolerance):
        a, b = np.asarray(a), np.asarray(b)
        assert a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all()
        maximum = float(abs(a-b).max())
        checks.append(dict(field=name, maximum_absolute_difference=maximum,
                           tolerance=tolerance, passed=maximum < tolerance))

    assert len(actual['frames']) == len(reference['frames']) == 4
    saved_frames = sorted((study/'frames').glob('*.json'))
    assert len(saved_frames) == 4
    for i, (a, b) in enumerate(zip(actual['frames'], reference['frames'])):
        name = a['name']
        assert name == b['name'] and a['active_elements'] == b['active_elements']
        assert read(saved_frames[i])['frame'] == a, 'accepted frame snapshot differs'
        field(name+'/displacement_m', a['displacement_m'], b['displacement_m'], 1e-9)
        av, bv = np.asarray(a['material_state']), np.asarray(b['material_state'])
        field(name+'/stress_kpa', av[..., :6], bv[..., :6], 1e-6)
        field(name+'/pc_kpa', np.exp(av[..., 6]), np.exp(bv[..., 6]), 1e-6)
        field(name+'/void_ratio', np.exp(av[..., 7])-1, np.exp(bv[..., 7])-1, 1e-8)
        field(name+'/plastic_strain', av[..., 8:], bv[..., 8:], 1e-8)
        assert np.count_nonzero(np.linalg.norm(av[..., 8:], axis=-1)>1e-10) == np.count_nonzero(np.linalg.norm(bv[..., 8:], axis=-1)>1e-10)
        assert a['relative_free_force_residual'] < 1e-8
        assert a['sparse_solver_diagnostics']['maximum_true_linear_residual'] < 1e-10
    inactive = ~np.asarray(actual['frames'][-1]['active_elements'])
    assert np.array_equal(np.asarray(actual['frames'][-1]['material_state'])[inactive],
                          np.asarray(actual['frames'][-2]['material_state'])[inactive])
    passed = all(c['passed'] for c in checks)
    report = dict(schema='soil-fem-compatible-nonlinear-reference-equivalence/1', passed=passed,
                  full_field_checks=checks, accepted_frame_snapshots_exact=True,
                  removed_material_history_frozen=True, goal_complete=False, p6_complete=False,
                  spatial_convergence_accepted=False,
                  sha256={str(p): sha(p) for p in (manifest_path, reference_path, actual_path, Path(__file__))},
                  limitations='Same-mesh compatible reference backend regression, not projected/compatible continuum equivalence')
    with output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps(dict(passed=passed, checks=len(checks), failed=[c for c in checks if not c['passed']])), flush=True)
    assert passed


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ('accepted', 'study', 'output'):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    compare(Path(__file__).resolve().parents[1], args.accepted.resolve(), args.study.resolve(), args.output.resolve())
