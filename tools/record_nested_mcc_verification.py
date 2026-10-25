"""Seal mesh-refinement/input evidence and remaining discretization differences."""
import argparse
import json
from pathlib import Path
import freeze_sources
from mcc_projected_volume_diagnostic import read, sha
from record_mcc_spatial_step import archive_check


def verify_bindings(report, key):
    for name, digest in report[key].items():
        assert sha(Path(name)) == digest, name


def record(root, output):
    output.mkdir(parents=True, exist_ok=False)
    freeze_sources.freeze(output)
    geometry = root/'runs/tet10-nested-refinement-verification-20261025-v2'
    prepared = root/'runs/projected-mcc-nested-inputs-20261025-v1'
    differences = root/'runs/mcc-discretization-spatial-20261025-v2'
    measured = root/'runs/mcc-discretization-spatial-20261025-v1'
    geometry_report = read(geometry/'verification.json')
    assert geometry_report['passed'] and len(geometry_report['cases']) == 4
    verify_bindings(geometry_report, 'sha256')
    assert sha(root/'tools/refine_tet10.py') == read(geometry/'source-snapshot.json')['tools/refine_tet10.py']
    assert read(geometry/'actual_solidcraft/refined.json') == read(prepared/'inputs/mesh.json')
    preparation = read(prepared/'preparation.json')
    assert preparation['input_preparation_passed'] and preparation['elements'] == 1792
    assert preparation['integration_points'] == 48384
    verify_bindings(preparation, 'original_input_sha256')
    verify_bindings(preparation, 'prepared_sha256')
    identity, ground = read(prepared/'inputs/mesh-identity.json'), read(prepared/'inputs/ground.json')
    assert ground['initialization']['geostatic']['seed']['mesh_sha256'] == identity['mesh_sha256']
    points = ground['initialization']['geostatic']['seed']['points']
    assert len(points) == len(identity['integration_points']) == 48384
    assert all(p['element_index'] == q['element_index'] and p['local_index'] == q['local_index'] and
               p['position_m'] == q['position_m'] for p, q in zip(points, identity['integration_points']))
    job = prepared/'ultrafine-job'
    manifest, status = read(job/'job.json'), read(job/'status.json')
    assert status['state'] == 'prepared' and status['process_id'] is None
    assert manifest['binding']['identity']['mesh_sha256'] == identity['mesh_sha256']
    for name, digest in manifest['binding']['files_sha256'].items():
        assert sha(job/name) == digest, name
    with (output/'prepared-job-status.json').open('x', encoding='utf-8') as stream:
        json.dump(status, stream, indent=2)
    models = [read(differences/name) for name in ('coarse.json', 'fine.json')]
    for model in models:
        verify_bindings(model, 'sha256')
        assert not model['spatial_convergence_verified']
        for frame in model['frames']:
            integrated = frame['physical_volume_displacement_difference']
            assert max(integrated['rms_change_under_quadrature_refinement_m'].values()) < 1e-12
            assert frame['projected_force_residual'] < 1e-8 and frame['compatible_force_residual'] < 1e-8
    assert read(root/'runs/compatible-sparse-spatial-3m-20261025-v1/progress.json')['status'] == 'completed'
    archives = {str(p): archive_check(p) for p in (geometry, prepared, differences)}
    build = read(root/'runs/projected-mcc-development-20261025-v1/build-sources/source-snapshot.json')
    compiled = {name: digest for name, digest in build.items() if name.startswith(('src/', 'tests/', 'examples/')) or name in ('Cargo.toml', 'Cargo.lock', 'rust-toolchain.toml')}
    assert all(sha(root/name) == digest for name, digest in compiled.items())
    evidence = [geometry/'verification.json', prepared/'preparation.json', job/'job.json',
                output/'prepared-job-status.json', differences/'coarse.json', differences/'fine.json',
                measured/'area-measures.json', measured/'probe-method-difference.json',
                measured/'probe-method-difference-refined.json',
                root/'runs/compatible-sparse-spatial-3m-20261025-v1/cycle.json',
                root/'runs/tet10-nested-refinement-verification-20261025-v1.log']
    assert read(output/'source-snapshot.json') == freeze_sources.manifest()
    report = dict(schema='soil-fem-nested-mcc-verification-record/1', evidence_integrity_passed=True,
                  scope='Nested geometry/trace tests, fresh initial-input preparation, completed compatible fine reference and unaccepted formulation differences',
                  nested_geometry_cases_passed=4, fresh_initial_points=48384,
                  nested_job_state_at_capture='prepared', nested_nonlinear_solution_verified=False,
                  source_archives_verified=archives, unchanged_rust_build_and_examples=len(compiled),
                  displacement_difference_volume_integral_quadrature_checked=[5, 7],
                  post_cut_relative_formulation_differences=[model['frames'][-1]['physical_volume_displacement_difference']['order7']['fields']['from_seed']['relative_l2'] for model in models],
                  spatial_convergence_verified=False, production_backend_accepted=False, goal_complete=False, p6_complete=False,
                  evidence_sha256={str(p): sha(p) for p in evidence},
                  limitations='1792-cell calculation not started; 723-cell comparison, 413-cell reference integration refinement and all other P0-P6 requirements remain pending')
    with (output/'verification.json').open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps(dict(evidence_integrity_passed=True, nested_geometry_cases_passed=4,
                          fresh_initial_points=48384, nested_job_state_at_capture='prepared',
                          spatial_convergence_verified=False)), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    record(Path(__file__).resolve().parents[1], args.output.resolve())
