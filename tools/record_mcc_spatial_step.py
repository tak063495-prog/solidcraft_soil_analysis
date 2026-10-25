"""Seal completed spatial-step evidence while leaving convergence unaccepted."""
import argparse
import ast
import json
from pathlib import Path
import zipfile
import freeze_sources
from mcc_projected_volume_diagnostic import read, sha
import mcc_projected_cycle_diagnostic as cycle
from mcc_projected_sparse_reference import SparseProjectedSystem


def archive_check(folder):
    expected = read(folder/'source-snapshot.json')
    with zipfile.ZipFile(folder/'source-snapshot.zip') as archive:
        assert archive.testzip() is None and set(archive.namelist()) == set(expected)
        import hashlib
        for name, digest in expected.items():
            assert hashlib.sha256(archive.read(name)).hexdigest() == digest, name
    return len(expected)


def symbol(text, name):
    return next(ast.dump(n, include_attributes=False) for n in ast.parse(text).body if getattr(n, 'name', None) == name)


def record(root, output):
    output.mkdir(parents=True, exist_ok=False)
    freeze_sources.freeze(output)
    fine = root/'runs/projected-sparse-spatial-3m-20261025-v1'
    job = root/'runs/projected-mcc-spatial-inputs-20261025-v1/mesh-3m/ultrafine-job'
    compatible = root/'runs/compatible-sparse-cycle-20261025-v1'
    rates = root/'runs/compatible-sparse-rate-equivalence-20261025-v1'
    assessment_path = root/'runs/projected-mcc-spatial-3m-assessment-20261025-v1/assessment.json'
    assessment = read(assessment_path)
    assert assessment['field_accuracy_passed'] and assessment['reference_refinement_passed'] is None
    assert not assessment['integration_gates_passed'] and not assessment['spatial_convergence_verified']
    for name, digest in assessment['sha256'].items():
        assert sha(Path(name)) == digest, name
    assert read(job/'status.json')['state'] == 'succeeded'
    assert read(fine/'progress.json')['status'] == 'completed'
    reference = read(fine/'cycle.json')
    for name, digest in reference['input_and_source_sha256'].items():
        assert sha(Path(name)) == digest, name
    regression_path = compatible/'backend-equivalence.json'
    regression = read(regression_path)
    assert regression['passed'] and len(regression['full_field_checks']) == 20
    for name, digest in regression['sha256'].items():
        assert sha(Path(name)) == digest, name
    assert read(rates/'rate-equivalence.json')['passed']
    archives = {str(p): archive_check(p) for p in (fine, compatible, rates)}
    current_rate = (root/'tools/mcc_rate_fem.py').read_text(encoding='utf-8')
    for folder in (fine, rates):
        with zipfile.ZipFile(folder/'source-snapshot.zip') as archive:
            old = archive.read('tools/mcc_rate_fem.py').decode('utf-8')
        for name in ('System', 'terms'):
            assert symbol(old, name) == symbol(current_rate, name), (str(folder), name)
    assert sha(root/'tools/mcc_compatible_sparse_reference.py') == sha(compatible/'reference-sources/mcc_compatible_sparse_reference.py')
    assert sha(root/'tools/mcc_compatible_sparse_run.py') == sha(compatible/'reference-sources/mcc_compatible_sparse_run.py')
    assert sha(root/'tools/mcc_rate_fem.py') == sha(compatible/'reference-sources/mcc_rate_fem.py')
    bad_case = read(fine/'inputs/case.json')
    bad_case['settings']['kinematics'] = 'explicit'
    wrong = output/'changed-kinematics.json'
    wrong.write_text(json.dumps(bad_case), encoding='utf-8')
    refused = False
    try:
        cycle.solve(fine/'inputs/mesh.json', fine/'inputs/ground.json', wrong,
                    fine/'inputs/mesh-identity.json', fine/'initial.json', output/'must-not-exist.json',
                    2, system_factory=SparseProjectedSystem)
    except AssertionError as error:
        refused = str(error) == 'initial physical case differs'
    assert refused and not (output/'must-not-exist.json').exists()
    files = [assessment_path, regression_path, rates/'rate-equivalence.json', fine/'initial.json', fine/'cycle.json',
             fine/'progress.json', job/'job.json', job/'status.json', job/'result.json',
             root/'runs/projected-mcc-spatial-inputs-20261025-v1/area-measures-step2.json',
             root/'runs/projected-mcc-spatial-probes-20261025-v1/points.json',
             root/'runs/projected-mcc-spatial-probes-20261025-v1/step2.json',
             root/'runs/projected-mcc-spatial-probes-20261025-v1/rust-independent-413.json',
             root/'runs/projected-mcc-spatial-probes-20261025-v1/independent-224-to-413.json', wrong]
    assert read(output/'source-snapshot.json') == freeze_sources.manifest()
    report = dict(schema='soil-fem-mcc-spatial-step-record/1', evidence_integrity_passed=True,
                  scope='Completed 413-cell fixed-mesh field comparison, coarse compatible backend regression and first spatial sensitivity',
                  compatible_rate_cases=12, compatible_nonlinear_field_checks=20,
                  fine_fixedmesh_field_accuracy_passed=True, reference_refinement_passed=None,
                  integration_gates_passed=False, spatial_convergence_verified=False,
                  native_ui_connected=False, production_backend_accepted=False, goal_complete=False, p6_complete=False,
                  reference_material_and_system_ast_unchanged=True, changed_kinematics_refused=True,
                  source_archives_verified=archives, evidence_sha256={str(p): sha(p) for p in files},
                  limitations='Independent 100-to-200 integration refinement, 723-cell results, compatible fine mesh, further spatial levels and all remaining P0-P6 work are pending')
    with (output/'verification.json').open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps(dict(evidence_integrity_passed=True, fine_fixedmesh_field_accuracy_passed=True,
                          reference_refinement_passed=None, spatial_convergence_verified=False)), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    record(Path(__file__).resolve().parents[1], args.output.resolve())
