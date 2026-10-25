"""Seal the actual 2D/3D confined-column comparison after all final gates."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
from freeze_sources import freeze


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record(base, assessment, log, output):
    root = Path(__file__).resolve().parents[1]
    report = read(assessment)
    assert report['passed'] and report['bounds_unchanged_from_first_assessment']
    assert report['accepted_native_resolution'] == '3d-n4-fine'
    assert report['native'][-1]['passed']
    assert not report['native'][0]['passed'] and not report['native'][1]['passed']
    assert all(r['all_spatial_errors_reduced'] for r in report['refinement'])
    assert sha(root/'tools/verify_2d_column_reference.py') == report['verifier_sha256']
    assert log.read_text(encoding='utf-8').rstrip().endswith('Exit code: 0')
    for name, expected in report['input_and_evidence_sha256'].items():
        assert sha(base/name) == expected, name
    revision = read(base/'provenance.json')['revision']
    external = root/'external/excel-vba-soil-fem'
    assert subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=external, text=True).strip() == revision
    assert revision == '440d246d9e32f5892c79539d4a7b797ad29442ee'
    solver = root/'target/release/soil-fem-3d.exe'
    parent = root/'runs/default-solver-equivalence-20261025-v1'
    accepted = read(parent/'final-validation.json')
    assert accepted['passed'] and sha(solver) == accepted['solver_sha256']
    assert read(parent/'integrity-verification.json')['passed']
    # The current verification additions do not replace the build sources of
    # the accepted executable. Require every Rust/build file to remain exact.
    for name, expected in read(parent/'source-snapshot.json').items():
        if name.startswith('src/') or name in ('Cargo.toml', 'Cargo.lock', 'rust-toolchain.toml'):
            assert sha(root/name) == expected, name
    output.mkdir(exist_ok=False)
    (output/'solver').mkdir()
    shutil.copyfile(solver, output/'solver/soil-fem-3d.exe')
    shutil.copytree(base, output/'data')
    shutil.copyfile(log, output/'assessment.log')
    (output/'additional-logs').mkdir()
    for name in ('2d-3d-column-yield-front-20261025-v1.log',
                 '2d-3d-column-yield-front-assessment-20261025-v1.log',
                 '2d-3d-column-yield-front-ny2-excel-20261025-v1.log',
                 '2d-3d-column-yield-front-ny4-excel-20261025-v1.log',
                 '2d-3d-column-yield-front-excel-assessment-20261025-v1.log'):
        source = root/'runs'/name
        assert source.read_text(encoding='utf-8').rstrip().endswith('Exit code: 0'), name
        shutil.copyfile(source, output/'additional-logs'/name)
    (output/'published-2d-source').mkdir()
    for name in ('LICENSE', 'README.md', 'docs/MODEL_INPUTS.md', 'docs/INCONSISTENT_REPAIR.md',
                 'src/vba/FEM.bas', 'src/vba/FEMIo.bas', 'src/vba/FEMPlasticLog.bas',
                 'src/vba/FEMEngine.bas', 'src/vba/FEMSolver.bas'):
        target = output/'published-2d-source'/name
        target.parent.mkdir(exist_ok=True, parents=True)
        shutil.copyfile(external/name, target)
    shutil.copyfile(parent/'source-snapshot.json', output/'build-source-manifest.json')
    freeze(output)
    final = dict(schema='soil-fem-actual-2d-3d-column-final-validation/1', passed=True,
        goal_complete=False, p6_complete=False,
        scope=report['scope'], published_2d_revision=revision,
        accepted_3d_resolution='384 Tet10 / 729 nodes / 10368 GP',
        solver_sha256=sha(solver), source_snapshot_sha256=sha(output/'source-snapshot.zip'),
        parent_default_build_manifest_sha256=sha(parent/'final-validation.json'),
        coarse_and_intermediate_resolution_failures_retained=True,
        numerical_bounds_not_relaxed=True,
        excel_program_code_unchanged=True,
        actual_2d_workbook_and_csv_all_node_gp_fields_compared=True,
        all_3d_nodes_gp_stress_and_cumulative_plastic_strain_compared=True,
        cumulative_plastic_history_exact_on_unloading=True,
        partially_yielding_front_diagnostics_included=True,
        partially_yielding_3d_path_is_not_accepted_by_final_endpoint_gates=True,
        refinement=report['refinement'],
        limitations=report['limitations'],
        sha256={p.relative_to(root).as_posix(): sha(p) for p in output.rglob('*') if p.is_file()})
    with (output/'final-validation.json').open('x', encoding='utf-8') as stream:
        json.dump(final, stream, indent=2)
    print(json.dumps({k: final[k] for k in ('passed', 'scope', 'published_2d_revision', 'accepted_3d_resolution')}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('base', type=Path)
    parser.add_argument('assessment', type=Path)
    parser.add_argument('log', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    record(*(p.resolve() for p in (args.base, args.assessment, args.log, args.output)))
