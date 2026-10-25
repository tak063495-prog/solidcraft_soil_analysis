"""Seal hydraulic-core proof against sources frozen before compilation."""
import argparse
import hashlib
import json
from pathlib import Path
import re
from freeze_sources import manifest


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def record(run):
    root=Path(__file__).resolve().parents[1];target=run/'final-validation.json'
    assert not target.exists() and read(run/'source-snapshot.json')==manifest(),'sources changed after freeze'
    sources=read(run/'source-snapshot.json');assert 'tests/data/retention-decimal-reference.csv' in sources
    digest=sha(run/'solver/soil-fem-3d.exe');suite=read(run/'verification.json');assert suite['passed'] and suite['solver_sha256']==digest
    logs=['freeze','format','rust-tests','lint','build','cad-tests','split-tests','decimal','cad-capture','cases','convergence','legacy']
    for name in logs:assert (run/f'{name}.log').read_text(encoding='utf-8').rstrip().endswith('Exit code: 0'),name
    rust=(run/'rust-tests.log').read_text(encoding='utf-8');count=sum(map(int,re.findall(r'test result: ok\. (\d+) passed;',rust)));assert count>=164
    for name in ('tight_flow_regression_uses_original_residual_increments_and_conductivity_tangent',
        'rain_pause_restart_checkpoints_restore_all_heads_storage_faces_and_step_history',
        'nonlinear_exhaustion_keeps_the_initial_field_and_steady_nonuniqueness_is_not_success',
        'rt0_mass_matches_exact_simplex_moments_and_piola_preserves_curved_face_fluxes',
        'zero_storage_saturated_closed_domain_does_not_certify_an_undetermined_pressure_reference'):
        assert name+' ... ok' in rust,name
    assert 'Ran 5 tests' in (run/'cad-tests.log').read_text(encoding='utf-8-sig') and 'Ran 2 tests' in (run/'split-tests.log').read_text(encoding='utf-8-sig')
    cases=read(run/'cases/verification.json');assert cases['passed'] and cases['solver_sha256']==digest and cases['cad_curved_hydrostatic_exercised']
    assert cases['zero_storage_missing_pressure_reference_rejected']
    assert len(cases['restarts'])==8 and all(r['exact_all_states_frames_and_attempts'] for r in cases['restarts'])
    for name in ('wet','dry','utf8','rain-tet4','rain-tet10','cutback','clean','tight-flow','cad-curved-hydrostatic'):
        data=read(run/f'cases/{name}/result.json');assert data['complete'] and data['solver_binary_sha256']==digest
        assert read(run/f'cases/{name}/verification.json')['passed']
        if name not in ('cutback','clean','tight-flow'):
            export=read(run/f'cases/{name}/results/export-verification.json');assert export['passed'] and export['cell_and_raw_point_roundtrip']=='passed'
    convergence=read(run/'convergence/verification.json');assert convergence['passed'] and convergence['solver_sha256']==digest
    assert len(convergence['comparisons'])==9 and len(convergence['spatial_refinement'])==3 and len(convergence['temporal_refinement'])==3
    for label in ('rt0-32-dt50','rt0-64-dt50','rt0-128-dt50','rt0-64-dt25','rt0-64-dt12.5'):
        result=read(run/f'convergence/{label}/result.json');assert result['complete'] and result['solver_binary_sha256']==digest
        assert read(run/f'convergence/{label}/verification.json')['passed']
    assert read(run/'legacy/seepage-legacy-verification.json')['passed']
    assert sha(run/'cad/solver/soil-fem-3d.exe')==digest
    cad=read(run/'cad/mesh.json');assert cad['tetrahedra10'] and read(run/'cad/bodies.json')['excavation_tools']
    files=[p for p in run.rglob('*') if p.is_file()]
    report=dict(schema='soil-fem-richards-core-final-validation/1',passed=True,solver_sha256=digest,rust_tests=count,cad_tests=5,joint_split_tests=2,
        separate_process_restarts=8,exact_legacy_cases=6,cad_curved_nodes=len(cad['nodes']),cad_curved_tet10=len(cad['tetrahedra10']),
        genuine_cutback_and_exhaustion_verified=cases['exact_cutback_accepted_history'] and cases['genuine_exhaustion_preserves_initial_state'],
        independent_column_comparisons=convergence['comparisons'],spatial_refinement=convergence['spatial_refinement'],temporal_refinement=convergence['temporal_refinement'],
        maximum_native_integrated_mass_error_m3=convergence['maximum_native_integrated_mass_error_m3'],
        maximum_reference_integrated_mass_error_m=convergence['maximum_reference_integrated_mass_error_m'],
        scope='Hydraulic core: 3D RT0 geometry, VG/Mualem, conservative BE, no-ponding atmosphere, exact restart, independent symmetry-reduced 1D FV and space/time refinement. Not arbitrary finite-width 3D free-surface/seepage-length validation, pond storage/routing, or suction-to-soil mechanics.',
        source_snapshot_sha256=sha(run/'source-snapshot.zip'),sha256={p.relative_to(root).as_posix():sha(p) for p in sorted(files)})
    with target.open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k not in ('sha256','independent_column_comparisons','spatial_refinement','temporal_refinement')},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);record(p.parse_args().run.resolve())
