"""Seal the full current engine plus physical-time ONE_WAY evidence."""
import argparse
import json
from pathlib import Path
from freeze_sources import manifest
from verify_acceptance_packet import verify as packet,sha


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def record(run):
    root=Path(__file__).resolve().parents[1];assert not (run/'final-validation.json').exists()
    assert read(run/'source-snapshot.json')==manifest(),'source changed during acceptance'
    baseline=packet(run/'baseline');digest=sha(run/'solver/soil-fem-3d.exe');assert digest==baseline['solver_sha256']==sha(root/'target/release/soil-fem-3d.exe')
    for name in ('baseline','cad-slope','cases','convergence'):assert (run/f'{name}.log').read_text(encoding='utf-8').rstrip().endswith('Exit code: 0'),name
    old=read(run/'baseline/final-validation.json');assert old['rust_tests']>=174 and old['separate_process_restarts']==8 and old['exact_legacy_cases']==6
    rust=(run/'baseline/core/rust-tests.log').read_text(encoding='utf-8')
    for name in ('genuine_mechanical_exhaustion_does_not_commit_the_successful_hydraulic_trial',
        'plastic_coupled_cutbacks_restore_both_histories_and_same_accepted_times_replay_exactly',
        'checkpoints_replay_initialization_and_physical_periods_and_reject_resealed_changes',
        'declared_compressible_water_mass_includes_positive_storage_and_closed_head_is_preserved',
        'bishop_zero_one_and_separate_exterior_pressure_match_exact_elastic_limits'):
        assert name+' ... ok' in rust,name
    cases=read(run/'cases/verification.json');convergence=read(run/'convergence/verification.json')
    assert cases['passed'] and cases['solver_sha256']==digest and cases['cad_slope_exercised']
    assert cases['derived_point_omission_preserves_exact_primary_state']
    assert len(cases['restarts'])==29 and all(r['all_states_frames_and_attempts_exact'] for r in cases['restarts'])
    assert cases['genuine_mechanical_exhaustion_preserves_both_states'] and cases['plastic_rejection_and_clean_accepted_times_exact']
    for name in ('mc-tet4','mc-tet10','mcc-tet4','mcc-tet10','cad-slope','omitted-points'):
        data=read(run/f'cases/{name}/result.json');assert data['complete'] and data['solver_binary_sha256']==digest
        assert read(run/f'cases/{name}/verification.json')['passed']
        assert read(run/f'cases/{name}/results/export-verification.json')['all_vtu_and_pvd_roundtrip']=='passed'
    assert convergence['passed'] and convergence['solver_sha256']==digest
    assert len(convergence['comparisons'])==9 and len(convergence['native_spatial_refinement'])==3 and len(convergence['native_temporal_refinement'])==3
    for label in ('rt0-16-dt50','rt0-32-dt50','rt0-64-dt50','rt0-32-dt25','rt0-32-dt12.5'):
        assert read(run/f'convergence/{label}/result.json')['solver_binary_sha256']==digest
        assert read(run/f'convergence/{label}/verification.json')['passed']
    cad=read(run/'cad/mesh.json');assert cad['tetrahedra10'] and len(read(run/'cad/bodies.json')['regions'])==1
    assert read(run/'cases/cad-slope/cad-input-link.json')['mesh_contents_equal']
    assert sha(run/'cad/solver/soil-fem-3d.exe')==digest
    files=[p for p in run.rglob('*') if p.is_file()]
    report=dict(schema='soil-fem-one-way-final-validation/1',passed=True,solver_sha256=digest,source_snapshot_sha256=sha(run/'source-snapshot.zip'),
        rust_tests=old['rust_tests'],cad_tests=5,joint_split_tests=2,coupled_separate_process_restarts=29,standalone_hydraulic_restarts=8,exact_legacy_cases=6,
        cad_slope_nodes=len(cad['nodes']),cad_slope_tet10=len(cad['tetrahedra10']),baseline_manifest_sha256=baseline['manifest_sha256'],
        genuine_mechanical_exhaustion_preserves_pair=True,plastic_pair_cutback_and_clean_time_replay_exact=True,derived_point_omission_exact=True,
        independent_elastic_comparisons=cases['elastic_comparisons'],independent_coupled_column_comparisons=convergence['comparisons'],native_spatial_refinement=convergence['native_spatial_refinement'],native_temporal_refinement=convergence['native_temporal_refinement'],independent_reference_refinement=convergence['independent_reference_refinement'],
        maximum_independent_mechanical_residual=cases['maximum_independent_mechanical_residual'],maximum_column_hydraulic_integrated_mass_error_m3=convergence['maximum_hydraulic_integrated_mass_error_m3'],
        scope='Fixed active-domain, frozen hydraulic porosity Richards to MC/MCC mechanics, explicit dry reference weight and water-storage mass, separate physical u/chi*u/exterior loads, initial equilibrium and physical-time paired rollback/restart. Independent elastic 3D mechanics given water heads plus independent 1D FV/continuous elastic column. No arbitrary 3D plastic displacement, excavation/BIRTH water transfer, wet JOINT, pond dynamics or slope stability certification.',
        sha256={p.relative_to(root).as_posix():sha(p) for p in sorted(files)})
    with (run/'final-validation.json').open('x',encoding='utf-8') as f:json.dump(report,f,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k not in ('sha256','independent_elastic_comparisons','independent_coupled_column_comparisons','native_spatial_refinement','native_temporal_refinement','independent_reference_refinement')},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);record(p.parse_args().run.resolve())
