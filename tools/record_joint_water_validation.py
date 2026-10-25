"""Seal the complete prior suite plus fresh wet-interface integration proof."""
import argparse
import json
from pathlib import Path
from freeze_sources import manifest
from verify_acceptance_packet import verify,sha

def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def record(run):
    root=Path(__file__).resolve().parents[1];assert read(run/'source-snapshot.json')==manifest(),'sources changed during acceptance'
    parent=verify(run/'baseline');old=read(run/'baseline/final-validation.json');cases=read(run/'cases/verification.json');digest=sha(run/'solver/soil-fem-3d.exe');assert cases['passed'] and cases['cad_exercised'] and digest==parent['solver_sha256']==cases['solver_sha256']==sha(root/'target/release/soil-fem-3d.exe')
    assert old['rust_tests']>=206 and old['lifecycle_separate_process_restarts']==33 and old['construction_separate_process_restarts']==28 and old['fixed_domain_coupled_restarts']==29 and old['standalone_hydraulic_restarts']==8 and old['exact_legacy_cases']==6
    assert len(cases['profiles'])==15 and len(cases['restarts'])==34 and all(c['all_states_frames_attempts_flow_maps_and_water_histories_exact'] for c in cases['restarts'])
    assert len(cases['dry_joint_legacy_cases'])==2 and all(c['exact_frames_and_committed_state'] for c in cases['dry_joint_legacy_cases'])
    for n in ('baseline','reference-tests','cases'):assert (run/f'{n}.log').read_text(encoding='utf-8').rstrip().endswith('Exit code: 0')
    assert 'Ran 3 tests' in (run/'reference-tests.log').read_text(encoding='utf-8')
    rust=(run/'baseline/baseline/baseline/baseline/core/rust-tests.log').read_text(encoding='utf-8')
    for n in ('hydrostatic_contact_uses_effective_spring_pressure_and_total_equals_effective_plus_water','pressure_survives_open_contact_and_does_not_enter_friction_twice','partial_submergence_integrates_actual_sloped_tri3_and_tri6_force_and_moment','continuous_heads_pass_series_flow_through_split_displacement_nodes','soil_structure_is_impermeable_and_prescribed_heads_are_not_averaged','wet_stage_cutback_restores_pressure_contact_history_and_checkpoint','staged_darcy_interface_policies_preserve_primary_pressure_and_exact_midstage_restart','death_after_wet_seepage_releases_interface_water_and_archives_contact_history','wet_global_sliding_opening_recontact_preserves_water_load_and_rejected_slip_history'):assert n+' ... ok' in rust,n
    for n in ('flow_cutback_clean_path_exact','wet_contact_open_pressure_and_rejected_slip_history_verified','negative_saturated_flow_target_preserves_committed_state','missing_pressure_policy_rejected_before_output'):assert cases[n]
    for name,checked in cases['profiles'].items():
        folder=run/'cases'/name
        if name=='cad-local':folder=folder/'analysis'
        data=read(folder/'stage-results.json');assert data['complete'] and data['provenance']['solver_binary_sha256']==digest and checked['passed'];exports=read(folder/'results/export-verification.json');assert exports['vtu_roundtrip']=='passed' and exports['joint_point_roundtrip']=='passed'
    cad=read(run/'cases/cad-local/contact-domain-verification.json');assert cad['passed'] and cad['removed_wet_contact_points']>0 and cad['retained_wet_contact_points']>0
    files=sorted(p for p in run.rglob('*') if p.is_file())
    result=dict(schema='soil-fem-wet-joint-final-validation/1',passed=True,solver_sha256=digest,source_snapshot_sha256=sha(run/'source-snapshot.zip'),parent_manifest_sha256=parent['manifest_sha256'],rust_tests=old['rust_tests'],cad_tests=5,joint_split_tests=2,independent_reference_tests=3,wet_joint_profiles=15,wet_joint_separate_process_restarts=34,lifecycle_separate_process_restarts=33,construction_separate_process_restarts=28,fixed_domain_coupled_restarts=29,standalone_hydraulic_restarts=8,exact_legacy_cases=6,exact_dry_joint_cases=2,cad_contact_domain=cad,maximum_independent_force_residual=max(f['relative_force_residual'] for c in cases['profiles'].values() for f in c['frames']),maximum_independent_physical_load_component_error_kn=max(c['maximum_physical_load_component_error_kn'] for c in cases['profiles'].values()),maximum_independent_hydraulic_head_difference_m=max(c['maximum_independent_hydraulic_head_difference_m'] for c in cases['profiles'].values()),scope='Hydrostatic and steady saturated Darcy JOINT: SOIL/CONTINUOUS/PRESCRIBED, independently split displacements and tied heads, explicit impermeable total-stress MC structure role, separate physical-pressure loads and effective contact friction, partial wet integration, partial DEATH, pressure/contact history rollback and exact restart. Independent Darcy FE and force/contact reassembly with analytic column/resultant checks. No unsaturated joint transport, BIRTH contact reactivation, general independent plastic 3D mechanics, plates/piles or full P6 completion.',sha256={p.relative_to(root).as_posix():sha(p) for p in files})
    with (run/'final-validation.json').open('x',encoding='utf-8') as stream:json.dump(result,stream,indent=2)
    print(json.dumps({k:v for k,v in result.items() if k!='sha256'},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',required=True,type=Path);record(p.parse_args().run.resolve())
