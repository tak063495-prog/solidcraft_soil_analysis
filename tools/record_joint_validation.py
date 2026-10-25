"""Bind raw interface, global equilibrium, cutback/restart and CAD evidence."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import zipfile


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def record(run,checks):
    root=Path(__file__).resolve().parents[1];target=run/'final-validation.json'
    if target.exists():raise ValueError('acceptance record exists')
    digest=sha(run/'solver/soil-fem-3d.exe');suite=read(run/'verification.json')
    assert suite['passed'] and suite['solver_sha256']==digest
    if not suite['rust_format_lint_build_included']:assert all(p is not None for p in checks.values()),'external Rust checks must be explicit'
    defaults={'rust':'rust-tests.log','clippy':'lint.log','format':'format.log','build':'build.log'}
    checks={key:(path or run/defaults[key]).resolve() for key,path in checks.items()}
    for path in checks.values():assert path.read_text(encoding='utf-8').rstrip().endswith('Exit code: 0'),path
    names=['cad-tests','split-tests','mc','mcc','mc-history','mcc-history',
        'mc-checkpoint','mcc-checkpoint','independent','local-integration-refinement','excavation','excavation-checkpoint','local-excavation','local-excavation-checkpoint','normal-column','normal-column-theory','normal-column-checkpoint','explicit-normal-column','explicit-normal-column-theory','explicit-normal-column-checkpoint','cad','cad-independent']
    for name in names:assert (run/(name+'.log')).read_text(encoding='utf-8').rstrip().endswith('Exit code: 0'),name
    rust=checks['rust'].read_text(encoding='utf-8');count=sum(map(int,re.findall(r'test result: ok\. (\d+) passed;',rust)))
    assert count>=114
    for name in ('global_joint_closure_sliding_opening_recontact_cutback_and_restart','excavation_of_contact_side_archives_history_and_releases_only_its_old_force','two_direction_reversal_preserves_residual_slip_and_positive_dissipation'):
        assert f'test {name} ... ok' in rust
    for name in ('contact_support_counts_open_normal_stick_and_one_sided_limit_modes', 'shared_edge_and_vertex_retain_real_hinge_modes',
                 'inactive_contact_histories_are_validated_before_reporting', 'explicit_contact_supported_column_releases_upper_side_without_artificial_clamps',
                 'rigid_audit_rejects_malformed_geometry_and_nonlocal_or_nonorthogonal_transforms'):
        assert f'test {name} ... ok' in rust
    assert 'Ran 5 tests' in (run/'cad-tests.log').read_text(encoding='utf-8')
    assert 'Ran 2 tests' in (run/'split-tests.log').read_text(encoding='utf-8')
    summaries={};restarts=0
    for name in ('mc','mcc','excavation','local-excavation','normal-column','explicit-normal-column','cad/analysis'):
        folder=run/name;data=read(folder/'stage-results.json');assert data['complete'] and data['schema']=='soil-fem-stage-result/9'
        assert data['provenance']['solver_binary_sha256']==digest
        assert all(f['joint_support']['unrestrained_modes']==0 for f in data['frames'])
        package=read(folder/'accepted-stage-run.json')
        for path,expected in package['inputs_and_results_sha256'].items():assert sha(Path(path))==expected,path
        audit=read(folder/'stage-verification.json');assert audit['joint_geometry_force_energy_audit']=='passed'
        assert audit['joint_rigid_support_svd_audit']=='passed'
        assert read(folder/'results/export-verification.json')['joint_point_roundtrip']=='passed'
        assert max(f['independent_force_residual'] for f in audit['frames'])<1e-8
        assert max(f['relative_moment_error'] for f in audit['frames'])<1e-7
        if name not in ('excavation','local-excavation','normal-column','explicit-normal-column'):
            opened=next(f for f in data['frames'] if f['name']=='open');closed=next(f for f in data['frames'] if f['name']=='recontact')
            assert all(p['opening_m']>0 and p['normal_pressure_kpa']==0 for p in opened['joint_results']['points'])
            assert all(p['normal_pressure_kpa']>0 and p['opening_m']==0 for p in closed['joint_results']['points'])
            assert max(p['accumulated_plastic_slip_m'] for f in data['frames'] for p in f['state']['joints']['points'])>1e-5
        elif name=='excavation':
            assert data['excavation_audits'][0]['zero_release_relative_force_residual']<1e-8
            assert not any(data['last_committed']['excavation_release_load_kn'])
            assert data['frames'][0]['state']['joints']==data['frames'][1]['state']['joints']
            assert all(not p['active'] for p in data['frames'][1]['joint_results']['points'])
        elif name=='local-excavation':
            assert data['excavation_audits'][0]['zero_release_relative_force_residual']<1e-8
            assert abs(data['excavation_audits'][0]['removed_volume_m3']-.25)<1e-10
            before=data['frames'][0]['state']['joints']['points'];after=data['frames'][1]['state']['joints']['points']
            points=data['frames'][1]['joint_results']['points'];layout=data['joint_integration_points']
            assert any(p['active'] for p in points) and any(not p['active'] for p in points)
            assert abs(sum(layout[i]['weight_m2'] for i,p in enumerate(points) if p['active'])-.5)<1e-10
            assert all(before[i]==after[i] for i,p in enumerate(points) if not p['active'])
            assert any(before[i]!=after[i] for i,p in enumerate(points) if p['active'])
        summaries[name]={'frames':len(data['frames']),'joint_points':len(data['joint_integration_points']),
            'maximum_independent_force_residual':max(f['independent_force_residual'] for f in audit['frames'])}
    for material in ('mc','mcc'):
        h=read(run/(material+'-history/verification.json'))
        assert h['exact_final_and_all_frames'] and h['rejected_trial_joint_slip_changed'] and h['solver_sha256']==digest
        for path,expected in h['sha256'].items():assert sha(Path(path))==expected,path
    for name in ('mc-checkpoint','mcc-checkpoint','excavation-checkpoint','normal-column-checkpoint','explicit-normal-column-checkpoint','local-excavation-checkpoint'):
        cp=read(run/name/'verification.json');assert cp['solver_sha256']==digest
        assert len(cp['restarts'])==({'excavation-checkpoint':2,'local-excavation-checkpoint':2,'normal-column-checkpoint':4,'explicit-normal-column-checkpoint':4}.get(name,7))
        assert all(r['exact_states_frames_material_ids_and_audits'] for r in cp['restarts']);restarts+=len(cp['restarts'])
        for path,expected in cp['sha256'].items():assert sha(run/name/path)==expected,path
    comparisons={}
    unseating=read(run/'unseating/verification.json')
    assert unseating['passed'] and unseating['solver_sha256']==digest
    assert unseating['zero_pressure_vertical_mechanism_rejected'] and unseating['all_rejected_trials_preserved_committed_state']
    assert unseating['retained_state_agrees_with_clean_to_roundoff']
    assert unseating['clean_maximum_normalized_component_error']<=128*2.220446049250313e-16
    mode_theory=unseating['rigid_mode_theory'];assert mode_theory['passed'] and len(mode_theory['cases'])==28
    assert all(c['unrestrained_modes']==c['expected_free_modes'] for c in mode_theory['cases'])
    for path,expected in mode_theory['sha256'].items():assert sha(Path(path))==expected,path
    assert (run/'unseating.log').read_text(encoding='utf-8').rstrip().endswith('Exit code: 0')
    for path,expected in unseating['sha256'].items():assert sha(run/'unseating'/path)==expected,path
    for name in ('normal-column','explicit-normal-column'):
        theory=read(run/(name+'-theory.json'));assert theory['solver_sha256']==digest and len(theory['frames'])==4
        for path,expected in theory['sha256'].items():assert sha(Path(path))==expected,path
        assert all(f['relative_all_active_node_3d_displacement_error']<1e-7 for f in theory['frames'])
        comparisons[name]={'maximum_relative_3d_displacement_error':max(f['relative_all_active_node_3d_displacement_error'] for f in theory['frames'])}
    for name in ('independent','cad-independent','refinement'):
        comparison=read(run/(name+'.json'));assert comparison['solver_sha256']==digest
        assert len(comparison['frames'])==7 and comparison['reference_subincrements']==[100,200]
        assert comparison['maximum_independent_force_residual']<1e-8
        for path,expected in comparison['sha256'].items():assert sha(Path(path))==expected,path
        comparisons[name]={'maximum_relative_3d_displacement_error':max(f['relative_3d_displacement_error'] for f in comparison['frames']),
            'maximum_reference_refinement_error':max(f['reference_refinement_error'] for f in comparison['frames'])}
        if name=='refinement':
            refined=comparison['native_contact_refinement'];assert refined['same_global_accepted_path']
            before,after=refined['maximum_relative_3d_displacement_error_before_after'];assert after<before/2
            comparisons[name]['native_maximum_increment_m']=refined['maximum_increment_m']
            comparisons[name]['native_maximum_displacement_error_before_after']=[before,after]
    original=read(run/'cad/conformal/mesh.json');split=read(run/'cad/split/mesh.json');provenance=read(run/'cad/split/mesh.provenance.json')['joint_split']
    assert split['material_ids']==original['material_ids'] and split['region_ids']==original['region_ids']
    assert provenance['source_mesh_sha256']==sha(run/'cad/conformal/mesh.json') and provenance['added_displacement_nodes']>0
    compatibility=run/'mcc-regression.json'
    if compatibility.exists():
        c=read(compatibility);assert c['exact_committed_state'] and c['exact_frames'] and c['current_solver_sha256']==digest
        for path,expected in c['files_sha256'].items():assert sha(Path(path))==expected,path
    joint_compatibility=run/'contact-state-compatibility.json'
    if joint_compatibility.exists():
        c=read(joint_compatibility);assert c['passed'] and len(c['comparisons'])==2
        for case in c['comparisons']:
            assert case['exact_committed_state'] and case['exact_frames_except_added_support_report']
            assert case['current_solver_sha256']==digest
        for path,expected in c['sha256'].items():assert sha(Path(path))==expected,path
    sources=[root/name for name in ('Cargo.lock','Cargo.toml','rust-toolchain.toml','requirements.lock.txt','sources.lock.json','README.md')]
    sources += [p for name in ('src','tests','examples','scripts','docs') for p in (root/name).rglob('*') if p.is_file() and p.suffix in ('.rs','.py','.json','.ps1','.md')]
    sources += list((root/'tools').glob('*.py'));manifest={str(p.relative_to(root)):sha(p) for p in sorted(set(sources))}
    archive=run/'source-snapshot.zip'
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for name,expected in manifest.items():
            payload=(root/name).read_bytes();assert hashlib.sha256(payload).hexdigest()==expected;z.writestr(name,payload)
    (run/'source-snapshot.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    report={'schema':'soil-fem-joint-final-validation/1','passed':True,'solver_sha256':digest,'rust_tests':count,'cad_and_split_tests':7,
        'separate_process_exact_restarts':restarts,'mc_mcc_cutback_exact':True,'joint_side_excavation_force_release_and_frozen_history':True,
        'explicit_normal_contact_support_and_one_sided_rigid_mode_audit':True,
        'zero_pressure_rigid_mechanism_rejected_with_history_preserved':True,
        'examples':summaries,'comparisons':comparisons,'legacy_mcc_exact':compatibility.exists(),
        'limitations':'Dry, small-strain, fixed CAD normals, cohesion=0/psi=0, fully matching Tri3/Tri6 interfaces. No tied crack fronts or interface BIRTH reactivation. Global independent reference uses elastic MC bulk with nonlinear frictional JOINT; plastic MCC bulk has history/force checks, not an independent plastic 3D displacement reference. P6 overall remains incomplete.',
        'rust_checks_external':not suite['rust_format_lint_build_included'],
        'source_snapshot_sha256':sha(archive),'sha256':{str(p.relative_to(root)):sha(p) for p in sorted(set(p for p in run.rglob('*') if p.is_file())|set(checks.values()))}}
    with target.open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True)
    for name in ('rust','clippy','format','build'):p.add_argument('--'+name+'-log',type=Path)
    a=p.parse_args();record(a.run.resolve(),{name:getattr(a,name+'_log') for name in ('rust','clippy','format','build')})
