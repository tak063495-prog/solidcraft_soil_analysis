"""Bind local-cut CAD, mechanics, independent rates and history evidence.

No global plastic 3D reference or spatial convergence is inferred from these
checks. Each acceptance run is immutable; source files are frozen in a ZIP.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import zipfile
from freeze_sources import manifest


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def record(run, checks):
    root=Path(__file__).resolve().parents[1]
    target=run/'final-validation.json'
    if target.exists():
        raise ValueError('acceptance record exists')
    solver=run/'solver/soil-fem-3d.exe';digest=sha(solver)
    suite=read(run/'verification.json')
    assert suite['passed'] and suite['solver_sha256']==digest
    for path in checks.values():
        assert path.read_text(encoding='utf-8').rstrip().endswith('Exit code: 0'),path
    rust=checks['rust'].read_text(encoding='utf-8')
    count=sum(map(int,re.findall(r'test result: ok\. (\d+) passed;',rust)))
    assert count>=87
    for name in ('local_cutback_restores_plastic_history_loads_and_active_elements',
                 'sparse_newton_matches_dense_for_plastic_local_excavation_and_reports_true_residual',
                 'submerged_cut_has_hydrostatic_wall_and_floor_tractions',
                 'quadratic_mcc_excavation_cutback_and_checkpoint_restore_all_quadrature_history'):
        assert f'test {name} ... ok' in rust
    cad=(run/'cad-tests.log').read_text(encoding='utf-8')
    assert 'Ran 5 tests' in cad and cad.rstrip().endswith('Exit code: 0')
    examples={}
    for name in ('flat','curved'):
        folder=run/name;accepted=read(folder/'accepted-stage-run.json')
        assert accepted['solver_binary_sha256']==digest==sha(folder/'solver/soil-fem-3d.exe')
        for path,expected in accepted['inputs_and_results_sha256'].items():
            assert sha(Path(path))==expected,path
        mesh=read(folder/'ground/mesh.json');data=read(folder/'stage-results.json')
        assert data['complete'] and len(data['frames'])==2 and mesh['tetrahedra10']
        audit,=data['excavation_audits']
        assert audit['zero_release_relative_force_residual']<1e-8
        assert not any(data['last_committed']['excavation_release_load_kn'])
        if name=='flat':
            assert abs(audit['removed_volume_m3']-32)<1e-9
        independent=read(folder/'stage-verification.json')
        assert independent['vtu_roundtrip']=='passed' and independent['raw_integration_point_roundtrip']=='passed'
        end=independent['frames'][-1]
        assert end['independent_force_residual']<1e-8 and end['relative_moment_error']<1e-7
        monitors=data['frames'][-1]['monitors']
        assert monitors['excavation_floor']['mean_displacement_z_m']>0
        assert monitors['surrounding_surface']['minimum_displacement_z_m']<0
        examples[name]={
            'nodes':len(mesh['nodes']),'tetrahedra10':len(mesh['tetrahedra10']),
            'integration_points':len(data['integration_points']),
            'removed_elements':audit['removed_element_count'],'removed_volume_m3':audit['removed_volume_m3'],
            'zero_release_relative_force_residual':audit['zero_release_relative_force_residual'],
            'final_independent_force_residual':end['independent_force_residual'],
            'relative_moment_error':end['relative_moment_error'],
            'floor_nodal_mean_heave_mm':1000*monitors['excavation_floor']['mean_displacement_z_m'],
            'surrounding_surface_minimum_uz_mm':1000*monitors['surrounding_surface']['minimum_displacement_z_m'],
            'surrounding_surface_nodal_mean_uz_mm':1000*monitors['surrounding_surface']['mean_displacement_z_m']}
    comparison=read(run/'independent-3d/independent-local-verification.json')
    assert comparison['solver_sha256']==digest
    assert comparison['all_active_nodes_three_components_relative_error']<.005
    assert comparison['maximum_normalized_yield']<-1e-5
    assert comparison['reference_refinement_error']<1e-4
    for path,expected in comparison['input_and_reference_sha256'].items():
        assert sha(Path(path))==expected,path
    history=read(run/'plastic-history/checkpoint-verification.json')
    assert history['solver_sha256']==digest
    for key in ('separate_process_resume','exact_final_committed_state','exact_frames','exact_increment_path',
                'plastic_history_preserved','cutback_exact_vs_clean_accepted_path',
                'cutback_rejected_changed_material_load_and_active_set','changed_inputs_rejected',
                'corrupted_history_rejected','existing_output_unchanged',
                'rejected_trial_plastic_history_changed'):
        assert history[key],key
    assert {r['stage'] for r in history['restarts']}=={'preload','unload','excavate_pit'}
    assert all(r['exact'] and 0<r['saved_progress']<1 for r in history['restarts'])
    assert history['removed_integration_point_histories_frozen']>0
    for path,expected in history['files_sha256'].items():
        assert sha(run/'plastic-history'/path)==expected,path
    continuous=read(run/'plastic-history/continuous.json')
    clean=read(run/'plastic-history/clean.json')
    resumed=read(run/'plastic-history/resumed.json')
    assert continuous['last_committed']==clean['last_committed']==resumed['last_committed']
    assert continuous['frames']==resumed['frames']
    assert [f['state'] for f in continuous['frames']]==[f['state'] for f in clean['frames']]
    for name in ('changed-input','tampered','no-overwrite'):
        assert (run/f'plastic-history/{name}.log').read_text().rstrip().endswith('Exit code: 1')
    for name in ('cad-tests','flat','curved','independent-3d','plastic-history'):
        assert (run/(name+'.log')).read_text(encoding='utf-8').rstrip().endswith('Exit code: 0')
    regression=run/'mcc-regression.json';legacy_files=[]
    if regression.exists():
        compatibility=read(regression)
        assert compatibility['exact_committed_state'] and compatibility['exact_frames']
        assert compatibility['current_solver_sha256']==digest
        for path,expected in compatibility['files_sha256'].items():
            path=Path(path);assert sha(path)==expected;legacy_files.append(path)
    # Freeze before building/testing, then reject any source change during the
    # suite. An archive captured only at the end cannot bind tested sources.
    source_manifest=read(run/'source-snapshot.json')
    assert source_manifest==manifest(),'sources changed during acceptance'
    archive=run/'source-snapshot.zip'
    with zipfile.ZipFile(archive) as zipped:
        assert set(zipped.namelist())==set(source_manifest)
        assert zipped.testzip() is None
        for name,expected in source_manifest.items():
            assert hashlib.sha256(zipped.read(name)).hexdigest()==expected,name
    files=[p for p in run.rglob('*') if p.is_file()]+list(checks.values())+legacy_files
    report={
        'schema':'soil-fem-local-excavation-final-validation/1','passed':True,
        'solver_sha256':digest,'rust_tests':count,'cad_tests':5,'examples':examples,
        'independent_3d_elastic_branch':{
            'all_active_nodes_three_components_relative_error':comparison['all_active_nodes_three_components_relative_error'],
            'reference_refinement_error':comparison['reference_refinement_error'],
            'maximum_normalized_yield':comparison['maximum_normalized_yield']},
        'plastic_cutback_exact_vs_clean':True,'separate_process_restart_exact':True,
        'restart_phases':[r['stage'] for r in history['restarts']],
        'rejected_trial_plastic_history_changed':True,
        'removed_integration_point_histories_frozen':history['removed_integration_point_histories_frozen'],
        'legacy_mcc_exact':regression.exists(),
        'limitations':'Small strain, fixed pore pressure. Independent global 3D comparison covers the elastic MCC unloading branch at pc/p=10; plastic local 3D displacement is not independently certified. Spatial convergence and P6 overall are not established.',
        'source_snapshot_sha256':sha(archive),
        'sha256':{str(p.relative_to(root)):sha(p) for p in sorted(set(files))}}
    with target.open('x',encoding='utf-8') as stream:
        json.dump(report,stream,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--run',type=Path,required=True)
    for name in ('rust','clippy','format','build'):
        parser.add_argument('--'+name+'-log',type=Path)
    args=parser.parse_args();run=args.run.resolve()
    defaults={'rust':'rust-tests.log','clippy':'lint.log','format':'format.log','build':'build.log'}
    checks={name:(getattr(args,name+'_log') or run/path).resolve() for name,path in defaults.items()}
    record(run,checks)
