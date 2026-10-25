"""Seal IDH material, global construction and actual SolidCraft evidence."""
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
    assert not target.exists() and manifest()==read(run/'source-snapshot.json'),'sources changed since freeze'
    digest=sha(run/'solver/soil-fem-3d.exe');suite=read(run/'verification.json')
    assert suite['passed'] and suite['solver_sha256']==digest
    reused=(run/'rust-test-reuse.json').exists()
    if reused:
        reuse=read(run/'rust-test-reuse.json');assert reuse['passed'] and reuse['expected_solver_sha256']==digest
        assert sha(run/'reused-rust-check/source-snapshot.zip')==reuse['previous_source_snapshot_sha256']
        for path,expected in reuse['original_check_logs_sha256'].items():assert sha(Path(path))==expected
        sources=read(run/'source-snapshot.json')
        assert all(sources[n]==h for n,h in reuse['identical_compilation_and_test_inputs'].items())
    names=['freeze','format','rust-tests','lint','build','cad-tests','split-tests','elastic-cycle','material','triaxial',
        'tet4','tet4-history','tet4-checkpoint','tet10','tet10-history','tet10-checkpoint','column-reference','local','local-history','local-checkpoint','legacy']
    for name in names:assert (run/(name+'.log')).read_text(encoding='utf-8').rstrip().endswith('Exit code: 0'),name
    rust=(run/'rust-tests.log').read_text(encoding='utf-8');count=sum(map(int,re.findall(r'test result: ok\. (\d+) passed;',rust)))
    assert count>=130
    for name in ('closed_nonproportional_elastic_cycles_leave_no_residual_stress_or_work',
        'a_dilative_return_can_recover_a_negative_elastic_predictor_without_pressure_clipping',
        'idh_preload_unload_reload_and_excavation_satisfy_current_force_and_sparse_dense_gates',
        'idh_cap_shear_void_ratio_and_release_histories_restore_after_cutback_and_all_phase_restarts',
        'local_cutback_restores_plastic_history_loads_and_active_elements'):
        assert name+' ... ok' in rust,name
    assert 'Ran 5 tests' in (run/'cad-tests.log').read_text(encoding='utf-8')
    assert 'Ran 2 tests' in (run/'split-tests.log').read_text(encoding='utf-8')
    reports={}
    for name in ('elastic-cycle','material','triaxial','column-reference','legacy'):
        d=read(run/name/'verification.json');assert d['passed'] and d['solver_sha256']==digest
        reports[name]=d
    assert len(reports['material']['cases'])==6 and reports['material']['native_return_types']==['cap','elastic','shear','shear_cap']
    assert reports['material']['excessive_expansion_domain_certificate']['maximum_possible_final_mean_pressure_kpa']<0
    assert len(reports['triaxial']['cases'])==3 and reports['triaxial']['reference_loading_increments']==[200,400]
    column=reports['column-reference'];assert column['reference_subincrements']==[64,128,256] and column['maximum_reference_force_residual']<1e-8
    cases={};restarts=0
    for name in ('tet4','tet10','local'):
        folder=run/name;data=read(folder/'stage-results.json');audit=read(folder/'stage-verification.json');accepted=read(folder/'accepted-stage-run.json')
        assert data['complete'] and data['schema']=='soil-fem-stage-result/10' and data['provenance']['solver_binary_sha256']==digest
        assert accepted['solver_binary_sha256']==digest
        for path,expected in accepted['inputs_and_results_sha256'].items():assert sha(Path(path))==expected
        assert audit['vtu_roundtrip']=='passed' and audit['raw_integration_point_roundtrip']=='passed'
        assert max(f['independent_force_residual'] for f in audit['frames'])<=data['case']['settings']['relative_force_tolerance']
        assert all(a['zero_release_relative_force_residual']<1e-8 for a in data['excavation_audits'])
        assert not any(data['last_committed']['excavation_release_load_kn'])
        history=read(run/(name+'-history')/'verification.json');assert history['passed'] and history['solver_sha256']==digest
        assert history['exact_accepted_path_states_frames_and_audits'] and history['all_rejected_trials_preserved_committed_state']
        cp=read(run/(name+'-checkpoint')/'verification.json');assert cp['solver_sha256']==digest
        assert len(cp['restarts'])==(3 if name=='local' else 4)
        assert any(r['saved_nonzero_excavation_holding'] for r in cp['restarts'])
        assert all(r['exact_states_frames_material_ids_and_audits'] for r in cp['restarts'])
        for path,expected in cp['sha256'].items():assert sha(run/(name+'-checkpoint')/path)==expected
        restarts+=len(cp['restarts'])
        final=data['frames'][-1];cases[name]={'accepted_increments':sum(a['accepted'] for a in data['attempts']),
            'rejected_increments':sum(not a['accepted'] for a in data['attempts']),
            'maximum_independent_force_residual':max(f['independent_force_residual'] for f in audit['frames']),
            'excavation_audits':data['excavation_audits'],'final_monitors':final['monitors'],
            'maximum_absolute_engineering_strain':max(f['maximum_absolute_engineering_strain_from_k0'] for f in data['frames'])}
        if name=='local':
            assert final['monitors']['surrounding_surface']['mean_displacement_z_m']<0
            assert final['monitors']['excavation_floor']['stage_increment_z_m']>0
            assert sha(folder/'solver/soil-fem-3d.exe')==digest
    assert restarts==11
    for shape in ('flat','curved'):
        d=read(run/'legacy'/(shape+'-comparison.json'));assert d['current_solver_sha256']==digest and d['exact_committed_state'] and d['exact_frames']
        for path,expected in d['files_sha256'].items():assert sha(Path(path))==expected
    report={'schema':'soil-fem-hardening-final-validation/1','passed':True,'solver_sha256':digest,
        'rust_tests':count,'cad_and_split_tests':7,'separate_process_exact_restarts':restarts,
        'rust_test_evidence_reused_for_identical_sources':reused,
        'source_freeze_before_build':True,'frozen_source_files':len(read(run/'source-snapshot.json')),
        'cases':cases,'independent_column_comparison':column,'drained_triaxial':reports['triaxial'],
        'rejected_elastic_formulation_preserved_and_fixed':True,'mcc_local_box_and_cylinder_exact_regression':True,
        'limitations':'Own invariant DP shear/elliptic cap model, distinct from HS/PH/CYSoil. Independent global plastic displacement comparison is a 3D mesh with oedometer plane constraints; general unconstrained 3D plastic fields are not independently certified. Local cut uses illustrative cohesive c*=30 kPa, dry small strain. No field calibration, variable dilatancy, liquefaction or coupled water. P6 remains incomplete.',
        'source_snapshot_sha256':sha(run/'source-snapshot.zip'),
        'sha256':{p.relative_to(root).as_posix():sha(p) for p in sorted(run.rglob('*')) if p.is_file()}}
    with target.open('x',encoding='utf-8') as f:json.dump(report,f,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k not in ('sha256','independent_column_comparison','drained_triaxial')},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();record(a.run.resolve())
