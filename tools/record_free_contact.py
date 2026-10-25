"""Seal the cold-contact experiment only after checking raw results and sources."""
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
    assert not target.exists()
    assert manifest()==read(run/'source-snapshot.json'),'source changed after pre-build freeze'
    digest=sha(run/'solver/soil-fem-3d.exe');suite=read(run/'verification.json')
    assert suite['passed'] and suite['solver_sha256']==digest
    reused=(run/'rust-test-reuse.json').exists()
    if reused:
        from reuse_rust_test_evidence import sha as reused_sha
        reuse=read(run/'rust-test-reuse.json')
        assert reuse['passed'] and reuse['expected_solver_sha256']==digest
        assert sha(run/'reused-rust-check/source-snapshot.zip')==reuse['previous_source_snapshot_sha256']
        for path,expected in reuse['original_check_logs_sha256'].items():assert reused_sha(Path(path))==expected
        sources=read(run/'source-snapshot.json')
        assert all(sources[n]==h for n,h in reuse['identical_compilation_and_test_inputs'].items())
    for name in ('freeze','format','rust-tests','lint','build','cad-tests','split-tests','cad','free-contact','mc','mcc','compatibility',
                 'local-flat','local-flat-regression','local-curved','local-curved-regression'):
        assert (run/(name+'.log')).read_text(encoding='utf-8').rstrip().endswith('Exit code: 0'),name
    rust=(run/'rust-tests.log').read_text(encoding='utf-8')
    count=sum(map(int,re.findall(r'test result: ok\. (\d+) passed;',rust)));assert count>=121
    for name in ('free_upper_block_closes_without_physical_pins_and_restores_predictor_history',
                 'free_frictionless_or_overcapacity_blocks_never_commit_a_false_equilibrium',
                 'cold_elastic_prediction_matches_a_closed_sticking_physical_increment',
                 'temporary_sparse_and_dense_gauges_do_not_certify_the_original_equations',
                 'local_cutback_restores_plastic_history_loads_and_active_elements'):
        assert f'{name} ... ok' in rust,name
    assert 'Ran 5 tests' in (run/'cad-tests.log').read_text(encoding='utf-8')
    assert 'Ran 2 tests' in (run/'split-tests.log').read_text(encoding='utf-8')
    free=read(run/'free/verification.json');assert free['passed'] and free['solver_sha256']==digest
    assert set(free['cases'])=={'horizontal','inclined','near-cap','stiff','mcc-oc','mcc-nc','cad'}
    assert free['separate_process_exact_restarts']==20
    assert set(free['failures'])=={'frictionless','overcapacity','critical','opening'}
    assert free['cases']['mcc-nc']['plastic_bulk_rejected_trial_changed']
    assert free['sparse_dense']['same_accepted_path']
    for name in free['cases']:
        folder=run/'cad/analysis' if name=='cad' else run/'free'/name
        data=read(folder/'stage-results.json');audit=read(folder/'stage-verification.json')
        assert data['complete'] and data['provenance']['solver_binary_sha256']==digest
        assert all(f['joint_support']['unrestrained_modes']==0 for f in data['frames'])
        assert all(c['node_set']=='bottom' for c in data['case']['boundary']['constraints'])
        assert all(a['committed_state_unchanged_on_rejection'] for a in data['geostatic']['attempts']+data['attempts'] if not a['accepted'])
        assert max(f['independent_force_residual'] for f in audit['frames'])<1e-8
        assert audit['joint_rigid_support_svd_audit']=='passed'
        clean=read(folder/'clean/result.json')
        assert clean['last_committed']==data['last_committed'] and clean['frames']==data['frames']
        assert data['frames'][-2]['state']['joints']==data['frames'][-1]['state']['joints']
        assert data['excavation_audits'][0]['zero_release_relative_force_residual']<1e-8
    for name in ('horizontal','inclined','near-cap','mcc-nc','cad'):
        cp=read(run/'free'/(name+'-checkpoint')/'verification.json')
        assert cp['solver_sha256']==digest and len(cp['restarts'])==4
        assert all(r['exact_states_frames_material_ids_and_audits'] for r in cp['restarts'])
        for path,expected in cp['sha256'].items():assert sha(run/'free'/(name+'-checkpoint')/path)==expected
    for name,comparison in free['comparisons'].items():
        assert comparison['reference_subincrements']==[100,200] and comparison['maximum_independent_force_residual']<1e-8
        assert all(f['relative_all_node_3d_displacement_error']<.005 and f['reference_refinement_error']<1e-4 for f in comparison['frames'])
    for name in ('flat','curved'):
        comparison=read(run/f'local-{name}-regression.json')
        assert comparison['exact_committed_state'] and comparison['exact_frames'] and comparison['current_solver_sha256']==digest
        for path,expected in comparison['files_sha256'].items():assert sha(Path(path))==expected
        same_input=read(run/f'local-{name}-same-input/verification.json')
        assert same_input['passed'] and same_input['exact_states_and_frames']
        for path,expected in same_input['inputs_sha256'].items():assert sha(Path(path))==expected
    compatibility=read(run/'contact-state-compatibility.json');assert compatibility['passed']
    for path,expected in compatibility['sha256'].items():assert sha(Path(path))==expected
    for item in compatibility['comparisons']:
        assert item['exact_committed_state'] and item['exact_frames_except_added_support_report'] and item['current_solver_sha256']==digest
    sources=read(run/'source-snapshot.json')
    report={'schema':'soil-fem-free-contact-final-validation/1','passed':True,'solver_sha256':digest,
        'rust_tests':count,'cad_and_split_tests':7,'separate_process_exact_restarts':20,
        'source_freeze_before_build':True,'frozen_source_files':len(sources),
        'rust_test_evidence_reused_for_identical_sources':reused,
        'cases':free['cases'],'independent_comparisons':free['comparisons'],'sparse_dense':free['sparse_dense'],
        'physical_failure_gates':free['failures'],'local_box_and_cylinder_mcc_exact_regression':True,
        'warm_mc_mcc_exact_regression':True,
        'limitations':free['reference_scope']+' '+free['nc_loading_scope']+' Dry fixed normals, small strain; no wet contact, crack fronts or interface BIRTH. P6 remains incomplete.',
        'source_snapshot_sha256':sha(run/'source-snapshot.zip'),
        'sha256':{p.relative_to(root).as_posix():sha(p) for p in sorted(run.rglob('*')) if p.is_file()}}
    with target.open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k not in ('sha256','independent_comparisons')},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();record(a.run.resolve())
