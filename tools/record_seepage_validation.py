"""Seal steady-flow acceptance only against pre-build frozen sources."""
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
    assert not target.exists() and read(run/'source-snapshot.json')==manifest(),'sources changed during acceptance'
    digest=sha(run/'solver/soil-fem-3d.exe');suite=read(run/'verification.json');assert suite['passed'] and suite['solver_sha256']==digest
    columns=('mc-tet4','mc-tet10','mcc-tet10');cad=('cad-mc','cad-mcc','cad-curved-mc')
    logs=['freeze','format','rust-tests','lint','build','cad-tests','split-tests','darcy-cli','legacy','mc-tet4-theory','mc-tet10-theory']
    logs += [name+tail for name in columns for tail in ('','-history','-checkpoint')]
    logs += list(cad)+['cad-mc-independent','cad-curved-mc-independent','cad-mcc-history','cad-mcc-checkpoint']
    for name in logs:assert (run/(name+'.log')).read_text(encoding='utf-8').rstrip().endswith('Exit code: 0'),name
    rust=(run/'rust-tests.log').read_text(encoding='utf-8');count=sum(map(int,re.findall(r'test result: ok\. (\d+) passed;',rust)));assert count>=151
    for name in ('original_rows_are_refined_after_prescribed_head_elimination_on_a_cad_mesh','dense_and_sparse_backends_match_the_original_source_problem',
        'mcc_seepage_cutback_and_each_phase_checkpoint_restore_pressure_and_plastic_history','unsaturated_candidate_does_not_change_the_committed_mechanical_state',
        'balanced_neumann_and_volume_exchange_use_absolute_flow_for_water_balance'):
        assert name+' ... ok' in rust,name
    assert 'Ran 5 tests' in (run/'cad-tests.log').read_text(encoding='utf-8')
    assert 'Ran 2 tests' in (run/'split-tests.log').read_text(encoding='utf-8')
    assert read(run/'darcy-cli/verification.json')['passed'] and read(run/'darcy-cli/verification.json')['absolute_flow_normalization_exercised']
    assert read(run/'darcy-cli/verification.json')['independent_domains_and_datum_exercised']
    rows=[];maximum_head_error=0.;maximum_balance=0.;refinements=0
    for name in columns+cad:
        folder=run/name;data=read(folder/'stage-results.json');assert data['complete'] and data['schema']=='soil-fem-stage-result/12'
        assert data['provenance']['solver_binary_sha256']==digest
        audit=read(folder/'stage-verification.json');assert audit['raw_integration_point_roundtrip']=='passed'
        assert max(f['independent_force_residual'] for f in audit['frames'])<1e-8
        for row in audit['independent_darcy_solutions']:maximum_head_error=max(maximum_head_error,row['maximum_head_error_m'])
        for flow in data['flow_audits']:
            solution=flow['solution'];assert solution['complete'] and solution['mesh_sha256']==data['provenance']['mesh_sha256']
            assert solution['relative_free_head_residual']<=1.01e-11
            for balance in solution['balances']:
                assert balance['relative_free_head_residual']<=1.01e-11
                maximum_balance=max(maximum_balance,balance['relative_nodal_balance_error']);assert balance['relative_nodal_balance_error']<=2e-10
            refinements+=len(solution.get('refinement_solves',[]))
        rows.append({'case':name,'active_elements':sum(data['last_committed']['active_elements']),'flow_stages':len(data['flow_audits']),
            'final_monitors':data['frames'][-1]['monitors'],'maximum_independent_force_residual':max(f['independent_force_residual'] for f in audit['frames'])})
        if name in cad:
            assert sha(folder/'solver/soil-fem-3d.exe')==digest
            assert len(read(folder/'ground/bodies.json')['excavation_tools'])==1
    assert refinements>0
    for name in ('mc-tet4','mc-tet10'):assert read(run/name/'column-theory.json')['passed']
    restarts=0
    for name in columns+('cad-mcc',):
        history=read(run/(name+'-history')/'verification.json');assert history['passed'] and history['exact_accepted_path_and_all_states']
        if name in ('mcc-tet10','cad-mcc'):assert history['rejected_plastic_trial_exercised']
        saved=read(run/(name+'-checkpoint')/'verification.json');assert saved['solver_sha256']==digest
        assert len(saved['restarts'])==5 and all(r['exact_states_frames_material_ids_and_audits'] for r in saved['restarts']);restarts+=len(saved['restarts'])
    linear=[]
    for name in ('cad-mc','cad-curved-mc'):
        data=read(run/(name+'-independent')/'verification.json');assert data['passed'] and data['solver_sha256']==digest
        linear.append({'case':name,'frames':data['frames']})
    assert read(run/'legacy/seepage-legacy-verification.json')['passed']
    files=[p for p in run.rglob('*') if p.is_file()]
    report={'schema':'soil-fem-seepage-final-validation/1','passed':True,'solver_sha256':digest,'rust_tests':count,'cad_tests':5,'joint_split_tests':2,
        'cases':rows,'separate_process_restarts':restarts,'maximum_independent_head_error_m':maximum_head_error,
        'maximum_component_nodal_water_balance_error':maximum_balance,'original_row_refinement_solves':refinements,
        'independent_3d_linear_mechanics':linear,'exact_legacy_cases':6,
        'scope':'Steady saturated hydraulic-head FEM and explicit ONE_WAY small-strain mechanics. No unknown phreatic surface, seepage-face complementarity, Richards time, local conservative flux reconstruction, or general independent plastic 3D displacement validation is inferred.',
        'source_snapshot_sha256':sha(run/'source-snapshot.zip'),'sha256':{str(p.relative_to(root)):sha(p) for p in sorted(files)}}
    with target.open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k not in ('sha256','independent_3d_linear_mechanics')},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);record(p.parse_args().run.resolve())
