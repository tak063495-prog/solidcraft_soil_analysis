"""Evidence identities and complete accepted/state field inventories; no tolerance updates."""
import argparse,json,copy
from pathlib import Path
from b03_cases import read,write
from study import exact


def verify(formal,scope,supplemental,mutations,output):
    measure=read(formal/'measurements.json');restarts=read(formal/'restarts.json');checks=[]
    for c in measure['cases']:
        f=formal/c['case'];native=[json.loads(s) for s in (f/'native/accepted.jsonl').read_text().splitlines()]
        independent=[json.loads(s) for s in (f/'independent.jsonl').read_text().splitlines()]
        steps=c['load_steps'];expected=[(i,j/steps) for i in range(5) for j in range(1,steps+1)]
        assert [(r['stage_index'],r['progress']) for r in native]==expected
        assert [(r['stage_index'],r['progress']) for r in independent]==expected
        assert [r['accepted_index'] for r in native]==list(range(1,5*steps+1))
        comparison=read(f/'comparison.json');assert len(comparison['all_point_component_comparisons'])==len(expected)
        run=read(f/'native/run.json');assert exact(run,read(f/'plain/run.json'))
        endpoints=[r for r in native if r['progress']==1.]
        assert len(endpoints)==len(run['frames'])==5
        for index,(frame,observed) in enumerate(zip(run['frames'],endpoints,strict=True)):
            assert exact(frame['beam_results'],observed['beam_results'])
            state=copy.deepcopy(observed['state'])
            if index==0:
                assert state['stage_progress']==1. and frame['state']['stage_progress']==0.
                state['stage_progress']=0.
                release=state.pop('geostatic_release_load_kn')
                assert len(release)==3*len(observed['state']['displacement_m']) and all(v==0. for v in release)
                assert 'geostatic_release_load_kn' not in frame['state']
                # stages.rs:5132-5133 resets progress and clears this zero vector.
            assert exact(frame['state'],state)
        m=read(f/'mesh.json');case=read(f/'stages.json');nplate=len(case['beam_system']['plates'])
        for n,p in zip(native,independent,strict=True):
            assert len(n['state']['material_states'])==27*len(m['tetrahedra'])==len(p['stress_kpa'])
            assert len(n['state']['displacement_m'])==len(m['nodes'])
            assert len(n['beam_results'].get('plates',[]))==sum(n['state']['beams']['active_plates'])
            assert all(len(v['samples'])==6 for v in n['beam_results'].get('plates',[]))
            assert len(n['state']['beams']['active_plates'])==nplate
        checks.append(dict(case=c['case'],accepted_points=len(expected),soil_gp_states=len(expected)*27*len(m['tetrahedra']),
                           all_point_inventory_verified=True,all_frame_and_plain_run_verified=True,
                           actual_saved_frame_recovery_matches_observation=True,
                           initial_frame_progress_transition='stages.rs:5132-5133: geostatic endpoint progress 1.0 -> frame 0.0 and exact-zero release vector cleared; every other field exact; restart comparisons exclude neither'))
    assert all(r['passed'] for r in restarts)
    s=read(scope);supp=read(supplemental);mut=read(mutations)
    mismatched=[r['name'] for r in s['records'] if r['expected_diagnostic'] and not r['refused_as_expected']]
    assert set(mismatched)=={'curved-midpoint-in-plane','curved-midpoint-out-of-plane'}
    for r in s['records']:
        if r['expected_diagnostic']:
            assert r['process']['exit_code']==1 and not r['result_file_exists']
            if r['name'] in mismatched:assert 'geostatic mesh/CAD volume difference' in r['diagnostic']
        else:assert r['accepted']
    assert supp['refusal_contract_passed'] and supp['legacy_controls_completed'] and mut['refusal_contract_passed']
    assert read(formal/'source-binding-before.json')==read(formal/'source-binding-after.json')
    result=dict(passed=False,b03_complete=False,evidence_contract_checks_passed=True,cases=checks,
        total_accepted_points=sum(r['accepted_points'] for r in checks),soil_gp_states_compared=sum(r['soil_gp_states'] for r in checks),
        restarts=len(restarts),later_states_compared=sum(r['later_accepted_states'] for r in restarts),
        all_runs_and_later_states_exact=True,original_probe_expectation_mismatches_retained=mismatched,
        isolated_midpoint_rejections_confirmed=True,all_content_mutations_refused=True,
        statement='Evidence completeness and refusal audit only. No B-03 full numerical acceptance or new all-field tolerance.')
    write(output,result)
if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ('formal','scope','supplemental','mutations','output'):p.add_argument(k,type=Path)
    a=p.parse_args();verify(a.formal,a.scope,a.supplemental,a.mutations,a.output)
