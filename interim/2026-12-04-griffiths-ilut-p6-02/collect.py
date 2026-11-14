"""Remeasure saved raw fields. This command performs no FEM solves."""
import argparse
import copy
import json
from pathlib import Path
from study import bits, compare, fields, read, sha, summary, write


def load_lines(p):
    with p.open(encoding='utf-8') as f:return [json.loads(x) for x in f]


def collect(raw,out):
    data=read(raw/'final/measurements.json')
    for case in data['cases']:
        folder=raw/'final'/f"stride{case['stride']}"
        runs={k:read(folder/(k+'.json')) for k in ['base','default','ilut']}
        assert {k:summary(v) for k,v in runs.items()}==case['summaries']
        assert compare(runs['default'],runs['ilut'])==case['field_differences']
        before=copy.deepcopy(runs['base']);after=copy.deepcopy(runs['default'])
        assert before['provenance']['solver_binary_sha256']==sha(raw/'base/soil-fem-3d')
        assert after['provenance']['solver_binary_sha256']==sha(raw/'frozen/soil-fem-3d')
        before['provenance']['solver_binary_sha256']=after['provenance']['solver_binary_sha256']
        assert bits(before)==bits(after) and case['default_unchanged'] is True
        for r in runs.values():
            assert r['complete'] and r['geostatic']['complete']
            attempts=r['attempts']+r['geostatic']['attempts']
            assert all(a['accepted'] and a['relative_force_residual']<=1e-8 for a in attempts)
            assert all(s['true_relative_residual']<=1e-10 for a in attempts for s in a['linear_solves'])
    restart=read(raw/'restarts/restarts.json')
    direct=read(raw/'restarts/ilut/run.json');trace=load_lines(raw/'restarts/ilut/accepted.jsonl')
    for row in restart['restarts']:
        n=int(Path(row['checkpoint']).stem.split('-')[-1]);folder=raw/'restarts'/f'resume-{n}'
        assert row['separate_pid'] and row['all_later_states_equal'] and row['whole_run_equal']
        assert bits(read(folder/'run.json'))==bits(direct)
        assert bits(load_lines(folder/'accepted.jsonl'))==bits([x for x in trace if x['accepted_index']>n])
    assert restart['ordinary_cli']['separate_pid'] and restart['ordinary_cli']['whole_run_equal']
    assert bits(read(raw/'restarts/cli-resumed.json'))==bits(read(raw/'final/stride5/ilut.json'))
    assert all(r['rejected_as_expected'] for r in restart['input_rejections']+restart['mutations'])
    assert restart['default_only']==restart['ilut_only']==0
    left=load_lines(raw/'restarts/default/accepted.jsonl')
    expected=[dict(stage_index=a['stage_index'],progress=a['progress'],**fields(a['state'],b['state']))
              for a,b in zip(left,trace,strict=True)]
    assert expected==restart['full_state_field_differences']
    bindings=[]
    for mode,file in [('default','default'),('ilut','ilut')]:
        a=read(raw/'restarts'/mode/'run.json');b=read(raw/'final/stride5'/(file+'.json'))
        assert a['provenance']['solver_binary_sha256']==sha(raw/'frozen/p602-observer')
        assert b['provenance']['solver_binary_sha256']==sha(raw/'frozen/soil-fem-3d')
        a['provenance']['solver_binary_sha256']=b['provenance']['solver_binary_sha256']
        bindings.append(dict(mode=mode,all_run_except_clock_and_authenticated_binary_sha_equal=bits(a)==bits(b)))
    assert all(x['all_run_except_clock_and_authenticated_binary_sha_equal'] for x in bindings)
    data.update(p602_decision=dict(limited_acceptance=['追加オプションの入力/保存/拒否契約',
        '48/300要素の既存数値受理条件と計測整合', '小型別PIDの全状態/全Run一致'],
        unmet=['元5%基準','正確な斜面Fs','実断面SRM','xy2','反復上限付近の改善']),
        restart_summary=dict(checkpoints=len(restart['restarts']),later_states=sum(r['later_states'] for r in restart['restarts']),
            ordinary_cli=restart['ordinary_cli'],input_rejections=len(restart['input_rejections']),
            mutations=len(restart['mutations']),accepted_states_per_preconditioner=len(trace)),
        observer_to_cli_binding=bindings,measurements_reproduced_from_raw=True)
    data['user_pc_report_unreproduced_in_cloud']=dict(elements=4800,geostatic_progress=.4717,
        relative_force_residual=1.109e-7,force_limit=1e-8,gmres_limit=3000,newton_limit=35,
        gmres_limit_reported_from_progress_approximately=.44,
        subdivision_residuals_reported=[1.33e-8,3.47e-8,1.11e-7],source='user instruction; no xy2 Run generated here')
    data['software_contracts_verified_in_this_small_scope']=True
    data['small_case_all_accepted_field_difference_maxima']={field:{key:max(x[field][key] for x in expected)
        for key in ('maximum_absolute','relative_l2')} for field in ('displacement_m','stress_kpa')}
    out.mkdir(parents=True,exist_ok=False)
    write(out/'measurements.json',data);write(out/'restarts.json',restart)
    return data


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('raw',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();collect(a.raw,a.output)
