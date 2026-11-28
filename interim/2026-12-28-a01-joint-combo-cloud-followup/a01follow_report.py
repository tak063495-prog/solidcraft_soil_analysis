"""Fixed-order P/R comparisons; no interpolation, extrapolation or pass gate."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
from a01follow_binding import load, compare, norm, budget, fixtures
from a01follow_run import write

def read(path):return json.loads(path.read_text(encoding='utf-8'))

def trace_summary(folder,engine):
    path=folder/('native/progress.jsonl' if engine=='native' else 'trace.jsonl')
    counts={};residuals=[];partition_events=[]
    for line in path.read_text(encoding='utf-8').splitlines():
        v=json.loads(line);kind=v['kind'];counts[kind]=counts.get(kind,0)+1
        if kind=='iterate':residuals.append(dict(stage=v['stage'],progress=v['progress'],iteration=v['iteration'],relative=v.get('relative',v.get('relative_force_residual'))))
        if kind=='repartition':partition_events.append(v)
    return dict(event_counts=counts,residual_history=residuals,repartition_events=partition_events,
        trace_sha256=hashlib.sha256(path.read_bytes()).hexdigest())

def main():
    p=argparse.ArgumentParser();p.add_argument('--formal',type=Path,required=True);p.add_argument('--extra',type=Path,required=True);p.add_argument('--separation',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False);cases={};paths=[]
    for case in ('shared-stress_free','shared-transfer_force','joint-plate'):
        def folder(engine,pred,part,boost=False):return a.formal/(f'{case}__{engine}__P{pred}__R{part}__g8'+('__boost' if boost else ''))
        values={(engine,pred,part):load(folder(engine,pred,part)/'fields.jsonl') for engine in ('native','independent') for pred in ('native','independent') for part in ('native','independent')}
        baseline=values['native','native','native'];compared={};histories={};paired={}
        for key,value in values.items():
            name='__'.join(key);compared[name]=compare(baseline,value)
            histories[name]=trace_summary(folder(*key),key[0])
            result=read(folder(*key)/'result.json');paths.append(dict(case=case,engine=key[0],predictor=key[1],repartition=key[2],complete=result['complete'],accepted_states=result['accepted_states']))
        for pred in ('native','independent'):
            for part in ('native','independent'):
                paired[f'P{pred}__R{part}']=compare(values['native',pred,part],values['independent',pred,part])
        boosted=load(folder('independent','independent','independent',True)/'fields.jsonl')
        olddiff=baseline[1]['endpoint_work']-boosted[1]['endpoint_work']
        remaining=baseline[1]['endpoint_work']-values['independent','native','native'][1]['endpoint_work']
        factorial={}
        for engine in ('native','independent'):
            w=lambda p,r:values[engine,p,r][1]['endpoint_work']
            p_effect=.5*(w('native','native')-w('independent','native')+w('native','independent')-w('independent','independent'))
            r_effect=.5*(w('native','native')-w('native','independent')+w('independent','native')-w('independent','independent'))
            interaction=w('native','native')-w('native','independent')-w('independent','native')+w('independent','independent')
            factorial[engine]=dict(predictor_main_effect_l2=norm(p_effect),partition_main_effect_l2=norm(r_effect),interaction_l2=norm(interaction))
        cases[case]=dict(paired_work_relative_l2={k:v['all_fields']['endpoint_work']['relative_l2'] for k,v in paired.items()},
            paired_open_mismatches={k:v['all_fields']['contact_open']['mismatches'] for k,v in paired.items()},
            old_partial_boost_comparison=compare(baseline,boosted)['all_fields']['endpoint_work'],
            unexplained_work_relative_to_boost_reference=norm(remaining)/max(norm(boosted[1]['endpoint_work']),1e-12),
            residual_vector_budget=budget(olddiff,remaining),factorial=factorial)
        write(a.output/(case+'-all-comparisons.json'),dict(compared_to_native_default=compared,paired=paired,histories=histories))
    stacked=[]
    for root in (a.formal,a.extra):
        for folder in sorted(root.glob('stacked__*')):
            if not folder.is_dir():continue
            r=read(folder/'result.json');stacked.append(dict(name=folder.name,**{k:v for k,v in r.items() if k not in ('rows','traceback')}))
    mode=read(a.separation/'mode-separation.json');descent=read(a.separation/'descent-decomposition.json')
    summary=dict(passed=False,a01_complete=False,diagnostic_only=True,new_limited_acceptance=False,
        cases=cases,exchange_paths=paths,stacked_failures=stacked,
        mode=dict(branch_gp_indices=mode['branch_gp_indices'],
            source_smallest=mode['source_pre_predictor_smallest_singular_values'][-2:],
            rigid_contact_rank={k:v['middle_rigid_contact_rank'] for k,v in mode['linearizations'].items()},
            physical_instability_proved=False),
        descent={r['label']:dict(actual_merit_slope=r['source_actual_merit_slope'],
            tangent_defect_merit_slope=r['tangent_defect_merit_slope'],initial_projected_force_norm=r['initial_projected_force_norm'],
            birth_reference_norm=r['birth_reference_norm'],holding_force_norm=r['holding_force_norm'],
            extra_grid=r['outside_candidate_grid_observations']) for r in descent['records']},
        acceptance_tolerance_changed=False,old_values_rewritten=False,cloud_exhausted=False,
        a02_a04_b01_b03_c01_c04_nc01_p6_complete=False)
    write(a.output/'measurements.json',summary)

if __name__=='__main__':main()
