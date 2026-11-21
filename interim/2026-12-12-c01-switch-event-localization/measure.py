"""Compare actual shared input coordinates; never extrapolate or fit a tolerance."""
import argparse
from collections import Counter
from functools import lru_cache
from decimal import Decimal, localcontext
from math import fsum
import importlib.util
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('c01_global_study',HERE/'study.py')
st=importlib.util.module_from_spec(spec);spec.loader.exec_module(st)

@lru_cache(maxsize=3)
def load(folder):
    rows=[];fields={key:[] for key in st.FIELDS}
    for row,f in st.states(Path(folder)):
        rows.append(row)
        for key in fields: fields[key].append(f[key])
    return rows,{key:np.stack(values) for key,values in fields.items()}

def keys(n):
    return [(s,i/n) for s in st.STAGES for i in range(1,n+1)]

def select(folder,coordinates):
    rows,fields=load(str(folder));index={(r['stage'],r['progress']):i for i,r in enumerate(rows)}
    assert len(index)==len(rows) and set(coordinates)<=set(index), 'missing actual coordinate'
    ids=[index[c] for c in coordinates]
    return np.asarray([rows[i]['status'] for i in ids]),{k:v[ids] for k,v in fields.items()}

def squared_sum(values):
    # Fixed exact summation avoids BLAS/SIMD reduction-order differences on Windows.
    return fsum((float(x) for x in values.flat))

def metric(a,b):
    assert a.shape==b.shape and np.isfinite(a).all() and np.isfinite(b).all()
    d=a-b;d2=squared_sum(d*d);r2=squared_sum(b*b)
    with localcontext() as context:
        context.prec=80
        relative=float(Decimal.from_float(d2).sqrt()/max(Decimal.from_float(r2).sqrt(),Decimal.from_float(1e-12)))
    return dict(difference_squared=d2,reference_squared=r2,pooled_relative_l2=relative,maximum_absolute=float(np.max(abs(d),initial=0)))

def compare(a,b,n):
    coordinates=keys(n);sa,fa=select(a,coordinates);sb,fb=select(b,coordinates)
    category=np.zeros((3,3),dtype=int)
    np.add.at(category,(sa.ravel(),sb.ravel()),1)
    return dict(common_states=len(coordinates),fields={k:metric(fa[k],fb[k]) for k in st.FIELDS},status_mismatches=int((sa!=sb).sum()),status_confusion_open_stick_slip=category.tolist(),stage_fields={s:{k:metric(fa[k][i*n:(i+1)*n],fb[k][i*n:(i+1)*n]) for k in st.FIELDS} for i,s in enumerate(st.STAGES)})

def orders(folders,n):
    values=[select(p,keys(n))[1] for p in folders];result={}
    for key in st.FIELDS:
        a=values[0][key]-values[1][key];b=values[1][key]-values[2][key]
        with localcontext() as context:
            context.prec=80
            na=Decimal.from_float(squared_sum(a*a)).sqrt();nb=Decimal.from_float(squared_sum(b*b)).sqrt()
            order=float((na/nb).ln()/Decimal(2).ln()) if na and nb else None
            cosine=float(Decimal.from_float(squared_sum(a*b))/(na*nb)) if na and nb else None
        result[key]=dict(coarse_absolute_norm=float(na),fine_absolute_norm=float(nb),observed_order=order,difference_direction_cosine=cosine)

    return dict(common_states=7*n,fields=result)

def switches(folder):
    rows,_=load(str(folder));previous=np.ones(32,dtype=int);old_stage=None;progress=0.;events=[]
    for row in rows:
        if row['stage']!=old_stage: progress=0.
        now=np.asarray(row['status']);assert now.shape==(32,) and np.isin(now,[0,1,2]).all()
        for point in np.flatnonzero(now!=previous):
            events.append(dict(stage=row['stage'],point=int(point),from_status=int(previous[point]),to_status=int(now[point]),progress_before=progress,progress_after=row['progress']))
        previous=now;progress=row['progress'];old_stage=row['stage']
    return dict(counts=dict(Counter(f"{x['stage']}/{x['from_status']}->{x['to_status']}" for x in events)),events=events)

def legacy_status(n,events):
    names={'open':0,'stick':1,'slip':2};status=np.ones(32,dtype=int);result={}
    by_coordinate={}
    for event in events[f'native/{n}']['native_endpoint_status_events']:
        by_coordinate.setdefault((event['stage'],event['progress_after']),[]).append(event)
    for coordinate in keys(n):
        for event in by_coordinate.get(coordinate,[]):
            assert status[event['point']]==names[event['from_status']]
            status[event['point']]=names[event['to_status']]
        result[coordinate]=status.copy()
    return result

def measure(out,levels=st.STEPS):
    pairs={};observed={};matched={};events={};audits={};estimators={}
    for kind in ('native','independent'):
        for n in levels:
            folder=out/f'g{n}'/kind
            audits[f'{kind}/{n}']=st.read(folder.parent/(kind+'-audit.json'))
            events[f'{kind}/{n}']=switches(folder)
            if kind=='native':
                run=st.read(folder.parent/'native-raw/run.json')
                attempts=run['geostatic']['attempts']+run['attempts']
            else: attempts=st.read(folder.parent/'independent-result.json')['attempts']
            estimated=[a for a in attempts if a.get('joint_path_error') is not None]
            accepted=[a['joint_path_error'] for a in estimated if a['accepted']]
            estimators[f'{kind}/{n}']=dict(attempts_with_estimate=len(estimated),
                accepted_estimates=len(accepted),rejected_estimates=sum(a['joint_path_error']['reject'] for a in estimated),
                midpoint_hidden_switch_point_samples=sum(a['joint_path_error']['hidden_midpoint_switch_points'] for a in estimated),
                coarse_fine_endpoint_disagreement_point_samples=sum(a['joint_path_error']['endpoint_disagreement_points'] for a in estimated),
                maximum_accepted_slip_ratio=max(a['maximum_slip_error_m']/a['slip_budget_m'] for a in accepted),
                maximum_accepted_dissipation_ratio=max(a['maximum_dissipation_error_kpa_m']/a['dissipation_budget_kpa_m'] for a in accepted),
                all_accepted_estimates_within_budget=all(not a['reject'] for a in accepted))
            events[f'{kind}/{n}']['interior_probe_samples']=[dict(stage=a['stage'],progress_from=a['progress_from'],progress_to=a['progress_to'],accepted=a['accepted'],estimate=a['joint_path_error'])
                for a in estimated if a['joint_path_error']['hidden_midpoint_switch_points'] or a['joint_path_error']['endpoint_disagreement_points']]
        for a,b in zip(levels,levels[1:]): pairs[f'{kind}/{a}__{b}']=compare(out/f'g{a}'/kind,out/f'g{b}'/kind,a)
        for a,b,c in zip(levels,levels[1:],levels[2:]): observed[f'{kind}/{a}__{b}__{c}']=orders([out/f'g{n}'/kind for n in (a,b,c)],a)
    for n in levels: matched[str(n)]=compare(out/f'g{n}'/'native',out/f'g{n}'/'independent',n)
    old=st.read(st.ROOT/'interim/2026-11-18-wet-joint-step-convergence/measurements.json')
    old_events=st.read(st.ROOT/'interim/2026-11-18-wet-joint-step-convergence/switch-events.json')
    legacy_statuses={n:legacy_status(n,old_events) for n in levels if n in st.STEPS}
    legacy_native_status_pairs={f'{a}__{b}':sum(int((legacy_statuses[a][k]!=legacy_statuses[b][k]).sum()) for k in keys(a)) for a,b in zip(levels,levels[1:]) if a in legacy_statuses and b in legacy_statuses}
    old_new_status={str(n):int((np.stack([legacy_statuses[n][k] for k in keys(n)])-select(out/f'g{n}/native',keys(n))[0]!=0).sum()) for n in legacy_statuses}
    legacy=dict(native_adjacent_status_mismatches=legacy_native_status_pairs,old_vs_new_native_status_mismatches=old_new_status,source_commit=st.BASE,source='interim/2026-11-18-wet-joint-step-convergence/measurements.json',pairs={k:v['fields']['plastic_slip'] for k,v in old['pairs'].items() if any(k.startswith(kind+'/') for kind in ('native','independent','fixed8'))},observed_orders={k:v['fields']['plastic_slip'] for k,v in old['observed_orders'].items()},switch_event_counts={k:dict(Counter(e['kind'] for e in v['events'])) for k,v in old_events.items()})
    prior_folder=st.ROOT/'interim/2026-12-11-c01-global-step-control'
    prior_report=st.read(prior_folder/'measurements.json')
    prior_events=st.read(prior_folder/'switch-events.json')
    switch_status_differences={}
    for n in levels:
        if n not in st.STEPS: continue
        actual=select(out/f'g{n}/native',keys(n))[0]
        previous=np.ones(32,dtype=int); restored=[]
        for stage in st.STAGES:
            stream=iter(e for e in prior_events[f'native/{n}']['events'] if e['stage']==stage)
            event=next(stream,None)
            for i in range(1,n+1):
                while event is not None and event['progress_after'] <= i/n:
                    assert previous[event['point']]==event['from_status']
                    previous[event['point']]=event['to_status']; event=next(stream,None)
                restored.append(previous.copy())
            assert event is None
        switch_status_differences[str(n)]=int((actual!=np.stack(restored)).sum())
    previous_switch_control=dict(source_commit=st.BASE,source_sha256=st.sha(prior_folder/'measurements.json'),
                                 pairs=prior_report['pairs'],observed_orders=prior_report['observed_orders'],
                                 new_vs_previous_native_status_mismatches=switch_status_differences)
    legacy_all_fields=dict(source='interim/2026-11-18-wet-joint-step-convergence/measurements.json',
                           pairs=old['pairs'],observed_orders=old['observed_orders'])
    return dict(local_error_budget_rule='base_increment_squared_times_scale',local_error_estimator=estimators,previous_switch_control=previous_switch_control,legacy_all_fields=legacy_all_fields,schema='c01-opt-in-global-path-error/1',starting_commit=st.BASE,passed=False,c01_decision='未達',limited_confirmation='fixed 12 Tet10 / 32 contact GP / seven phases; global step doubling estimator, input-only reference and real restart; finite midpoint sampling',global_steps=list(levels),path_error_control=st.PATH_ERROR,force_tolerance=1e-8,true_linear_tolerance=1e-10,independent_force_tolerance=1e-10,physical_model_changed=False,force_tolerances_changed=False,legacy_default_changed=False,extrapolation_used=False,c02_c03_c04_accepted=False,nc01_accepted=False,p6_complete=False,ab_items_completed=False,pairs=pairs,observed_orders=observed,matched_native_vs_independent=matched,path_audits=audits,legacy=legacy,separate_process_restart=st.read(out/f'g{levels[-1]}/restart-verification.json')) , events

def validate_report(report):
    for key in ('passed','physical_model_changed','force_tolerances_changed','legacy_default_changed','extrapolation_used','c02_c03_c04_accepted','nc01_accepted','p6_complete','ab_items_completed'):
        if report.get(key) is not False:
            raise ValueError('unsupported acceptance or scope change: '+key)
    if (report.get('c01_decision')!='未達' or report.get('force_tolerance')!=1e-8
            or report.get('true_linear_tolerance')!=1e-10 or report.get('independent_force_tolerance')!=1e-10):
        raise ValueError('original decision or force gates changed')
    for kind in ('native','independent'):
        limit=1e-8 if kind=='native' else 1e-10
        for n in report['global_steps']:
            if report['path_audits'][f'{kind}/{n}']['physical']['maximum_reassembled_force_residual']>limit:
                raise ValueError('original force gate exceeded')
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('input',type=Path);p.add_argument('output',type=Path);p.add_argument('--steps',type=int,nargs='+',default=st.STEPS);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False);report,events=measure(a.input.resolve(),tuple(a.steps))
    st.write(a.output/'measurements.json',validate_report(report));st.write(a.output/'switch-events.json',events)
if __name__=='__main__': main()
