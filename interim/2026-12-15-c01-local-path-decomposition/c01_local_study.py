"""Analytical local paths compared with the production Rust point law."""
import argparse
import hashlib
import json
from math import fsum, log2, sqrt
from pathlib import Path
import subprocess
from c01_local_oracle import advance_exact, initial, refine
from c01_path_replay import integrate
import numpy as np

KN, KT, MU = 1e6, 2e5, .4
PATHS = {
    'rising_pressure': [[0,-.0001],[.003,-.0005]],
    'reverse_variable_pressure': [[0,-.0005],[.0015,-.0001],[-.0018,-.00035],[.0025,-.00015]],
    'opening_recontact_reverse': [[0,-.0002],[.002,-.0002],[.0035,.00015],[-.002,.0001],[-.003,-.0002],[.002,-.0003]],
    'pressure_release_at_fixed_slip': [[0,-.0005],[.002,-.0005],[.002,-.0001],[-.001,-.0003]],
}


def write(path,value):
    with path.open('x',encoding='utf-8',newline='\n') as f:
        json.dump(value,f,indent=2,ensure_ascii=False,allow_nan=False); f.write('\n')


def make_request():
    runs=[]; metadata={}
    for name,vertices in PATHS.items():
        for n in (1,2,4,8,16,32):
            points=refine(vertices,n)
            for policy,count,overshoot in [('fixed',c,False) for c in (1,2,4,8)]+[('endpoint_lower',1,False),('frozen_rule_replay',1,False),('frozen_rule_replay',1,True)]:
                key=f'{name}/n{n}/{policy}-{count}'+('-overshoot' if overshoot else '')
                steps=[]
                for _,_,(x,g) in points:
                    trial=[[x+.0006,0.,g]] if overshoot else []
                    steps.append(dict(jump_m=[x,0.,g],trial_jumps_m=trial))
                runs.append(dict(id=key,policy=policy,fixed_count=count,steps=steps))
                metadata[key]=(name,n,policy,count,overshoot,points)
    return dict(material=dict(normal_stiffness_kpa_m=KN,tangential_stiffness_kpa_m=KT,friction_coefficient=MU),normal=[0.,0.,1.],maximum_increment_m=5e-5,runs=runs),metadata


def compare(response):
    _,metadata=make_request(); result={}
    for run in response:
        name,n,policy,count,overshoot,points=metadata[run['id']]
        exact=initial(); before=initial(); q=np.zeros((1,3)); p=np.zeros(3); d=0.
        previous_jump=np.zeros((1,3)); previous_native_d=0.
        rows=[]; increments=[]; defects=[]; slips=[]; tractions=[]; endpoint_errors=[]
        min_dd=0.; min_num=0.; max_balance=0.; partitions=[]
        for (segment,i,target),sample in zip(points,run['samples'],strict=True):
            exact=advance_exact(exact,target,KN,KT,MU)
            native=sample['state']; counts=sample['actual_count']; partitions.append(counts)
            dest=np.array([[target[0],0.,target[1]]])
            qq,dp,dd,_,_=integrate(previous_jump,q,dest,np.array([[0.,0.,1.]]),np.array([KN]),np.array([KT]),np.array([MU]),counts)
            p+=dp[0]; d+=float(dd[0]); q=qq; previous_jump=dest
            slips.append(abs(native['plastic_slip_m'][0]-exact['slip']))
            tractions.append(abs(native['traction_kpa'][0]-exact['tau']))
            defects.append(native['plastic_dissipation_kpa_m']-exact['dissipation'])
            increments.append(sample['plastic_dissipation_increment_kpa_m']-(exact['dissipation']-before['dissipation']))
            endpoint_errors.append(max(abs(d-native['plastic_dissipation_kpa_m']),*(abs(float(a)-b) for a,b in zip(p,native['plastic_slip_m'])),*(abs(float(a)-b) for a,b in zip(q[0],native['traction_kpa']))))
            min_dd=min(min_dd,sample['plastic_dissipation_increment_kpa_m'])
            min_num=min(min_num,native['numerical_dissipation_kpa_m'])
            balance=sample['elastic_energy_kpa_m']+native['plastic_dissipation_kpa_m']+native['numerical_dissipation_kpa_m']-native['endpoint_work_kpa_m']
            max_balance=max(max_balance,abs(balance))
            rows.append(dict(segment=segment,increment=i,exact=exact,native_dissipation=native['plastic_dissipation_kpa_m'],native_dissipation_increment=sample['plastic_dissipation_increment_kpa_m'],exact_dissipation_increment=exact['dissipation']-before['dissipation'],count=counts,lower=sample['endpoint_lower_bound']))
            before=exact;previous_native_d=native['plastic_dissipation_kpa_m']
        result[run['id']]=dict(path=name,outer_subdivisions_per_segment=n,policy=policy,fixed_count=count,prescribed_rejected_trial=overshoot,maximum_slip_error_m=max(slips),maximum_traction_error_kpa=max(tractions),maximum_increment_dissipation_error_kpa_m=max(map(abs,increments)),maximum_cumulative_dissipation_error_kpa_m=max(map(abs,defects)),final_dissipation_error_kpa_m=defects[-1],exact_final_dissipation_kpa_m=exact['dissipation'],native_final_dissipation_kpa_m=previous_native_d,minimum_increment_dissipation_kpa_m=min_dd,minimum_cumulative_numerical_dissipation_kpa_m=min_num,maximum_discrete_balance_residual_kpa_m=max_balance,maximum_independent_endpoint_difference=max(endpoint_errors),partition_counts=partitions,history=rows)
    order={}
    for name in PATHS:
        for policy in ('fixed-1','fixed-8','endpoint_lower-1','frozen_rule_replay-1','frozen_rule_replay-1-overshoot'):
            keys=[f'{name}/n{n}/{policy}' for n in (1,2,4,8,16,32)]
            if any(k not in result for k in keys):continue
            series=[result[k]['final_dissipation_error_kpa_m'] for k in keys]
            order[f'{name}/{policy}']=dict(signed_errors=series,absolute_errors_monotone=all(abs(b)<=abs(a) for a,b in zip(series,series[1:])),observed_orders=[log2(abs(a/b)) if a*b !=0 else None for a,b in zip(series,series[1:])])
    return dict(schema='c01-local-analytical-paths/1',c01_passed=False,new_global_equilibrium_solved=False,analytic_solution='scalar piecewise affine sweep; exact yield/open roots and integral of linear friction cap',runs=result,orders=order)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    request,_=make_request();write(a.output/'request.json',request)
    subprocess.run([str(a.solver.resolve()),str((a.output/'request.json').resolve()),str((a.output/'response.json').resolve())],check=True)
    summary=compare(json.loads((a.output/'response.json').read_text(encoding='utf-8')))
    write(a.output/'local-results.json',summary)
    print(json.dumps(dict(runs=len(summary['runs']),maximum_independent_endpoint_difference=max(v['maximum_independent_endpoint_difference'] for v in summary['runs'].values()))))

if __name__=='__main__':main()
