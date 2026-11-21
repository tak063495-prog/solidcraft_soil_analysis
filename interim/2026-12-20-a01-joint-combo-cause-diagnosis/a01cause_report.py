"""Deterministic finite-discretization diagnosis; original acceptance stays false."""
import argparse
import json
import math
from pathlib import Path
import numpy as np
from a01cause_binding import fixtures,study
from a01cause_driver import write
from a01cause_metrics import load,compare,metric,norm,cosine
from a01cause_work import analyze


def aligned(data,coarse_steps):
    rows,f=data;phases=list(dict.fromkeys(r['stage'] for r in rows));chosen=[]
    for stage in phases:
        lookup={r['progress']:i for i,r in enumerate(rows) if r['stage']==stage}
        chosen.extend(lookup[i/coarse_steps] for i in range(1,coarse_steps+1))
    return [rows[i] for i in chosen],{k:v[chosen] for k,v in f.items()}


def main():
    p=argparse.ArgumentParser();p.add_argument('--matrix',type=Path,required=True);p.add_argument('--old',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    work={name:analyze(folder) for folder in sorted((a.old/'selected').iterdir()) if (folder/'comparison.json').exists() for name in [folder.name]}
    write(a.output/'work-accounting.json',dict(passed=False,diagnostic_only=True,local_replay_is_global_equilibrium=False,cases=work))
    controls={};fresh={}
    for folder in sorted(a.matrix.iterdir()):
        if not folder.is_dir() or not (folder/'result.json').exists():continue
        r=fixtures.read(folder/'result.json')
        if not r['complete']:controls[folder.name]=dict(complete=False,error=r.get('error'));continue
        case=folder.name.split('-g')[0]
        for prefix in ('adaptive-','boost-'):
            if case.startswith(prefix):case=case.removeprefix(prefix)
        if case.startswith('fixed'):case='shared-stress_free'
        steps=round(1/fixtures.read(folder/'case.json')['settings']['initial_increment'])
        baseline=(a.matrix/('native-'+case+'-g'+str(steps))/'native-fields.jsonl')
        if not baseline.exists():baseline=a.old/'selected'/case/'native-fields.jsonl'
        left=load(baseline);right=load(folder/'fields.jsonl')
        native_steps=sum(r['stage']=='geostatic' for r in left[0]);common=min(steps,native_steps)
        c=compare(aligned(left,common),aligned(right,common))
        c['native_steps']=native_steps;c['independent_steps']=steps;c['common_input_grid_steps']=common
        controls[folder.name]=dict(complete=True,uses_native_results=False,maximum_force_residual=r['maximum_force_residual'],comparison=c)
        if folder.name.startswith('adaptive-') and steps==8:
            old=load(a.old/'selected'/case/'independent-fields.jsonl');new=load(folder/'fields.jsonl')
            fresh[folder.name]=dict(old_fields_object_exact=all(np.array_equal(old[1][k],new[1][k]) for k in old[1]),comparison=compare(old,new))
    write(a.output/'global-controls.json',dict(passed=False,original_tolerances_unchanged=True,cases=controls,fresh_reference_binding=fresh))
    rates={}
    for side in ('native','adaptive'):
        for case in ('shared-stress_free', 'shared-transfer_force'):
            paths=[a.matrix/(side+'-'+case+'-g'+str(n))/('native-fields.jsonl' if side=='native' else 'fields.jsonl') for n in (4,8,16)]
            if not all(p.exists() for p in paths):continue
            data=[aligned(load(p),4) for p in paths];result={}
            for key in ('u','stress','traction','plastic_slip','endpoint_work'):
                first,second,third=(v[1][key] for v in data);d1=first-second;d2=second-third
                n1,n2=norm(d1),norm(d2)
                result[key]=dict(g4_to_g8=metric(first,second),g8_to_g16=metric(second,third),difference_vector_cosine=cosine(d1,d2),
                    observed_order=math.log2(n1/n2) if n1 and n2 else None)
            rates[side+'-'+case]=dict(common_coordinates='exact i/4 in every phase; no interpolation',fields=result)
    write(a.output/'reference-assessment.json',dict(passed=False,finest_is_exact_reference=False,extrapolation_used_for_acceptance=False,pc_executed=False,cases=rates))
    spectra=fixtures.read(a.matrix/'linear/native-spectra.json');independent=fixtures.read(a.matrix/'linear/independent-spectra.json')
    replays=fixtures.read(a.matrix/'linear/replay.json')
    last=spectra[-1];last_ind=independent[-1]
    linear=dict(passed=False,native_complete=False,independent_complete=False,numerical_noncompletion_is_scope_refusal=False,
        physical_failure_proved=False,last_native=last,last_independent=last_ind,replays=replays,
        observed_unobserved_and_old_binary_numerical_equality=fixtures.read(a.matrix/'stacked/software-observation-equality.json'))
    write(a.output/'linear-diagnosis.json',linear)
    summary=dict(passed=False,a01_complete=False,diagnostic_measurement_complete=True,new_limited_acceptance=False,
        original_acceptance_all_met=False,production_defaults_unchanged=True,refusal_scope_unchanged=True,
        pc_executed=False,cloud_diagnostics_remaining=True,old_values_rewritten=False,
        old_archive_binding=fixtures.read(a.old/'verified.json'),
        work_cases={name:dict(actual=val['actual_difference'],budgets=val['budgets'],stage_endpoints=val['stage_endpoints'],common16=val['common_subdivision_saved_path_replays']['16']) for name,val in work.items()},
        independent_controls={name:{key:c['comparison']['all_fields'][key] for key in ('u','stress','traction','endpoint_work')} if c['complete'] else c for name,c in controls.items()},
        scope_refusals=fixtures.read(a.matrix/'scope-refusals/refusals/verification.json'),
        native_last_condition_2=last['assessment']['condition_2'],native_last_direct_residual=last['assessment']['direct_lu']['true_relative_residual'],
        independent_last_condition_2=last_ind['assessment']['condition_2'],independent_last_direct_residual=last_ind['assessment']['direct_lu']['true_relative_residual'],
        untouched_statuses={k:False for k in ('a02_complete','a03_complete','a04_complete','b_items_complete','c01_complete','c02_complete','c03_complete','c04_complete','nc01_complete','p6_complete')})
    write(a.output/'measurements.json',summary)

if __name__=='__main__':main()
