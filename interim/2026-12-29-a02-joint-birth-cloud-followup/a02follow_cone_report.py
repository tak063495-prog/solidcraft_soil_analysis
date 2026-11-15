"""Reproducible summary of new ULP/build/rounding evidence."""
import argparse
from pathlib import Path
import numpy as np
from a02follow_binding import study,metric
from a02follow_run import write

def report(cone,rounding,out):
    v=study.read(cone/'cone-results.json');rows=v['rows'];rr=study.read(rounding/'rounding-results.json')['rows'];local=[r for r in rows if r['kind']=='local'];projection=[r for r in rows if r['kind']=='projection']
    summary=dict(passed=False,a02_complete=False,diagnostic_only=True,local_requests=len(local),projection_requests=len(projection),rotated_requests=len(rr),build_profiles=['release','debug'],
        unchanged_formula_tolerance_clip=True,windows_reproduction_confirmed=False,windows_ci_after_push='未確定',
        native_build_mismatches=sum(not r['builds_exact'] for r in rows+rr),
        rejected_work_is_null=all(a.get('work') is None for r in rows+rr for a in r['answers'].values() if not a['accepted']),
        local={},projection={},rounding={})
    for name in sorted({r['case'] for r in local}):
        x=[r for r in local if r['case']==name];works=[r['answers']['release']['state']['endpoint_work_kpa_m'] for r in x if r['answers']['release'].get('state')]
        summary['local'][name]=dict(accepted_ulp=[r['ulp_steps'] for r in x if r['answers']['release']['accepted']],rejected_ulp=[r['ulp_steps'] for r in x if not r['answers']['release']['accepted']],
            work_minimum=min(works),work_maximum=max(works),work_span=max(works)-min(works),
            statuses=sorted({r['answers']['release']['state']['status'] for r in x if r['answers']['release'].get('state')}))
    for name in sorted({r['case'] for r in projection}):
        x=[r for r in projection if r['case']==name];raw=lambda r:np.array([q for _,q in sorted(r['answers']['release']['raw_projection'].items(),key=lambda z:int(z[0]))])
        def change(r):
            base=next(t for t in x if t['perturbation']==r['perturbation'] and t['ulp_steps']==0)
            return metric(raw(r),raw(base))['maximum_absolute']
        native_excess=[]
        for r in x:
            q=raw(r);pressure=-q[:,2];native_excess.append(float(max(np.linalg.norm(q[:,:2],axis=1)-.3*pressure)))
        summary['projection'][name]=dict(requests=len(x),native_accepted=sum(r['answers']['release']['accepted'] for r in x),
            errors=sorted({r['answers']['release'].get('error','accepted') for r in x}),
            maximum_input_change_kn=max(r['input_change']['maximum_absolute'] for r in x),maximum_native_projection_change_kpa=max(change(r) for r in x),
            native_maximum_friction_excess_minimum_kpa=min(native_excess),native_maximum_friction_excess_maximum_kpa=max(native_excess),
            maximum_gelsd_gelss_traction_change_kpa=max(r['linear_path_difference']['maximum_absolute'] for r in x),
            gelsd_gelss_cone_disagreements=sum(r['independent']['gelsd']['cone']!=r['independent']['gelss']['cone'] for r in x))
    changed=[dict(normal=r['normal'],axis=r['axis'],ulp_steps=r['ulp_steps'],orders=r['orders'],native_accepted=r['answers']['release']['accepted']) for r in rr if r['orders']['order_changes_gate']]
    summary['rounding']=dict(requests=len(rr),native_accepted=sum(r['answers']['release']['accepted'] for r in rr),same_input_order_gate_changes=len(changed),changed=changed,
        maximum_excess_order_span_kpa=max(max(o['excess'] for o in r['orders']['permutations'])-min(o['excess'] for o in r['orders']['permutations']) for r in rr),
        explicit_sequential_addition=True,python_builtin_sum_not_used_for_order_probe=True,diagnostic_order_variants_are_not_adopted=True)
    write(out,summary);return summary
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cone',type=Path,required=True);p.add_argument('--rounding',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();report(a.cone,a.rounding,a.output)
