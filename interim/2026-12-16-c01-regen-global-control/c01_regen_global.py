"""New input-only 3D equilibrium solves and diagnostic fixed-open controls.

The forced solve is a modified problem, excluded from C-01 acceptance.
Old native/independent results and defaults are never edited.
"""
import argparse
import copy
import importlib.util
import json
from pathlib import Path
from math import fsum, sqrt
import sys
import time
import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'interim/2026-12-13-c01-cause-diagnosis'))
import c01_cause_diagnostics as diag
sys.path.insert(0,str(ROOT/'interim/2026-12-15-c01-local-path-decomposition'))
from c01_decomposition import capture,contact_catalog,STAGES,dot,metric


def reference_module():
    # Dependency setup precedes this import; module names do not collide with
    # any of the archived study.py/reference.py/diagnose imports.
    spec=importlib.util.spec_from_file_location('c01_regeneration_equilibrium',HERE/'c01_regen_reference.py')
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def solve(out,steps,schedule=None):
    old,_=diag.dependencies();ref=reference_module()
    mesh=old.read(ROOT/'examples/joint-water-column-tet10-mesh.json')
    ground=old.fem.normalize_ground(old.read(ROOT/'examples/ground-joint-water-mc.json'))
    case=old.configuration(old.fem.normalize_stage_case(old.read(ROOT/'examples/stages-joint-water-cycle.json')),steps,5e-5)
    case['joint_system']['global_step_control']=dict(switch_bisections=4)
    case['settings']['minimum_increment']=min(case['settings']['minimum_increment'],1/(steps*32))
    out.mkdir(parents=True,exist_ok=False)
    for name,value in [('mesh',mesh),('ground',ground),('case',case)]:old.write(out/(name+'.json'),value)
    old.write(out/'diagnostic-options.json',dict(fixed_open_schedule=schedule is not None,used_for_acceptance=False))
    writer=old.Writer(out/'fields');start=time.monotonic();history=[];linear_audit=[]
    def sink(row,fields):
        writer(row,fields);history.append((copy.deepcopy(row),{k:fields[k].copy() for k in ('jump','plastic_slip','plastic_dissipation','elastic_energy','u','traction')}))
    try:
        value=ref._reference(mesh,ground,case,steps,1,7,False,True,sink,maximum_m=5e-5,switch_bisections=4,diagnostic_schedule=schedule,diagnostic_linear_audit=linear_audit)
        old.write(out/'linear-audit.json',linear_audit)
        old.write(out/'result.json',value)
        old.write(out/'summary.json',dict(complete=True,c01_passed=False,forced_open=schedule is not None,acceptance_eligible=False,steps=steps,states=len(value['path']),maximum_force_residual=value['maximum_relative_force_residual'],original_force_gate=1e-10,original_linear_gate=case['settings']['relative_linear_tolerance'],maximum_true_linear_residual=max((r['true_relative_linear_residual'] for r in linear_audit),default=0.),elapsed_seconds=time.monotonic()-start))
        return history
    except Exception as e:
        if not (out/'linear-audit.json').exists():old.write(out/'linear-audit.json',linear_audit)
        old.write(out/'failure.json',dict(complete=False,passed=False,forced_open=schedule is not None,states=len(history),error_type=type(e).__name__,error=str(e),elapsed_seconds=time.monotonic()-start))
        raise
    finally:writer.close()


def schedule_from(history,normal):
    data=dict(coords=[(r['stage'],r['progress']) for r,_ in history],fields=dict(jump=np.stack([f['jump'] for _,f in history])))
    data['index']={c:i for i,c in enumerate(data['coords'])}
    catalog=contact_catalog(data,normal)
    return dict(schema='c01-fixed-open-schedule/1',finite_fine_solution=True,intervals={s:[catalog[(s,g)] for g in range(len(normal))] for s in STAGES},source='fresh ordinary highest level; exact linear gap roots; not exact contact times',endpoint_open={s:(data['fields']['jump'][data['index'][(s,1.)]]*normal).sum(axis=1).astype(float).tolist() for s in STAGES})


def comparisons(series,levels):
    result={}
    for mode in ('ordinary','forced'):
        result[mode]={};deltas={}
        for a,b in zip(levels,levels[1:]):
            if (mode,a) not in series or (mode,b) not in series:continue
            aa,bb=series[(mode,a)],series[(mode,b)];ia={(r['stage'],r['progress']):f for r,f in aa};ib={(r['stage'],r['progress']):f for r,f in bb}
            coordinates=[(s,i/a) for s in STAGES for i in range(1,a+1)]
            fields={}
            for key in ('plastic_slip','plastic_dissipation','elastic_energy','u','jump'):
                x=np.stack([ia[c][key] for c in coordinates]);y=np.stack([ib[c][key] for c in coordinates]);fields[key]=metric(x,y);deltas[(a,key)]=(coordinates,y-x)
            result[mode][f'{a}__{b}']=fields
        a,b,c=levels
        for key in ('plastic_slip','plastic_dissipation','elastic_energy'):
            if (a,key) not in deltas or (b,key) not in deltas:continue
            coords,x=deltas[(a,key)];fine,y=deltas[(b,key)];index={v:i for i,v in enumerate(fine)};y=y[[index[v] for v in coords]]
            den=sqrt(dot(x,x))*sqrt(dot(y,y))
            result[mode][key+'_direction']=dict(dot=dot(x,y),cosine=dot(x,y)/den if den else None)
    # Intervention effect at the exact same base coordinates, without
    # interpolation or changing the original pooled norms.
    result['intervention']={}
    for n in levels:
        if ('ordinary',n) not in series or ('forced',n) not in series:continue
        a={(r['stage'],r['progress']):f for r,f in series[('ordinary',n)]};b={(r['stage'],r['progress']):f for r,f in series[('forced',n)]}
        coords=[(s,i/n) for s in STAGES for i in range(1,n+1)]
        result['intervention'][str(n)]={k:metric(np.stack([a[c][k] for c in coords]),np.stack([b[c][k] for c in coords])) for k in ('plastic_slip','plastic_dissipation','elastic_energy','u')}
    a,b,c=levels
    if all((mode,n) in series for mode in ('ordinary','forced') for n in levels):
        differences={}
        coordinates=[(s,i/a) for s in STAGES for i in range(1,a+1)]
        for mode in ('ordinary','forced'):
            maps={n:{(r['stage'],r['progress']):f for r,f in series[(mode,n)]} for n in levels}
            differences[mode]={k:[np.stack([maps[hi][p][k]-maps[lo][p][k] for p in coordinates]) for lo,hi in ((a,b),(b,c))] for k in ('plastic_slip','plastic_dissipation','elastic_energy')}
        budget={}
        for k in differences['ordinary']:
            ordinary=differences['ordinary'][k];forced=differences['forced'][k]
            effect=[ordinary[i]-forced[i] for i in range(2)]
            cs=[dict(fixed_timeline=forced[i],intervention_difference=effect[i]) for i in range(2)]
            cross={x+'/'+y:dot(v,w) for x,v in cs[0].items() for y,w in cs[1].items()}
            budget[k]=dict(original_dot=dot(*ordinary),forced_dot=dot(*forced),cross_inner_products=cross,identity_error=abs(dot(*ordinary)-fsum(cross.values())),effect_norm_over_original=[sqrt(dot(effect[i],effect[i])/dot(ordinary[i],ordinary[i])) for i in range(2)])
        result['difference_direction_intervention_budget']=budget
    return result


def intervention_artifacts(series,schedule,normal,levels):
    from c01_regen_schedule import open_at
    result={}
    for n in levels:
        if ('forced',n) not in series:continue
        disagreement=tensile=mismatch=0;max_tension=open_penetration=closed_opening=0.
        for row,f in series[('forced',n)]:
            mask=open_at(schedule,row['stage'],row['progress'])
            gap=((f['jump'][:,0]*normal[:,0]+f['jump'][:,1]*normal[:,1])+f['jump'][:,2]*normal[:,2])
            qn=((f['traction'][:,0]*normal[:,0]+f['traction'][:,1]*normal[:,1])+f['traction'][:,2]*normal[:,2])
            disagreement+=int(((mask & (gap<0)) | (~mask & (gap>0))).sum())
            tensile+=int((qn>0).sum());max_tension=max(max_tension,float(qn.max()))
            mismatch+=int(((np.asarray(row['status'])==0)!=mask).sum())
            open_penetration=max(open_penetration,float(np.max(np.where(mask,-gap,0.))))
            closed_opening=max(closed_opening,float(np.max(np.where(~mask,gap,0.))))
        result[str(n)]=dict(prescribed_open_status_mismatches=mismatch,gap_sign_disagreements=disagreement,tensile_samples=tensile,maximum_tensile_traction_kpa=max_tension,maximum_forced_open_penetration_m=open_penetration,maximum_forced_closed_opening_m=closed_opening)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--levels',nargs=3,type=int,default=[64,128,256]);a=p.parse_args()
    if a.levels not in ([16,32,64],[32,64,128],[64,128,256]):p.error('bounded nested cloud series only')
    a.output.mkdir(parents=True,exist_ok=False);old,_=diag.dependencies();series={};failures=[]
    for n in a.levels:
        print('ordinary',n,flush=True)
        try:series[('ordinary',n)]=solve(a.output/f'ordinary-{n}',n)
        except Exception as e:failures.append(dict(mode='ordinary',steps=n,error=str(e)))
    fine=series.get(('ordinary',a.levels[-1]))
    if fine:
        normal=np.array([p['normal'] for p in json.loads((ROOT/'interim/2026-12-13-c01-cause-diagnosis/geometry.json').read_text())['points']])
        schedule=schedule_from(fine,normal);old.write(a.output/'fine-open-schedule.json',schedule)
        for n in a.levels:
            print('forced',n,flush=True)
            try:series[('forced',n)]=solve(a.output/f'forced-{n}',n,copy.deepcopy(schedule))
            except Exception as e:failures.append(dict(mode='forced',steps=n,error=str(e)))
    if fine:old.write(a.output/'intervention-artifacts.json',intervention_artifacts(series,schedule,normal,a.levels))
    old.write(a.output/'comparison.json',dict(schema='c01-regen-global/1',c01_passed=False,original_problem_accepted=False,diagnostic_only=True,new_equilibrium_solves_completed=len(series),failures=failures,levels=a.levels,comparisons=comparisons(series,a.levels)))
    print('completed',len(series),'failures',len(failures),flush=True)

if __name__=='__main__':main()
