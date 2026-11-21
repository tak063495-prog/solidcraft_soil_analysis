"""Signed budgets on archived accepted jumps; replay is not equilibrium."""
import argparse
from collections import Counter
from decimal import Decimal, localcontext
import hashlib
import importlib.util
import json
from math import fsum, sqrt
from pathlib import Path
import shutil
import time
import sys
import tempfile
import zipfile
import numpy as np
from c01_path_replay import integrate, energy, dot3
from c01_local_study import write

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
PRIOR=ROOT/'interim/2026-12-13-c01-cause-diagnosis'
STAGES=('geostatic','shear','open','recontact','reverse','reload','return')
FIELDS=('plastic_slip','plastic_dissipation','elastic_energy')
COMPONENTS=('local_integration','displacement_path_same_contact','displacement_path_switch')


def read(path):return json.loads(path.read_text(encoding='utf-8'))
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def prepare(output):
    diag=module('c01_decomposition_prior',PRIOR/'c01_cause_diagnostics.py')
    for control,folder in [('switch','2026-12-11-c01-global-step-control'),('path','2026-12-12-c01-switch-event-localization')]:
        diag.extract_selected(ROOT/'interim'/folder,output/control)
    extract_prior(output)
    write(output/'provenance.json',dict(starting_commit='73e29c327745fbc045b84b23b405408a67e91fe8',source_index_sha256={name:sha(ROOT/'interim'/name/'archive-index.json') for name in ('2026-12-11-c01-global-step-control','2026-12-12-c01-switch-event-localization','2026-12-13-c01-cause-diagnosis')},prior_partition_reports_sha256={name:sha(PRIOR/name) for name in ('partition-switch.json','partition-path.json')},saved_reference_512='prior experiment, replay only; no new equilibrium solve',new_global_equilibrium_solved=False))


def extract_prior(output):
    index=read(PRIOR/'archive-index.json')
    if read(PRIOR/'acceptance.json')['archive_index_sha256'] != sha(PRIOR/'archive-index.json'):
        raise ValueError('prior index changed')
    with tempfile.TemporaryDirectory(dir=output) as tmp:
        archive=Path(tmp)/'prior.zip'
        with archive.open('xb') as out:
            for part in index['parts']:
                p=PRIOR/part['file']
                if p.stat().st_size!=part['bytes'] or sha(p)!=part['sha256']:raise ValueError('prior part integrity')
                with p.open('rb') as f:shutil.copyfileobj(f,out)
        if sha(archive)!=index['sha256']:raise ValueError('prior ZIP integrity')
        expected={x['path']:x for x in index['entries']}
        with zipfile.ZipFile(archive) as z:
            if set(z.namelist())!=set(expected) or len(z.namelist())!=len(expected):raise ValueError('prior inventory')
            for name in sorted(z.namelist()):
                if not (name.startswith('native-counts/') or name.startswith('experiments/switch-n512-q4/')):continue
                p=Path(name)
                if p.is_absolute() or '..' in p.parts:raise ValueError('unsafe archive')
                dest=output/p;dest.parent.mkdir(parents=True,exist_ok=True)
                if dest.exists():
                    if sha(dest)!=expected[name]['sha256']:raise ValueError('existing prior selection differs')
                    continue
                with z.open(name) as f,dest.open('xb') as out:shutil.copyfileobj(f,out)
                if sha(dest)!=expected[name]['sha256']:raise ValueError('prior counts integrity')



def dot(a,b):return fsum(float(x)*float(y) for x,y in zip(np.asarray(a).flat,np.asarray(b).flat))
def metric(a,b):
    d=np.asarray(a)-b;ds=dot(d,d);rs=dot(b,b)
    with localcontext() as context:
        context.prec=80
        relative=float(Decimal.from_float(ds).sqrt()/max(Decimal.from_float(rs).sqrt(),Decimal.from_float(1e-12)))
    return dict(difference_squared=ds,reference_squared=rs,relative_l2=relative)


def capture(saved,counts,old,normal,kn,kt,mu,kind,n):
    path=saved/f'g{n}'/kind if n!=512 else saved.parents[1]/'experiments/switch-n512-q4/independent'
    raw=list(old.states(path))
    coords=[(r['stage'],r['progress']) for r,_ in raw]
    if len(coords)!=len(set(coords)):raise ValueError('duplicate accepted coordinate')
    fields={name:np.stack([f[name] for _,f in raw]) for name in ('jump','traction','plastic_slip','plastic_dissipation','elastic_energy','u','physical_pressure')}
    if kind=='native':
        c=read(counts/f'n{n}/counts.json')
        if c['counts_bound_to_saved_old_and_new_displacement_bytes'] is not True or c['accepted_coordinates_exact'] is not True:raise ValueError('native counts not bound')
        rows=c['rows']
    else:rows=read(saved/f'g{n}/independent-result.json' if n!=512 else path.parent/'result.json')['path']
    mapping={(r['stage'],r['progress']):r['contact_partition_counts'] for r in rows}
    if set(mapping)!=set(coords):raise ValueError('counts coordinate mismatch')
    observed=np.array([mapping[c] for c in coords],dtype=int)
    del raw
    return dict(coords=coords,index={c:i for i,c in enumerate(coords)},fields=fields,counts=observed)


def replay(data,normal,kn,kt,mu,multiple):
    q=np.zeros_like(normal);jump=q.copy();plastic=q.copy();diss=np.zeros(len(q));out={k:[] for k in FIELDS}
    min_dd=0.;max_balance=0.;previous_e=np.zeros(len(q));endpoint_work=np.zeros(len(q));numerical=np.zeros(len(q))
    for target,count in zip(data['fields']['jump'],data['counts']):
        if multiple==0:
            use=np.maximum(1,np.ceil(np.sqrt(dot3(target-jump,target-jump))/5e-5)).astype(int)
        else:use=count*multiple
        q,dp,dd,dw,dn=integrate(jump,q,target,normal,kn,kt,mu,use)
        plastic+=dp;diss+=dd;endpoint_work+=dw;numerical+=dn
        e=energy(q,normal,kn,kt)
        min_dd=min(min_dd,float(np.min(dd)))
        max_balance=max(max_balance,float(np.max(abs(endpoint_work-e-diss-numerical))))
        out['plastic_slip'].append(plastic.copy());out['plastic_dissipation'].append(diss.copy());out['elastic_energy'].append(e)
        jump=target
    return {k:np.stack(v) for k,v in out.items()},dict(minimum_increment_dissipation=min_dd,maximum_discrete_work_balance_residual=max_balance)


def open_intervals(data,stage,lo,hi,gp,normal):
    """Exact linear gn roots on each saved segment, expressed in stage progress."""
    indices=[i for i,c in enumerate(data['coords']) if c[0]==stage and lo<c[1]<=hi]
    start_index=data['index'].get((stage,lo))
    if lo==0:
        s=STAGES.index(stage)
        start_index=None if s==0 else data['index'][(STAGES[s-1],1.)]
    g0=0. if start_index is None else float(dot3(data['fields']['jump'][start_index,gp],normal[gp]))
    t0=lo;opened=[]
    for i in indices:
        t1=data['coords'][i][1];g1=float(dot3(data['fields']['jump'][i,gp],normal[gp]))
        if g0>0 and g1>0:opened.append((t0,t1))
        elif (g0>0)!=(g1>0):
            root=t0+(t1-t0)*g0/(g0-g1)
            if g0>0:opened.append((t0,root))
            else:opened.append((root,t1))
        t0,g0=t1,g1
    # Coalesce adjacent open segments so additional saved knots do not change status history.
    merged=[]
    for a,b in opened:
        if merged and a==merged[-1][1]:merged[-1]=(merged[-1][0],b)
        else:merged.append((a,b))
    return merged


def contact_catalog(data,normal):
    return {(s,gp):open_intervals(data,s,0.,1.,gp,normal) for s in STAGES for gp in range(len(normal))}


def clipped(intervals,lo,hi):
    return [(max(a,lo),min(b,hi)) for a,b in intervals if a<hi and b>lo]


def pair_budget(a,b,n,normal):
    coords=[(s,i/n) for s in STAGES for i in range(1,n+1)]
    ia=[a['index'][c] for c in coords];ib=[b['index'][c] for c in coords]
    components={k:{c:[] for c in COMPONENTS} for k in FIELDS}
    previous_path={k:np.zeros_like(a['fields'][k][0]) for k in FIELDS}
    switch_acc={k:np.zeros_like(a['fields'][k][0]) for k in FIELDS}
    intervals=[]
    for row,(stage,t) in enumerate(coords):
        lo=(round(t*n)-1)/n
        mask=np.array([clipped(a['contacts'][(stage,gp)],lo,t)!=clipped(b['contacts'][(stage,gp)],lo,t) for gp in range(len(normal))])
        detail=dict(stage=stage,progress=t,point_count=len(normal),switch_timing_mismatch_gp=np.flatnonzero(mask).tolist(),fields={})
        for key in FIELDS:
            path=b['refined'][key][ib[row]]-a['refined'][key][ia[row]]
            local=(b['fields'][key][ib[row]]-b['refined'][key][ib[row]])-(a['fields'][key][ia[row]]-a['refined'][key][ia[row]])
            inc=path-previous_path[key]
            shaped=mask[:,None] if inc.ndim==2 else mask
            switch_acc[key]+=np.where(shaped,inc,0.)
            same=path-switch_acc[key]
            values=(local,same,switch_acc[key].copy())
            for c,v in zip(COMPONENTS,values):components[key][c].append(v)
            detail['fields'][key]=dict(path_increment_signed_sum=fsum(float(x) for x in inc.flat),switch_increment_signed_sum=fsum(float(x) for x in np.where(shaped,inc,0.).flat))
            previous_path[key]=path
        intervals.append(detail)
    summaries={}
    for key in FIELDS:
        x=a['fields'][key][ia];y=b['fields'][key][ib];delta=y-x
        cs={c:np.stack(v) for c,v in components[key].items()}
        ss=dot(delta,delta)
        summaries[key]=dict(original=metric(x,y),components={c:dict(squared_norm=dot(v,v),signed_projection_fraction=dot(v,delta)/ss if ss else None) for c,v in cs.items()},cross_inner_products={f'{c}/{d}':dot(cs[c],cs[d]) for c in COMPONENTS for d in COMPONENTS},maximum_identity_residual=float(np.max(abs(delta-sum(cs.values())))),stages={s:dict(difference_squared=dot(delta[j*n:(j+1)*n],delta[j*n:(j+1)*n]),component_projection={c:dot(v[j*n:(j+1)*n],delta[j*n:(j+1)*n]) for c,v in cs.items()}) for j,s in enumerate(STAGES)})
        components[key]=cs
    return dict(fields=summaries,displacement=metric(a['fields']['u'][ia],b['fields']['u'][ib]),jump=metric(a['fields']['jump'][ia],b['fields']['jump'][ib]),intervals=intervals),dict(coords=coords,components=components,delta={k:b['fields'][k][ib]-a['fields'][k][ia] for k in FIELDS})


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--prepare',action='store_true');p.add_argument('--refinement',type=int,default=32)
    args=p.parse_args();started=time.monotonic();args.output.mkdir(parents=True,exist_ok=False)
    if args.prepare:prepare(args.inputs)
    diag=module('c01_decomposition_prior_run',PRIOR/'c01_cause_diagnostics.py');old,_=diag.dependencies()
    saved=args.inputs/'switch/paths';mesh,ground,case=[read(saved/name) for name in ('mesh-h0.json','ground.json','stages.json')]
    system=old.operators(mesh,ground,case,'cut').system
    normal=np.array([p['normal'] for p in read(PRIOR/'geometry.json')['points']]);kn,kt,mu=system.kn,system.kt,system.mu
    report=dict(schema='c01-saved-path-decomposition/1',c01_passed=False,new_global_equilibrium_solved=False,refinement=args.refinement,coarser_refinement=args.refinement//2,components_interpretation=dict(local_integration='stored minus refined local replay; finite refinement and polygonal input only',displacement_path_same_contact='remaining saved displacement path increment, excludes switch-labelled increments; not an isolated equilibrium cause',displacement_path_switch='telescoped increment budget where open intervals differ; not an independent causal intervention',other='finite replay reference error, nonlinear coupling, unresolved causal identifiability'),audits={},pairs={},directions={})
    point_stream=(args.output/'point-increments.jsonl').open('x',encoding='utf-8',newline='\n')
    for control in ('switch','path'):
        for kind in ('native','independent'):
            series={};budgets={}
            levels=(16,32,64,128,256,512) if (control,kind)==('switch','independent') else (16,32,64,128,256)
            for n in levels:
                data=capture(args.inputs/control/'paths',args.inputs/'native-counts'/control,old,normal,kn,kt,mu,kind,n)
                matched,audit=replay(data,normal,kn,kt,mu,1)
                lower,_=replay(data,normal,kn,kt,mu,0)
                coarse,_=replay(data,normal,kn,kt,mu,args.refinement//2)
                fine,_=replay(data,normal,kn,kt,mu,args.refinement)
                data['refined']=fine
                data['contacts']=contact_catalog(data,normal)
                previous_jump=np.zeros_like(normal);previous_d=np.zeros(len(normal));previous_fine=previous_d.copy();previous_coarse=previous_d.copy()
                for j,(stage,progress) in enumerate(data['coords']):
                    target=data['fields']['jump'][j]
                    lower_count=np.maximum(1,np.ceil(np.sqrt(dot3(target-previous_jump,target-previous_jump))/5e-5)).astype(int)
                    record=dict(control=control,kind=kind,steps=n,stage=stage,progress=progress,actual_partition_counts=data['counts'][j].tolist(),endpoint_lower_bound=lower_count.tolist(),normal_pressure_old_kpa=(-kn*np.minimum(dot3(previous_jump,normal),0.)).tolist(),normal_pressure_new_kpa=(-kn*np.minimum(dot3(target,normal),0.)).tolist(),physical_water_pressure_kpa=data['fields']['physical_pressure'][j].tolist(),stored_dissipation_increment=(data['fields']['plastic_dissipation'][j]-previous_d).tolist(),refined_replay_dissipation_increment=(fine['plastic_dissipation'][j]-previous_fine).tolist(),coarser_replay_dissipation_increment=(coarse['plastic_dissipation'][j]-previous_coarse).tolist())
                    record['fields']={key:dict(stored=data['fields'][key][j].tolist(),refined_replay=fine[key][j].tolist(),coarser_replay=coarse[key][j].tolist()) for key in FIELDS}
                    point_stream.write(json.dumps(record,allow_nan=False,separators=(',',':'))+'\n')
                    previous_jump=target;previous_d=data['fields']['plastic_dissipation'][j];previous_fine=fine['plastic_dissipation'][j];previous_coarse=coarse['plastic_dissipation'][j]
                label=f'{control}/{kind}/{n}'
                report['audits'][label]=dict(states=len(data['coords']),gp_count=len(normal),partition_histogram={str(k):int(v) for k,v in sorted(Counter(data['counts'].flat).items())},matched_partition_maximum_field_errors={k:float(np.max(abs(matched[k]-data['fields'][k]))) for k in FIELDS},endpoint_lower_vs_actual_replay={k:metric(lower[k],matched[k]) for k in FIELDS},refinement_difference={k:metric(coarse[k],fine[k]) for k in FIELDS},local_defect={k:metric(data['fields'][k],fine[k]) for k in FIELDS},**audit)
                series[n]=data
                print(label,flush=True)
            for a,b in zip(levels,levels[1:]):
                summary,budget=pair_budget(series[a],series[b],a,normal)
                report['pairs'][f'{control}/{kind}/{a}__{b}']=summary;budgets[a]=budget
            for a,b in ([(32,64),(64,128),(128,256)] if len(levels)==6 else [(32,64),(64,128)]):
                first,second=budgets[a],budgets[b]
                idx={c:i for i,c in enumerate(second['coords'])};ii=[idx[c] for c in first['coords']]
                result={}
                for key in FIELDS:
                    x=first['delta'][key];y=second['delta'][key][ii]
                    with localcontext() as context:
                        context.prec=80
                        denom=Decimal.from_float(dot(x,x)).sqrt()*Decimal.from_float(dot(y,y)).sqrt()
                        cosine=float(Decimal.from_float(dot(x,y))/denom) if denom else None
                    matrix={f'{c}/{d}':dot(first['components'][key][c],second['components'][key][d][ii]) for c in COMPONENTS for d in COMPONENTS}
                    result[key]=dict(dot_difference=dot(x,y),difference_direction_cosine=cosine,component_cross_dot=matrix,maximum_dot_identity_error=abs(dot(x,y)-fsum(matrix.values())))
                report['directions'][f'{control}/{kind}/{a}_{b}_{2*b}']=result
    point_stream.close()
    report['elapsed_seconds']=time.monotonic()-started
    write(args.output/'decomposition.json',report)

if __name__=='__main__':main()
