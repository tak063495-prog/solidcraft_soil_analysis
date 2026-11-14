"""Read-only C-04 full-history measurement and local replay, separate from solves."""
import argparse
from collections import Counter
import json
import math
from pathlib import Path
import sys
from types import SimpleNamespace
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE));sys.path.insert(0,str(ROOT/'tools'))
import c04cause_metrics as metrics
from c04cause_driver import reference,write
from verify_joint_coupling import contact
from wet_joint_cycle_reference import contact_work


def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))


def subset(system,mask):
    return SimpleNamespace(**{k:getattr(system,k)[mask] for k in ('normal','kn','kt','mu')})


def partitions(folder,n):
    if (folder/'raw/constitutive.jsonl').exists():
        return np.asarray([json.loads(x)['partitions'][n:] for x in (folder/'raw/constitutive.jsonl').read_text().splitlines()])
    return np.asarray([r['contact_partitions'] for r in read(folder/'result.json')['rows']])


def last_update(system,start,oldq,target,count):
    before=start.copy();before_q=oldq.copy()
    if count>1:
        # Same equal subsegment positions as contact(); no n-1 repartition of the full target.
        before=start+(target-start)*(count-1)/count
        before_q=contact(system,start,oldq,before,count-1)[0][:,0]
    n=system.normal;dt=target-before;dt-=np.sum(dt*n,axis=1)[:,None]*n
    tau0=before_q-np.sum(before_q*n,axis=1)[:,None]*n
    trial=tau0+system.kt[:,None]*dt
    pressure=-system.kn*np.minimum(np.sum(target*n,axis=1),0.)
    radius=system.mu*pressure;trialnorm=np.linalg.norm(trial,axis=1)
    flow=contact(system,before,before_q,target,1)[1][:,0]
    return dict(trial_margin_kpa=radius-trialnorm,trial_margin_scaled=(radius-trialnorm)/(1+radius+trialnorm),
                last_flow_m=np.linalg.norm(flow,axis=1),radius_kpa=radius)


def margin_bin(value):
    value=abs(float(value))
    if value<=1e-14:return '0..1e-14'
    if value<=1e-12:return '1e-14..1e-12'
    if value<=1e-10:return '1e-12..1e-10'
    if value<=1e-8:return '1e-10..1e-8'
    return '>1e-8'


def mismatch_rows(a,b,system,ca,cb):
    rows,fa=a;_,fb=b;result=[]
    oldja=np.zeros_like(fa['jump'][0]);oldjb=oldja.copy();oldqa=oldja.copy();oldqb=oldja.copy()
    for i,row in enumerate(rows):
        mismatch=np.flatnonzero(fa['contact_status'][i]!=fb['contact_status'][i])
        for gp in mismatch:
            sa=subset(system,np.array([gp]));targeta=fa['jump'][i,gp:gp+1];targetb=fb['jump'][i,gp:gp+1]
            aa=last_update(sa,oldja[gp:gp+1],oldqa[gp:gp+1],targeta,int(ca[i,gp]))
            bb=last_update(sa,oldjb[gp:gp+1],oldqb[gp:gp+1],targetb,int(cb[i,gp]))
            entry=dict(index=i,stage=row['stage'],progress=row['progress'],gp=int(gp),
                       native_status=int(fa['contact_status'][i,gp]),independent_status=int(fb['contact_status'][i,gp]),
                       native_parts=int(ca[i,gp]),independent_parts=int(cb[i,gp]),
                       jump_difference_m=metrics.norm(targeta-targetb),traction_difference_kpa=metrics.norm(fa['traction'][i,gp]-fb['traction'][i,gp]))
            for label,data,fields,old in (('native',aa,fa,oldqa),('independent',bb,fb,oldqb)):
                q=fields['traction'][i,gp];n=system.normal[gp];tau=q-n*np.dot(q,n)
                radius=float(data['radius_kpa'][0]);margin=radius-metrics.norm(tau)
                entry[label]=dict(**{k:float(v[0]) for k,v in data.items()},returned_margin_kpa=margin,
                                  returned_margin_scaled=margin/(1+radius+metrics.norm(tau)))
            entry['margin_bin']=margin_bin(max(abs(entry[k]['trial_margin_scaled']) for k in ('native','independent')))
            entry['zero_saved_plastic_increment_both']=all(np.array_equal(f['plastic_slip'][i,gp],f['plastic_slip'][i-1,gp] if i else np.zeros(3)) for f in (fa,fb))
            result.append(entry)
        oldja=fa['jump'][i];oldjb=fb['jump'][i];oldqa=fa['traction'][i];oldqb=fb['traction'][i]
    return result


def dissipation_identity(a,b,system):
    """Exact algebra of endpoint products, with a separately measured within-step remainder."""
    fa=a[1];fb=b[1]
    pa=np.diff(fa['plastic_slip'],axis=0,prepend=np.zeros_like(fa['plastic_slip'][:1]))
    pb=np.diff(fb['plastic_slip'],axis=0,prepend=np.zeros_like(fb['plastic_slip'][:1]))
    n=system.normal
    qa=fa['traction']-np.sum(fa['traction']*n,axis=2)[:,:,None]*n
    qb=fb['traction']-np.sum(fb['traction']*n,axis=2)[:,:,None]*n
    traction=np.sum((qa-qb)*(pa+pb)*.5,axis=2)
    slip=np.sum((qa+qb)*(pa-pb)*.5,axis=2)
    # Fixed temporal order: fsum of every prefix, then GP order in norm.
    prefix=lambda x:np.asarray([[math.fsum(map(float,x[:i+1,g])) for g in range(x.shape[1])] for i in range(len(x))])
    traction=prefix(traction);slip=prefix(slip)
    actual=fa['plastic_dissipation']-fb['plastic_dissipation']
    remainder=actual-traction-slip
    return dict(actual=actual,endpoint_traction_contribution=traction,endpoint_plastic_increment_contribution=slip,
                within_increment_remainder=remainder,identity_residual_max=float(np.max(abs(actual-traction-slip-remainder))),
                explanation='The remainder is not an equilibrium residual or a bug verdict. Endpoint products omit changing traction and plastic slip within local substeps; vector increments can cancel.')


def replay(jumps,system,count):
    nq=len(system.normal);start=np.zeros((nq,3));q=start.copy();pd=np.zeros(nq);plastic=start.copy();out=[];qs=[];ps=[]
    for target in jumps:
        end,dp,_=contact(system,start,q,target,count)
        pd=pd+contact_work(system,start,q,target,count)[0]
        plastic=plastic+dp[:,0];q=end[:,0];start=target
        out.append(pd.copy());qs.append(q.copy());ps.append(plastic.copy())
    return dict(plastic_dissipation=np.asarray(out),traction=np.asarray(qs),plastic_slip=np.asarray(ps))


def common_pair(a,b):
    ia={(r['stage'],r['progress']):i for i,r in enumerate(a[0])}
    ib={(r['stage'],r['progress']):i for i,r in enumerate(b[0])}
    keys=[k for k in ia if k in ib]
    if not keys:raise ValueError('no exact common history coordinates')
    return (([a[0][ia[k]] for k in keys],{n:v[[ia[k] for k in keys]] for n,v in a[1].items()}),
            ([b[0][ib[k]] for k in keys],{n:v[[ib[k] for k in keys]] for n,v in b[1].items()}))


def audit_global(folder,root):
    model=reference.Model(read(root/'mesh.json'),read(root/'ground.json'),read(folder/'case.json'))
    rows,fields=metrics.load(folder/'fields.jsonl');case=model.case
    phases=[dict(name='geostatic',n=round(1/case['settings']['initial_increment']))]+[dict(name=x['name'],n=round(1/x['initial_increment'])) for x in case['stages']]
    expected=[(p['name'],i/p['n']) for p in phases for i in range(1,p['n']+1)]
    if [(r['stage'],r['progress']) for r in rows]!=expected:raise ValueError('input accepted schedule differs')
    gate=case['settings']['relative_force_tolerance'] if (folder/'raw').exists() else 1e-10
    maximum=0.;water=0.;reaction=0.
    for i,row in enumerate(rows):
        model.set_level(row['level_m'])
        states=[dict(effective_stress_compression_positive_kpa=x) for x in fields['stress'][i]]
        residual=model.internal(states,fields['traction'][i])-fields['applied_load'][i].ravel()
        relative=metrics.norm(model.t.T@residual)/max(metrics.norm(fields['applied_load'][i]),1.)
        if relative>gate:raise ValueError(('original equilibrium gate',folder.name,i,relative,gate))
        maximum=max(maximum,relative)
        water=max(water,float(np.max(abs(fields['physical_pressure'][i]-model.interface['pressure']))))
        reaction=max(reaction,float(np.max(abs(residual-fields['reaction'][i].ravel()))))
    mapping=None
    if (folder/'raw/run.json').exists():
        run=read(folder/'raw/run.json');native=run['joint_integration_points'];system=model.system
        from joint_reference import System
        System(model.mesh,case['joint_system'],native)  # Existing correspondence checks, no nearest-point matching.
        material=run['integration_points'];coordinates=[]
        for ei,(q,ids) in enumerate(model.ops):
            coordinates.extend((ei,li) for li in range(len(q['weight'])))
        if coordinates!=[(p['element_index'],p['local_index']) for p in material]:raise ValueError('material GP order differs')
        mapping=dict(material_gp=len(material),contact_gp=len(native),material_ids_order_exact=True,contact_discrete_ids_order_exact=True,
          maximum_material_position_difference_m=float(np.max(abs(np.asarray(model.positions)-np.asarray([p['position_m'] for p in material])))),
          maximum_material_weight_difference_m3=float(np.max(abs(np.asarray(model.weights)-np.asarray([p['weight_m3'] for p in material])))),
          maximum_contact_position_difference_m=float(np.max(abs(np.asarray([p['position_m'] for p in system.points])-np.asarray([p['position_m'] for p in native])))),
          maximum_contact_weight_difference_m2=float(np.max(abs(system.weight-np.asarray([p['weight_m2'] for p in native])))),
          independent_input_descriptors=system.points,native_contact_descriptors=native,native_material_descriptors=material)
    return dict(accepted_states=len(rows),input_schedule_exact=True,maximum_reassembled_force_residual=maximum,original_force_gate=gate,
                maximum_pressure_transfer_difference_kpa=water,maximum_stored_reaction_difference_kn=reaction,gp_correspondence=mapping)


def observe(baseline,output,controls=None):
    output.mkdir(parents=True,exist_ok=False)
    mesh=read(baseline/'mesh.json');ground=read(baseline/'ground.json')
    report=dict(schema='c04-cause-observation/1',passed=False,c04_complete=False,
                local_replay_is_global_equilibrium=False,original_l2_definition_retained=True,
                reduction_order='accepted time, GP/node, component C order; math.fsum',levels={},sensitivity={})
    loaded={}
    for water in (8,16,32):
        na=baseline/f'native-w{water}';ib=baseline/f'independent-w{water}'
        a=metrics.load(na/'fields.jsonl');b=metrics.load(ib/'fields.jsonl');loaded['native',water]=a;loaded['independent',water]=b
        model=reference.Model(mesh,ground,read(na/'case.json'));system=model.system
        ca=partitions(na,len(model.seeds));cb=partitions(ib,len(model.seeds))
        if ca.shape!=cb.shape or ca.shape!=(len(a[0]),len(system.normal)):raise ValueError('partition correspondence')
        mismatch=mismatch_rows(a,b,system,ca,cb)
        write(output/f'mismatch-w{water}.json',mismatch)
        identity=dissipation_identity(a,b,system)
        write(output/f'dissipation-identity-w{water}.json',identity)
        actual=identity['actual'];summary={}
        for k in ('endpoint_traction_contribution','endpoint_plastic_increment_contribution','within_increment_remainder'):
            summary[k]=metrics.budget(actual,identity[k])
        increments=np.diff(actual,axis=0,prepend=np.zeros_like(actual[:1]))
        stage_increment_norms={name:metrics.norm(increments[[i for i,r in enumerate(a[0]) if r['stage']==name]]) for name in dict.fromkeys(r['stage'] for r in a[0])}
        common={}
        for count in (1,2,4,8,16):
            ar=replay(a[1]['jump'],system,count);br=replay(b[1]['jump'],system,count)
            write(output/f'local-replay-w{water}-c{count}.json',dict(native=ar,independent=br,global_equilibrium_solved=False))
            delta=ar['plastic_dissipation']-br['plastic_dissipation']
            common[str(count)]=dict(difference=metrics.budget(actual,delta),
              native_stored_vs_replay=metrics.metric(a[1]['plastic_dissipation'],ar['plastic_dissipation']),
              independent_stored_vs_replay=metrics.metric(b[1]['plastic_dissipation'],br['plastic_dissipation']))
        report['levels'][str(water)]=dict(comparison=metrics.compare(a,b),mismatch_count=len(mismatch),
          margin_bins=dict(sorted(Counter(r['margin_bin'] for r in mismatch).items())),
          zero_saved_plastic_increment_both=sum(r['zero_saved_plastic_increment_both'] for r in mismatch),
          native_parts=dict(sorted(Counter(ca.ravel().tolist()).items())),independent_parts=dict(sorted(Counter(cb.ravel().tolist()).items())),
          dissipation_identity=summary,stage_dissipation_increment_difference_l2=stage_increment_norms,common_saved_jump_replay=common,
          native_global_audit=audit_global(na,baseline),independent_global_audit=audit_global(ib,baseline))
    for kind in ('native','independent'):
        for lo,hi in ((8,16),(16,32)):
            report['sensitivity'][f'{kind}-w{lo}-w{hi}']=metrics.compare(*common_pair(loaded[kind,lo],loaded[kind,hi]))
    if controls:
        report['global_controls']={}
        for folder in sorted(controls.iterdir()):
            if not (folder/'fields.jsonl').exists():continue
            data=metrics.load(folder/'fields.jsonl')
            report['global_controls'][folder.name]=dict(comparison=metrics.compare(*common_pair(loaded['native',32],data)),global_audit=audit_global(folder,controls))
    write(output/'measurements.json',report)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--baseline',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--controls',type=Path)
    a=p.parse_args();observe(a.baseline.resolve(),a.output.resolve(),a.controls.resolve() if a.controls else None)

if __name__=='__main__':main()
