"""Exhaustive field comparisons; no fitted numerical acceptance thresholds."""
import argparse
from collections import Counter
import itertools
import json
from pathlib import Path
import sys
import numpy as np
from c03diag_study import read,write,Model,STAGES,sha
from c03diag_reference import RETURNS,CONTACT
import mc_reference as mc


def fields(folder):
    with (folder/'fields.jsonl').open(encoding='utf-8') as stream:
        for line in stream:
            v=json.loads(line);yield v['row'],{k:np.asarray(x) for k,x in v['fields'].items()}

def pair(a,b,steps):
    accum={};phases={};count=0
    def add(store,key,x,y):
        xx=np.asarray(x,dtype=float);yy=np.asarray(y,dtype=float);assert xx.shape==yy.shape and np.isfinite(xx).all() and np.isfinite(yy).all(),key
        if key in ('material_return','contact_status'):
            s=store.setdefault(key,dict(mismatches=0,components=0));s['mismatches']+=int((xx!=yy).sum());s['components']+=xx.size
            return
        delta=xx-yy;s=store.setdefault(key,dict(difference_squared=0.,reference_squared=0.,maximum_absolute=0.,components=0))
        s['difference_squared']+=float(np.square(delta).sum());s['reference_squared']+=float(np.square(yy).sum());s['maximum_absolute']=max(s['maximum_absolute'],float(np.max(np.abs(delta),initial=0)));s['components']+=delta.size
    for va,vb in itertools.zip_longest(fields(a),fields(b)):
        assert va is not None and vb is not None,'missing accepted state'
        ra,fa=va;rb,fb=vb;assert ra['stage']==rb['stage']==STAGES[count//steps] and ra['progress']==rb['progress']==(count%steps+1)/steps
        assert set(fa)==set(fb)
        for k in fa:add(accum,k,fa[k],fb[k]);add(phases.setdefault(ra['stage'],{}),k,fa[k],fb[k])
        count+=1
    assert count==7*steps
    def finish(store):
        for k,v in store.items():
            if k in ('material_return','contact_status'):v['mismatch_fraction']=v['mismatches']/v['components']
            else:v['absolute_l2']=v['difference_squared']**.5;v['relative_l2']=v['absolute_l2']/max(v['reference_squared']**.5,1e-12)
    finish(accum)
    for store in phases.values():finish(store)
    return dict(accepted_pairs=count,all_fields=accum,stages=phases,numerical_acceptance=False,categories=dict(material_return=list(RETURNS),contact_status=list(CONTACT)))


def audit(out,folder,steps):
    model=Model(read(out/'mesh.json'),read(out/'ground.json'),read(folder/'stages.json'))
    n=len(model.seeds);count=0;unique=np.zeros(n,dtype=bool);rows=[];max_residual=0.;max_yield=-np.inf;returns=Counter();previous=None;water=None
    for row,f in fields(folder):
        assert row['stage']==STAGES[count//steps] and row['progress']==(count%steps+1)/steps
        assert all(np.isfinite(v).all() for v in f.values())
        eq=f['equivalent_plastic_strain'];work=f['plastic_work'];plast=np.max(np.abs(f['plastic_strain']),axis=1)>0
        yielded=(eq>0)&(work>0)&plast;unique|=yielded
        assert np.all(eq>=0) and np.all(work>=0)
        if previous is not None:
            assert np.all(eq>=previous['equivalent_plastic_strain']-1e-14) and np.all(work>=previous['plastic_work']-1e-11)
        for mat,stress in zip(model.mats,f['stress']):max_yield=max(max_yield,float(mc.normalized_yield(mat,stress)))
        np.testing.assert_allclose(f['strain'],model.strains(f['u'].ravel()),atol=1e-12,rtol=1e-10)
        np.testing.assert_allclose(f['pore'],model.pore,atol=1e-12,rtol=0)
        initial_void=np.asarray([p['void_ratio'] for p in model.seeds]);expected_void=(1+initial_void)*np.exp(f['strain'][:,:3].sum(axis=1))-1
        np.testing.assert_allclose(f['void'],expected_void,atol=1e-12,rtol=1e-10)
        phase=count//steps;factor=row['factor'];load=model.base+(1-row['progress'])*(model.seed_force-model.base) if phase==0 else model.base+factor*model.boundary
        residual=float(np.linalg.norm(model.t.T@f['reaction'].ravel())/max(np.linalg.norm(load),1));max_residual=max(max_residual,residual)
        # Existing native force gate and stricter independent gate remain distinct.
        gate=read(folder/'stages.json')['settings']['relative_force_tolerance'] if (folder/'raw').exists() else 1e-10
        assert residual<=gate,('equilibrium audit',folder.name,row,residual,gate)
        current=(f['physical_pressure'],f['water_load'])
        if water is not None:
            assert np.array_equal(water[0],current[0]) and np.array_equal(water[1],current[1])
        water=current
        np.testing.assert_allclose(current[0],model.interface['pressure'],rtol=1e-10,atol=1e-10)
        np.testing.assert_allclose(current[1],model.interface['load'],rtol=1e-9,atol=1e-9)
        rows.append(dict(stage=row['stage'],progress=row['progress'],plastic_material_gp=int(yielded.sum()),maximum_equivalent_plastic_strain=float(eq.max()),maximum_plastic_work_kpa=float(work.max()),volume_integrated_plastic_work_kn_m=float(work@model.weights),open_gp=int(f['open'].sum()),maximum_accumulated_joint_slip_m=float(f['accumulated'].max()),relative_reassembled_force_residual=residual))
        previous=f;count+=1
    assert count==7*steps and unique.any(),'plastic bulk was not exercised'
    if (folder/'raw/accepted.jsonl').exists():
        for line in (folder/'raw/accepted.jsonl').open():returns.update(p['last_return'] for p in json.loads(line)['material_states'])
    return dict(accepted_states=count,material_gp_states_checked=count*n,material_gp=n,contact_gp=len(model.system.points),unique_plastic_material_gp=int(unique.sum()),maximum_relative_reassembled_force_residual=max_residual,maximum_normalized_mc_yield=max_yield,last_returns=dict(returns),physical_pressure_and_water_load_time_invariant=True,all_accepted_points=rows)


def local_native_audit(out,folder,steps):
    """Replay local increments only, as a diagnostic distinct from fresh 3D solve."""
    from c03diag_reference import integrate
    model=Model(read(out/'mesh.json'),read(out/'ground.json'),read(folder/'stages.json'));old=model.seeds;previous_u=np.zeros(model.nd);maxima={k:0. for k in ('stress','plastic_strain','equivalent_plastic_strain','plastic_work','void')};counter=0;counts=[]
    contact_maxima={k:0. for k in ('traction','plastic_slip','accumulated','plastic_dissipation','numerical_dissipation','endpoint_work')};previous_joint=None;jcounts=[]
    from verify_joint_coupling import contact
    from wet_joint_cycle_reference import contact_work
    from types import SimpleNamespace
    keys=dict(stress='effective_stress_compression_positive_kpa',plastic_strain='plastic_strain_compression_positive',equivalent_plastic_strain='equivalent_plastic_strain',plastic_work='plastic_work_kpa',void='void_ratio')
    with (folder/'raw/accepted.jsonl').open() as states,(folder/'raw/constitutive.jsonl').open() as diagnostics:
        for line,diagnostic in itertools.zip_longest(states,diagnostics):
            assert line is not None and diagnostic is not None
            state=json.loads(line);d=json.loads(diagnostic);part=d['partitions'][:len(old)];u=np.asarray(state['displacement_m']).ravel();de=-model.strains(u-previous_u)
            np.testing.assert_allclose(de,d['strain_increment'],rtol=1e-10,atol=1e-12)
            for mat,p,e,c,actual in zip(model.mats,old,de,part,state['material_states']):
                expected=integrate(mat,p,e,c)
                for k,raw in keys.items():maxima[k]=max(maxima[k],float(np.max(np.abs(np.asarray(actual[raw])-expected[raw]))))
            joint=state['joints']['points'];system=model.system;target=np.asarray([p['jump_m'] for p in joint]);jp=d['partitions'][len(old):];jcounts.extend(jp)
            start=np.zeros_like(target) if previous_joint is None else np.asarray([p['jump_m'] for p in previous_joint])
            traction=np.zeros_like(target) if previous_joint is None else np.asarray([p['traction_kpa'] for p in previous_joint])
            for c in np.unique(jp):
                mask=np.asarray(jp)==c;sub=SimpleNamespace(**{k:getattr(system,k)[mask] for k in ('normal','kn','kt','mu')})
                q,flow,cum=contact(sub,start[mask],traction[mask],target[mask],int(c));pd,nd,work=contact_work(sub,start[mask],traction[mask],target[mask],int(c))
                for k,value,raw in [('traction',q[:,0],'traction_kpa'),('plastic_slip',flow[:,0],'plastic_slip_m'),('accumulated',cum[:,0],'accumulated_plastic_slip_m'),('plastic_dissipation',pd,'plastic_dissipation_kpa_m'),('numerical_dissipation',nd,'numerical_dissipation_kpa_m'),('endpoint_work',work,'endpoint_work_kpa_m')]:
                    actual=np.asarray([p[raw] for p in joint])[mask]
                    if k!='traction' and previous_joint is not None:actual=actual-np.asarray([p[raw] for p in previous_joint])[mask]
                    contact_maxima[k]=max(contact_maxima[k],float(np.max(np.abs(actual-value))))
            model.system.audit(state)
            old=state['material_states'];previous_joint=joint;previous_u=u;counter+=1;counts.extend(part)
    assert counter==7*steps
    return dict(diagnostic_only_not_fresh_3d_solve=True,increments=counter,point_updates=counter*len(old),maximum_absolute_local_integration_difference=maxima,material_partition_histogram=dict(sorted(Counter(counts).items())),contact_partition_histogram=dict(sorted(Counter(jcounts).items())),maximum_absolute_local_contact_difference=contact_maxima)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--steps',type=int,default=16);p.add_argument('--routes',nargs='+',default=['baseline','material-half','material-quarter','contact-half','contact-quarter']);p.add_argument('--report',type=Path,required=True);a=p.parse_args();out=a.output.resolve()
    report=dict(schema='wet-joint-plastic-bulk-integration/1',steps=a.steps,audits={},native_vs_independent={},sensitivity={},local_diagnostics={},C03=None,numerical_acceptance=False)
    for route in a.routes:
        for kind in ('native','independent'):
            name=route+'-'+kind;report['audits'][name]=audit(out,out/name,a.steps)
        report['native_vs_independent'][route]=pair(out/(route+'-native'),out/(route+'-independent'),a.steps)
        report['local_diagnostics'][route]=local_native_audit(out,out/(route+'-native'),a.steps)
    for kind in ('native','independent'):
        for family in ('material','contact'):
            for x,y in [('baseline',family+'-half'),(family+'-half',family+'-quarter')]:
                if x in a.routes and y in a.routes:report['sensitivity'][kind+'-'+x+'-to-'+y]=pair(out/(x+'-'+kind),out/(y+'-'+kind),a.steps)
    report['raw_files_sha256']={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file() and p.suffix in ('.json','.jsonl')}
    write(a.report,report)
if __name__=='__main__':main()
