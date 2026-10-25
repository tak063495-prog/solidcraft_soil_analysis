"""Independent global column solve using scalar brackets for IDH materials.

The mesh is 3D Tet4 with tied horizontal planes (oedometer kinematics). This
checks global plastic coupling, not general unconstrained 3D plastic fields.
Only native accepted load fractions define the common construction path.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import fem_reference as fem
from hardening_reference import integrate
from make_hardening_examples import initial
from verify_high_order import physical
from verify_checkpoint import invoke


def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,value):
    with p.open('x',encoding='utf-8') as stream:json.dump(value,stream,indent=2)


def reference(mesh,ground,case,native,subincrements):
    assert not mesh.get('tetrahedra10') and case['settings']['kinematics']=='oedometer' and not case.get('boundary')
    assert ground['initialization']['water_level_m'] is None and not ground['initialization'].get('geostatic')
    assert all(s['kind'] in ('pressure','excavate') for s in case['stages'])
    xyz=np.asarray(mesh['nodes']);levels=sorted(set(xyz[:,2])-{min(xyz[:,2])});t=np.zeros((xyz.size,len(levels)))
    for n,p in enumerate(xyz):
        if p[2] in levels:t[3*n+2,levels.index(p[2])]=1
    materials={m['id']:m['invariant_double_hardening'] for m in ground['materials']}
    groups=[];keys={};members=[];states=[]
    for ei in range(len(mesh['tetrahedra'])):
        op=fem.volume(mesh,ei);assert len(op['weight'])==1
        stress,pore,material,layer=fem.initial_fields(mesh,ground,ei,op['position']);assert not pore.any()
        ids=(3*op['nodes'][:,None]+np.arange(3)).reshape(-1)
        br=-op['b'][0]@t[ids];assert np.max(abs(br[[0,1,3,4,5]]))<1e-12
        a=br[2];mid=mesh['material_ids'][ei];sigma=stress[0];pc=layer['consolidation']['ratio']*sigma[:3].mean()
        key=(mesh['region_ids'][ei],mid,tuple(sigma),tuple(a))
        if key not in keys:
            keys[key]=len(groups);groups.append({'material':mid,'a':a,'weight':0.,'elements':[]})
            states.append(initial(sigma.tolist(),layer['void_ratio']['value'],float(pc)))
        group=keys[key];groups[group]['weight']+=op['weight'][0];groups[group]['elements'].append(ei);members.append(group)
    u=np.zeros(len(levels));active=np.ones(len(mesh['tetrahedra']),bool);pressure=ground['initialization']['surface_surcharge_kpa']
    applied=physical(mesh,ground,active,pressure,True)[0].reshape(-1);frames=[];max_residual=0.;zero_release=[]
    def active_groups(mask):
        result=[]
        for i,g in enumerate(groups):
            flags=mask[g['elements']];assert np.all(flags==flags[0]),'whole-body column excavation required'
            if flags[0]:result.append(i)
        return result
    def force(history,mask):
        return sum((g['weight']*g['a']*history[i]['effective_stress_compression_positive_kpa'][2] for i,g in enumerate(groups) if mask[g['elements'][0]]),np.zeros(len(levels)))
    def increment(load,mask):
        nonlocal u,states,max_residual
        live=active_groups(mask);columns=np.flatnonzero(np.any([np.abs(groups[i]['a'])>1e-12 for i in live],axis=0))
        old_u=u.copy();old=states;candidate_u=u.copy();scale=max(np.linalg.norm(load[columns]),1.)
        def evaluate(value,derivative):
            history=old.copy();internal=np.zeros(len(levels));jac=np.zeros((len(levels),len(levels)))
            for i in live:
                g=groups[i];par=materials[g['material']];strain=float(g['a']@(value-old_u));de=[0.,0.,strain,0.,0.,0.]
                state=integrate(par,old[i],de,subincrements);history[i]=state
                internal+=g['weight']*g['a']*state['effective_stress_compression_positive_kpa'][2]
                if derivative:
                    h=1e-7;plus=de.copy();minus=de.copy();plus[2]+=h;minus[2]-=h
                    a=integrate(par,old[i],plus,subincrements)['effective_stress_compression_positive_kpa'][2]
                    b=integrate(par,old[i],minus,subincrements)['effective_stress_compression_positive_kpa'][2]
                    jac+=g['weight']*(a-b)/(2*h)*np.outer(g['a'],g['a'])
            return internal,jac,history
        for iteration in range(40):
            internal,jac,history=evaluate(candidate_u,True);r=load-internal;relative=float(np.linalg.norm(r[columns])/scale)
            if relative<2e-11:
                u=candidate_u;states=history;max_residual=max(max_residual,relative);return
            correction=np.linalg.solve(jac[np.ix_(columns,columns)],r[columns]);accepted=False
            for power in range(18):
                trial=candidate_u.copy();trial[columns]+=correction/(2**power)
                try:trial_force,_,_=evaluate(trial,False)
                except (AssertionError,ValueError,OverflowError,FloatingPointError):continue
                if np.linalg.norm((load-trial_force)[columns])<np.linalg.norm(r[columns]):candidate_u=trial;accepted=True;break
            if not accepted:raise AssertionError(('independent column line search',iteration,relative))
        raise AssertionError('independent column iteration limit')
    def capture(name):
        frames.append({'name':name,'displacement_m':(t@u).reshape(-1,3).tolist(),
            'material_states':[copy.deepcopy(states[i]) for i in members],'active_elements':active.tolist()})
    assert np.linalg.norm(force(states,active)-t.T@applied)<1e-9
    capture('initial_k0')
    for stage in case['stages']:
        origin=pressure
        if stage['kind']=='excavate':
            retained=active.copy()
            for body in stage['bodies']:
                rid=next(r['id'] for r in mesh['regions'] if r['body']==body);retained[np.asarray(mesh['region_ids'])==rid]=False
            removed=active&~retained;pressure=stage['remaining_surface_pressure_kpa']
            final=physical(mesh,ground,retained,pressure,True)[0].reshape(-1)
            hold=t.T@(applied-final)-force(states,removed)
            columns=np.flatnonzero(np.any([np.abs(groups[i]['a'])>1e-12 for i in active_groups(retained)],axis=0))
            residual=force(states,retained)-(t.T@final+hold)
            zero_release.append(float(np.linalg.norm(residual[columns])/max(np.linalg.norm((t.T@final+hold)[columns]),1.)))
            assert zero_release[-1]<1e-9
            active=retained
        for attempt in native['attempts']:
            if attempt['stage']!=stage['name'] or not attempt['accepted']:continue
            fraction=attempt['progress_to']
            if stage['kind']=='pressure':
                pressure=origin+(stage['target_pressure_kpa']-origin)*fraction
                applied=physical(mesh,ground,active,pressure,True)[0].reshape(-1);load=t.T@applied
            else:load=t.T@final+(1-fraction)*hold
            increment(load,active)
        if stage['kind']=='excavate':applied=final
        capture(stage['name'])
    return {'frames':frames,'maximum_reduced_force_residual':max_residual,'zero_release_residuals':zero_release,
        'independent_subincrements_per_trial':subincrements,'exact_identical_point_groups':len(groups)}


def verify(solver,run,output):
    output.mkdir(parents=True,exist_ok=False);mesh=read(run/'ground/mesh.json');ground=read(run/'ground/ground.json');case=read(run/'stages.json');native=read(run/'stage-results.json')
    assert native['complete'] and native['provenance']['solver_binary_sha256']==sha(solver);levels=[]
    for count in (64,128,256):
        result=reference(mesh,ground,case,native,count);write(output/f'reference-{count}.json',result);levels.append(result)
        print('Independent IDH global column completed:',count,flush=True)
    refined_case=copy.deepcopy(case);refined_case['settings']['material_maximum_substep']/=4
    cp=output/'refined-stages.json';rp=output/'refined-native.json';write(cp,refined_case)
    invoke([str(solver),'solve-stages',str(run/'ground/mesh.json'),str(run/'ground/ground.json'),str(cp),str(rp)],output/'refined-native.log')
    refined=read(rp);assert refined['complete']
    # Audit the refined native solution with the same independent force and
    # raw integration-point export checks as the original fixture.
    write(output/'stages.json',refined_case)
    root=Path(__file__).resolve().parents[1]
    invoke([sys.executable,str(root/'tools/export_stages.py'),str(run/'ground/mesh.json'),str(rp),str(output/'results')],output/'export.log')
    invoke([sys.executable,str(root/'tools/verify_stage_package.py'),str(run/'ground/mesh.json'),str(run/'ground/ground.json'),str(rp),str(output/'results'),'--solver',str(solver)],output/'audit.log')
    path=lambda d:[(a['stage'],a['progress_to']) for a in d['attempts'] if a['accepted']]
    assert path(native)==path(refined),'material refinement changed accepted load fractions'
    scale=max(np.linalg.norm(f['displacement_m']) for f in levels[-1]['frames']);rows=[]
    for i,reference_frame in enumerate(levels[-1]['frames']):
        live=np.asarray(reference_frame['active_elements']);used=np.unique(fem.cells(mesh)[live]);target=np.asarray(reference_frame['displacement_m'])[used]
        errors=[float(np.linalg.norm(np.asarray(d['frames'][i]['state']['displacement_m'])[used]-target)/max(scale,1e-12)) for d in (native,refined)]
        a,b,c=[np.asarray(level['frames'][i]['displacement_m'])[used] for level in levels]
        differences=[float(np.linalg.norm(a-b)/max(scale,1e-12)),float(np.linalg.norm(b-c)/max(scale,1e-12))]
        assert errors[1]<.002,(reference_frame['name'],errors)
        assert differences[1]<.001,(reference_frame['name'],differences)
        fine_u=np.asarray(refined['frames'][i]['state']['displacement_m'])[used]
        absolute=float(np.max(np.abs(fine_u-target)))
        current=float(np.linalg.norm(fine_u-target)/max(np.linalg.norm(target),1e-12))
        increment_error=0.
        if i:
            fine_old=np.asarray(refined['frames'][i-1]['state']['displacement_m'])[used]
            reference_old=np.asarray(levels[-1]['frames'][i-1]['displacement_m'])[used]
            increment_error=float(np.linalg.norm((fine_u-fine_old)-(target-reference_old))/max(np.linalg.norm(target-reference_old),1e-12))
            assert increment_error<.005,(reference_frame['name'],increment_error)
        rows.append({'frame':reference_frame['name'],'native_coarse_fine_displacement_errors_scaled_by_maximum_path_norm':errors,'independent_refinement_differences':differences,
            'fine_maximum_active_node_absolute_displacement_error_m':absolute,'fine_current_frame_relative_displacement_error':current,'fine_stage_increment_relative_displacement_error':increment_error})
    coarse=max(r['native_coarse_fine_displacement_errors_scaled_by_maximum_path_norm'][0] for r in rows);fine=max(r['native_coarse_fine_displacement_errors_scaled_by_maximum_path_norm'][1] for r in rows)
    assert fine<coarse*.8,(coarse,fine)
    report={'schema':'soil-fem-independent-hardening-column/1','passed':True,'solver_sha256':sha(solver),'frames':rows,
        'reference_subincrements':[64,128,256],'native_maximum_substeps':[case['settings']['material_maximum_substep'],refined_case['settings']['material_maximum_substep']],
        'maximum_reference_force_residual':max(v['maximum_reduced_force_residual'] for v in levels),
        'scope':'Input-only 3D Tet4 mesh with oedometer plane constraints; independent quadratic shear elimination/scalar cap bracket, numerical material derivative and NumPy global Newton. Includes plastic preload, unloading, reloading and upper-body excavation with independently computed holding forces. It does not certify general unconstrained 3D plastic displacement.',
        'sha256':{str(p):sha(p) for p in [run/'ground/mesh.json',run/'ground/ground.json',run/'stages.json',run/'stage-results.json',Path(__file__).resolve(),Path(__file__).with_name('hardening_reference.py')]}}
    write(output/'verification.json',report);print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--run',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();verify(a.solver.resolve(),a.run.resolve(),a.output_dir.resolve())
