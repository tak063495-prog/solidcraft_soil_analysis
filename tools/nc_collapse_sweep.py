"""Generate separate NC diagnostic inputs and run time-bounded existing paths.

Original NC inputs, validation, solver and defaults remain unchanged. Weight=0
is attempted literally; the existing ground validator rejects it. That outcome
is retained as unsupported, never replaced by a small positive weight.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
import numpy as np
from probe_nc_midpoint import inputs
from nc_fold_path import load_driver, native_loads, CONTROL, digest
from nc_collapse_classify import read, write, rows
from nc_fixed_inputs import fixed_inputs_match

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'interim/2026-11-09-nc-pseudo-arclength'))
from native_arc import Native
from low_preload_nc_reference import ClosedStickSystem
from verify_high_order import prepare


def profile(m, side, gravity):
    if m not in (1.2,1.3,1.4,1.45,1.5,1.6) or side not in (0,50,100) or type(gravity) is not bool:
        raise ValueError('requested M, side pressure and boolean gravity required')
    mesh, ground, case = inputs('midpoint')
    xyz = np.asarray(mesh['nodes']); caps=[]; sides=[]
    for f in mesh['face_sets']['exterior']:
        (caps if np.ptp(xyz[f,2]) < 1e-12 else sides).append(f)
    mesh['face_sets']['collapse_caps']=caps; mesh['face_sets']['collapse_sides']=sides
    base=ground['initialization']['geostatic']['base_loads']['surface_loads']
    original=base.pop(0)
    if original != dict(kind='pressure', face_set='exterior', pressure_kpa=100.):
        raise ValueError('unexpected original exterior loading')
    base[:0]=[dict(kind='pressure', face_set='collapse_caps', pressure_kpa=100.),
              dict(kind='pressure', face_set='collapse_sides', pressure_kpa=float(side))]
    for material in ground['materials']:
        material['mcc']['critical_state_m']=float(m)
        material['unit_weight_kn_m3']=20. if gravity else 0.
    return mesh,ground,case


def generate(output):
    output.mkdir(parents=True,exist_ok=False); records=[]
    for m in (1.2,1.3,1.4,1.45,1.5,1.6):
        for side in (0,50,100):
            for gravity in (True,False):
                name=f'M{m:g}-side{side}-weight{int(gravity)}'; folder=output/name;folder.mkdir()
                values=profile(m,side,gravity)
                for key,value in zip(('mesh','ground','stages'),values):write(folder/(key+'.json'),value)
                records.append(dict(name=name,M=m,side_kpa=side,gravity_on=gravity,
                                    input_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in folder.iterdir()}))
    write(output/'profiles.json',dict(profiles=records,generated_only=True,passed=False,nc01_complete=False))
    return records


class SweepNative(Native):
    def __init__(self,folder,binary,out):
        self.out=out;out.mkdir(parents=True,exist_ok=False);self.binary=binary;self.call=0;self.processes=[]
        self.mesh,self.ground,self.case=[read(folder/(n+'.json')) for n in ('mesh','ground','stages')]
        history=folder/'history'; accepted=rows(history/'accepted.jsonl')
        saved=read(history/'checkpoint-stage-3.json')
        if not fixed_inputs_match((self.mesh,self.ground,self.case),(self.mesh,saved['run']['ground'],saved['run']['case'])):
            raise ValueError('trace input differs')
        observed=[r['state'] for r in accepted if r['state']['stage_index']==3]
        if len(observed)<2:raise ValueError('no pair of accepted removal states')
        self.cp=observed[-1]; previous=observed[-2]; self.active=self.cp['active_elements'];self.lam0=self.cp['stage_progress']
        self.stage_origin=saved['run']['continuation']['origin']
        system=ClosedStickSystem(self.mesh,self.ground,self.active,self.case)
        if system.live.any() or np.any(system.kjoint):raise ValueError('contacts remain active')
        self.R=np.zeros((system.xyz.size,len(system.free)))
        for j,gi in enumerate(system.free):
            v=np.zeros(system.xyz.size);v[gi]=1.;self.R[:,j]=system.lift(v)
        self.q,self.offsets=prepare(self.mesh);self.rb=[];self.w=[];self.pi=[];self.par=[]
        parameters={x['id']:dict(x['mcc'],integration=self.case['settings']['mcc_integration']) for x in self.ground['materials']}
        for ei in np.flatnonzero(self.active):
            q=self.q[ei]; gd=(q['nodes'][:,None]*3+np.arange(3)).ravel()
            for qi,weight in enumerate(q['weight']):
                self.rb.append(q['b'][qi]@self.R[gd]);self.w.append(weight);self.pi.append(int(self.offsets[ei])+qi);self.par.append(parameters[self.mesh['material_ids'][ei]])
        self.rb=np.asarray(self.rb);self.w=np.asarray(self.w);self.pi=np.asarray(self.pi)
        final,self.holding=native_loads(self.mesh,self.ground,self.case,self.stage_origin,self.active)
        self.load0=final+(1-self.lam0)*self.holding;self.L=-self.R.T@self.holding
        actual=np.asarray(self.cp['applied_nodal_load_kn']).ravel()
        if np.linalg.norm(self.R.T@(actual-self.load0))>1e-10:raise ValueError('common lambda reconstruction differs')
        du=np.linalg.lstsq(self.R,(np.asarray(self.cp['displacement_m'])-previous['displacement_m']).ravel(),rcond=None)[0]
        # A bounded directional derivative selects the current loading branch.
        # The large proportional increment may use many native partitions and
        # must not be replayed as one fixed observer partition. No probe state
        # becomes the seed or a committed state.
        probe_strain=-np.einsum('qai,i->qa',self.rb,du)
        probe_scale=min(1.,1e-6/max(float(np.abs(probe_strain).sum(axis=1).max()),1e-30))
        _,self.k0,_=self.evaluate(du*probe_scale,self.cp['material_states'],'current-loading-tangent-probe')
        self.source=dict(inputs=[self.mesh,self.ground,self.case],seed=self.cp,origin=self.stage_origin,
            accepted_sha256=hashlib.sha256((history/'accepted.jsonl').read_bytes()).hexdigest())


class SweepAdapter:
    kind='native'
    def __init__(self,folder,binary,out):
        self.m=SweepNative(folder,binary,out/'material');m=self.m
        self.old=m.cp['material_states'];self.u=np.asarray(m.cp['displacement_m']);self.lam=m.lam0
        self.t=np.r_[np.linalg.solve(m.k0,m.L),1.];self.t/=np.linalg.norm(self.t)
        self.contacts=m.cp['joints'];self.R=m.R;self.L=m.L;self.source=m.source
    def evaluate(self,*a):return self.m.evaluate(*a)
    def audit(self,*a):return self.m.audit(*a)
    def serial(self,s):return s
    def state(self,s):return s


def follow(folder,binary,budget=20.,count=128,label='fold'):
    module=load_driver(); out=folder/label;out.mkdir();model=SweepAdapter(folder,binary,out);m=model.m
    final=m.load(1.);holding=m.holding
    binding=dict(fold_control=CONTROL,experiment='nc-collapse-sweep/1',model_sha256=digest(model.source),
                 source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),step_size_m=1e-4)
    write(out/'controls.json',dict(binding=binding,count=count,budget_seconds=budget,diagnostic_only=True,passed=False))
    write(out/'geometry.json',dict(map='native',mesh=m.mesh,ground=m.ground,case=m.case,active_elements=m.active,
        transform=m.R.tolist(),reduced_b_rows=m.rb.tolist(),weights_m3=m.w.tolist()))
    write(out/'original-load.json',dict(final_load_kn=final.tolist(),holding_load_kn=holding.tolist()))
    evaluate=model.evaluate; latest={};append=module.append
    def observe(du,old,label):
        result=evaluate(du,old,label)
        response=read(m.out/f'call-{m.call-1:05d}-result.json')
        latest.update(local_tangents=[p['implicit_tangent'] for p in response['points']],reduced_increment_m=np.asarray(du).tolist())
        return result
    def observe_append(path,row):
        if path.name=='samples.jsonl':row['applied_nodal_load_kn']=m.load(row['progress']).tolist()
        append(path,row)
        if path.name=='samples.jsonl':append(out/'gp-tangents.jsonl',dict(step=row['step'],**latest))
    model.evaluate=observe;module.append=observe_append
    try:report=module.follow(model,out,1e-4,count,budget,binding)
    finally:module.append=append
    write(out/'evidence-sha256.json',{str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in out.rglob('*') if p.is_file()})
    from analyze_nc_fold_path import analyze
    diagnostic,_=analyze(out);write(out/'fold-analysis.json',diagnostic)
    return dict(report=report,verdict=diagnostic['verdict'],maximum_lambda=diagnostic['peak_lambda'])


def run_one(folder,trace_binary,observer,budget=20.,fold_budget=15.):
    started=time.monotonic();command=[str(trace_binary),'trace']+[str(folder/(n+'.json')) for n in ('mesh','ground','stages')]+[str(folder/'history')]
    timed_out=False
    with (folder/'trace.log').open('x') as log:
        child=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT)
        try:child.wait(timeout=budget)
        except subprocess.TimeoutExpired:
            timed_out=True;child.terminate()
            try:child.wait(timeout=3)
            except subprocess.TimeoutExpired:child.kill();child.wait(timeout=3)
    observation=dict(command=command,exit_code=child.returncode,timed_out=timed_out,elapsed_seconds=time.monotonic()-started,
                     lambda_fold=None,full_removal_reached=False,passed=False,nc01_complete=False)
    rp=folder/'history/run.json'
    if rp.exists():
        run=read(rp);observation.update(native_complete=run['complete'],native_failure=run['failure_reason'],
                                      last_stage=run['last_committed']['stage_index'],last_progress=run['last_committed']['stage_progress'])
        observation['full_removal_reached']=bool(run['complete'] and run['last_committed']['stage_index']==3 and run['last_committed']['stage_progress']==1)
    ap=folder/'history/accepted.jsonl'
    if ap.exists():
        accepted=rows(ap);observation['accepted_rows']=len(accepted)
        if accepted:
            observation.update(last_stage=accepted[-1]['state']['stage_index'],last_progress=accepted[-1]['state']['stage_progress'])
    write(folder/'process-observation.json',observation)
    if not observation['full_removal_reached'] and observation.get('last_stage')==3 and observation.get('last_progress',0)>.85:
        try:
            arc=follow(folder,observer,fold_budget);observation['arc']=arc
            if arc['verdict'].get('conclusion')=='local_fold_bracket':observation['lambda_fold']=arc['maximum_lambda']
        except (ValueError,KeyError,FileNotFoundError,np.linalg.LinAlgError,subprocess.SubprocessError) as exc:
            observation['arc_error']=str(exc)
    write(folder/'result.json',observation)
    return observation


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--generate-only',action='store_true');p.add_argument('--trace-binary',type=Path);p.add_argument('--observer',type=Path)
    p.add_argument('--budget-seconds',type=float,default=20.);p.add_argument('--fold-budget-seconds',type=float,default=15.)
    p.add_argument('--names',nargs='*');a=p.parse_args()
    profiles=generate(a.output) if not (a.output/'profiles.json').exists() else read(a.output/'profiles.json')['profiles']
    if not a.generate_only:
        for profile_record in profiles:
            if a.names and profile_record['name'] not in a.names:continue
            folder=a.output/profile_record['name']
            if (folder/'result.json').exists():continue
            result=run_one(folder.resolve(),a.trace_binary.resolve(),a.observer.resolve(),a.budget_seconds,a.fold_budget_seconds)
            print(json.dumps(dict(name=profile_record['name'],**result)),flush=True)


def measure_bounds(folder, baseline_geometry, baseline_mode, max_iterations=80):
    """Frozen strength at the actually available field; never invent a fold pc."""
    from nc_collapse_bounds import discrete_bound, mechanism_bound, yield_values
    mesh,ground,case=[read(folder/(n+'.json')) for n in ('mesh','ground','stages')]
    process=read(folder/'result.json')
    if process.get('last_stage')!=3:
        return dict(available=False,lambda_LB=None,lambda_UB=None,reason='no accepted removal field',passed=False)
    history=folder/'history'; observed=[r['state'] for r in rows(history/'accepted.jsonl') if r['state']['stage_index']==3]
    state=observed[-1];source='last accepted proportional state; not a fold'
    for label in ('fold-384','fold'):
        if (folder/label/'samples.jsonl').exists():
            samples=rows(folder/label/'samples.jsonl')
            if samples:
                peak=max(samples,key=lambda r:r['progress'])
                state=dict(state,material_states=peak['states'],stage_progress=peak['progress'],displacement_m=peak['displacement_m'],applied_nodal_load_kn=peak['applied_nodal_load_kn'])
                source='maximum sampled local path; fold only if independently bracketed'
                break
    # All sweep meshes here retain exactly the original geometry/constraint map.
    if mesh['nodes']!=baseline_geometry['mesh']['nodes'] or mesh['tetrahedra10']!=baseline_geometry['mesh']['tetrahedra10']:
        raise ValueError('sweep bounds require matching baseline Tet10 geometry')
    active=state['active_elements']
    if active!=baseline_geometry['active_elements']:
        raise ValueError('sweep bounds topology differs')
    origin=read(history/'checkpoint-stage-3.json')['run']['continuation']['origin']
    final,h=native_loads(mesh,ground,case,origin,active)
    r=np.asarray(baseline_geometry['transform']);b=np.asarray(baseline_geometry['reduced_b_rows']);w=np.asarray(baseline_geometry['weights_m3'])
    final=r.T@final;h=r.T@h
    q,offsets=prepare(mesh);ids=np.concatenate([np.arange(offsets[i],offsets[i+1]) for i in np.flatnonzero(active)])
    pc=np.array([state['material_states'][i]['preconsolidation_pressure_kpa'] for i in ids])
    m=np.full(len(ids),ground['materials'][0]['mcc']['critical_state_m'])
    sigma=np.asarray([state['material_states'][i]['effective_stress_compression_positive_kpa'] for i in ids])
    mu=1-state['stage_progress'];residual=-np.einsum('qai,qa,q->i',b,sigma,w)-final-mu*h
    if mu==0.:
        bound=dict(status='mu_zero_witness',mu_feasible=0.,mu_dual_lower=0.,lambda_LB=1.,lambda_optimum_upper=1.,
                   equilibrium_relative=float(np.linalg.norm(residual)/max(np.linalg.norm(final),1.)),
                   maximum_yield=float(np.max(yield_values(sigma,pc,m,'ellipse'))),
                   stress_kpa=sigma.tolist(),physical_bound_certified=False)
    else:
        bound=discrete_bound(b,w,final,h,pc,m,'ellipse',max_iterations=max_iterations)
    upper=mechanism_bound(b,w,baseline_mode,final,h,pc,m)
    return dict(available=True,state_lambda=state['stage_progress'],state_source=source,fold_lambda=process['lambda_fold'],
                ellipse=bound,baseline_mode_upper=upper,passed=False,nc01_complete=False,
                limitation='Each profile freezes pc at its own available state. Unbracketed states are not fold states; only qualitative comparison is possible.')
