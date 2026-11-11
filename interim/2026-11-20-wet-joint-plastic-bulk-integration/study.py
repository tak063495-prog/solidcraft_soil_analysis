"""Bounded C-03 fresh solves and real-file restarts; preserves every attempt."""
import os
for _name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[_name]='1'
import argparse
import copy
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'tools'))
import fem_reference as fem
from reference import Model,run,STAGES,CONTACT
BASE='f377e8eda41777c84c3bb849af556bf61683b4cb'


def read(path):return json.loads(path.read_text(encoding='utf-8'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def write(path,value):
    def convert(v):
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        raise TypeError(type(v).__name__)
    with path.open('x',encoding='utf-8') as f:json.dump(value,f,default=convert,ensure_ascii=False,allow_nan=False,indent=2);f.write('\n')
def without_clock(v):
    if isinstance(v,dict):return {k:without_clock(x) for k,x in v.items() if k!='elapsed_seconds'}
    if isinstance(v,list):return [without_clock(x) for x in v]
    return v

def invoke(cmd,log,timeout):
    with log.open('x',encoding='utf-8') as stream:
        stream.write(repr(cmd)+'\n');stream.flush();start=time.monotonic()
        process=subprocess.Popen(list(map(str,cmd)),stdout=stream,stderr=subprocess.STDOUT)
        try:code=process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            process.kill();process.wait();write(log.with_suffix('.failure.json'),dict(complete=False,passed=False,pid=process.pid,reason='timeout',command=list(map(str,cmd))));raise
    return dict(command=list(map(str,cmd)),pid=process.pid,exit_code=code,elapsed_seconds=time.monotonic()-start)

def initialize(out,solver,kind,cohesion,flow=None,dilation=None,friction=None,lower_cohesion=None,open_factor=None):
    out.mkdir(parents=True,exist_ok=False)
    mesh=read(ROOT/('examples/joint-water-column-mesh.json' if kind=='tet4' else 'examples/joint-water-column-tet10-mesh.json'))
    ground=fem.normalize_ground(read(ROOT/'examples/ground-joint-water-mc.json'))
    for mat in ground['materials']:
        mat['mohr_coulomb']['cohesion_kpa']=cohesion
        if flow is not None:mat['mohr_coulomb']['flow']=flow
        if dilation is not None:mat['mohr_coulomb']['dilation_deg']=dilation
        if friction is not None:mat['mohr_coulomb']['friction_deg']=friction
        if lower_cohesion is not None and mat['id']==1:mat['mohr_coulomb']['cohesion_kpa']=lower_cohesion
    case=fem.normalize_stage_case(read(ROOT/'examples/stages-joint-water-cycle.json'))
    if open_factor is not None:
        next(s for s in case['stages'] if s['name']=='open')['target_factor']=open_factor
    write(out/'mesh.json',mesh);write(out/'ground.json',ground);write(out/'stages.json',case)
    binary=out/('probe.exe' if solver.suffix=='.exe' else 'probe');shutil.copy2(solver,binary)
    write(out/'identity.json',dict(starting_commit=BASE,solver_sha256=sha(binary),inputs={p.name:sha(p) for p in out.glob('*.json')},element=kind,cohesion_kpa=cohesion,production_changed=False))

def configuration(out,steps,material,contact):
    case=read(out/'stages.json');case['settings']['initial_increment']=1/steps
    case['settings']['minimum_increment']=min(case['settings']['minimum_increment'],1/(2*steps))
    case['settings']['material_maximum_substep']=material;case['joint_system']['maximum_increment_m']=contact
    case['verification_reject_once_stages']=[]
    for stage in case['stages']:stage['initial_increment']=1/steps
    return case

def schedule(result,steps):
    assert result['complete'] and len(result['frames'])==7
    groups=[result['geostatic']['attempts']]+[[a for a in result['attempts'] if a['stage']==name] for name in STAGES[1:]]
    report=[]
    for name,attempts in zip(STAGES,groups):
        accepted=[a for a in attempts if a['accepted']]
        assert [a['progress_to'] for a in accepted]==[i/steps for i in range(1,steps+1)],('changed global schedule',name)
        assert not [a for a in attempts if not a['accepted']],('unexpected cutback',name)
        assert max(a['relative_force_residual'] for a in accepted)<=result['case']['settings']['relative_force_tolerance']
        solves=[s for a in accepted for s in a['linear_solves']]
        assert all(s['true_relative_residual']<=result['case']['settings']['relative_linear_tolerance'] for s in solves)
        report.append(dict(stage=name,accepted=len(accepted),maximum_force=max(a['relative_force_residual'] for a in accepted),maximum_linear=max((s['true_relative_residual'] for s in solves),default=0.)))
    return report

def binary(out):
    p=out/('probe.exe' if (out/'probe.exe').exists() else 'probe');assert sha(p)==read(out/'identity.json')['solver_sha256'];return p

def native(out,folder,steps,material,contact,timeout,whole=False):
    folder.mkdir(parents=True,exist_ok=False);write(folder/'stages.json',configuration(out,steps,material,contact))
    args=[binary(out),out/'mesh.json',out/'ground.json',folder/'stages.json',folder/'raw']+(['--whole'] if whole else [])
    result=invoke(args,folder/'solve.log',timeout);write(folder/'process.json',result)
    if result['exit_code']!=0:raise RuntimeError(('native incomplete',folder,result['exit_code']))
    result=read(folder/'raw/run.json');write(folder/'schedule.json',schedule(result,steps))
    if not whole:
        verification=read(folder/'raw/verification.json');assert verification['accepted_states_written']==7*steps
        model=Model(read(out/'mesh.json'),read(out/'ground.json'),read(folder/'stages.json'))
        with (folder/'fields.jsonl').open('x') as f:
            for i,line in enumerate((folder/'raw/accepted.jsonl').open()):
                state=json.loads(line);j=state['joints']['points'];v=lambda k:np.asarray([p[k] for p in j])
                fields=model.fields(np.asarray(state['displacement_m']).ravel(),state['material_states'],v('jump_m'),v('traction_kpa'),v('plastic_slip_m'),v('accumulated_plastic_slip_m'),v('plastic_dissipation_kpa_m'),v('numerical_dissipation_kpa_m'),v('endpoint_work_kpa_m'),np.asarray(state['applied_nodal_load_kn']),[CONTACT.index(p['status']) for p in j])
                fields['physical_pressure']=np.asarray(state['joint_water']['physical_pressure_kpa']);fields['water_load']=np.asarray(state['joint_water']['nodal_load_kn']).reshape(-1,3)
                row=dict(stage=STAGES[i//steps],progress=(i%steps+1)/steps,factor=state['boundary_factor'])
                f.write(json.dumps(dict(row=row,fields={k:v.tolist() for k,v in fields.items()}),allow_nan=False)+'\n')
        assert i+1==7*steps
    print('native complete',folder.name,flush=True)

def independent(out,folder,steps,material,contact):
    folder.mkdir(parents=True,exist_ok=False);case=configuration(out,steps,material,contact);write(folder/'stages.json',case)
    model=Model(read(out/'mesh.json'),read(out/'ground.json'),case);start=time.monotonic()
    with (folder/'fields.jsonl').open('x') as f:
        def sink(row,fields):
            f.write(json.dumps(dict(row=row,fields={k:v.tolist() for k,v in fields.items()}),allow_nan=False)+'\n');f.flush()
            print(row['stage'],row['progress'],row['iterations'],flush=True)
        result=run(model,steps,material,contact,sink)
    write(folder/'result.json',result);write(folder/'process.json',dict(elapsed_seconds=time.monotonic()-start,pid=os.getpid(),reference_uses_native_results=False))

def restart(out,folder,target,timeout,steps):
    target.mkdir(parents=True,exist_ok=False);raw=folder/'raw';base=without_clock(read(raw/'run.json'));basepid=read(folder/'process.json')['pid'];rows=[]
    cps=sorted(raw.glob('checkpoint-*.json'));assert len(cps)==7
    for phase,cp in enumerate(cps):
        saved=read(cp);state=saved['run']['last_committed'];assert 0<state['stage_progress']<1
        suffix=target/cp.stem
        process=invoke([binary(out),out/'mesh.json',out/'ground.json',folder/'stages.json',suffix,'--resume',cp],target/(cp.stem+'.log'),timeout)
        assert process['exit_code']==0 and process['pid']!=basepid
        actual=read(suffix/'run.json');assert actual['complete'] and without_clock(actual)==base
        count=0;stop=int(cp.stem.split('-')[-1])
        with (raw/'accepted.jsonl').open() as a,(suffix/'accepted.jsonl').open() as b:
            for _ in range(stop):next(a)
            for x,y in itertools.zip_longest(a,b):assert x is not None and y is not None and json.loads(x)==json.loads(y);count+=1
        assert count==7*steps-stop
        rows.append(dict(stage=STAGES[phase],checkpoint=cp.name,checkpoint_sha256=sha(cp),saved_progress=state['stage_progress'],baseline_pid=basepid,resume_pid=process['pid'],suffix_states_exact=count,entire_run_exact_except_elapsed=True,process=process))
    write(target/'verification.json',dict(complete=True,real_files_separate_process=True,checks=rows,solver_sha256=sha(binary(out))))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--mode',choices=['init','native','independent','restart','whole'],required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--folder',default='baseline');p.add_argument('--source',default='baseline')
    p.add_argument('--solver',type=Path,default=ROOT/'target/release/examples/wet_joint_plastic_probe');p.add_argument('--element',choices=['tet4','tet10'],default='tet4');p.add_argument('--cohesion',type=float,default=30.)
    p.add_argument('--open-factor',type=float);p.add_argument('--lower-cohesion',type=float);p.add_argument('--flow',choices=['non_associated','associated','davis']);p.add_argument('--dilation',type=float);p.add_argument('--friction',type=float);p.add_argument('--steps',type=int,default=16);p.add_argument('--material',type=float,default=1e-3);p.add_argument('--contact',type=float,default=5e-4);p.add_argument('--timeout',type=int,default=600);a=p.parse_args();out=a.output.resolve()
    try:
        if a.mode=='init':initialize(out,a.solver.resolve(),a.element,a.cohesion,a.flow,a.dilation,a.friction,a.lower_cohesion,a.open_factor)
        elif a.mode=='independent':independent(out,out/a.folder,a.steps,a.material,a.contact)
        elif a.mode=='restart':restart(out,out/a.source,out/a.folder,a.timeout,a.steps)
        else:native(out,out/a.folder,a.steps,a.material,a.contact,a.timeout,a.mode=='whole')
    except Exception as e:
        if out.exists():write(out/f'failure-{a.mode}-{time.time_ns()}.json',dict(passed=False,complete=False,error=repr(e)))
        raise
if __name__=='__main__':main()
