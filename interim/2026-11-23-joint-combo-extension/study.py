"""Fresh solves, actual checkpoint/PID recovery, input refusals and comparison."""
import argparse
import copy
import hashlib
import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import traceback
import zipfile
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
import numpy as np
from fixtures import ROOT,read,write,cases,separated,structure,specimen
from reference import Model,run
from fields import Sink,native_fields,compare,load_fields,CATEGORIES

HERE=Path(__file__).resolve().parent
BASE='c70fb37ea025b6e40dda4db4adf67379ea2a6ee3'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def exact(value):
    def clean(v):
        if isinstance(v,dict):return {k:clean(x) for k,x in v.items() if k!='elapsed_seconds'}
        if isinstance(v,list):return [clean(x) for x in v]
        return v
    return json.dumps(clean(value),sort_keys=True,allow_nan=False,separators=(',',':'))


def invoke(command,log,timeout=900):
    started=time.monotonic()
    with Path(log).open('x',encoding='utf-8') as stream:
        stream.write('Command: '+repr(list(map(str,command)))+'\n');stream.flush()
        p=subprocess.Popen(list(map(str,command)),stdout=stream,stderr=subprocess.STDOUT)
        timed_out=False
        try:code=p.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out=True;p.kill();code=p.wait()
        stream.write(f'\nexit_code={code}; timed_out={timed_out}; pid={p.pid}\n')
    return dict(command=list(map(str,command)),pid=p.pid,exit_code=code,timed_out=timed_out,elapsed_seconds=time.monotonic()-started)


def inputs(folder):return [folder/(key+'.json') for key in ('mesh','ground','case')]


def native(folder,probe,steps,timeout):
    result=invoke([probe,*inputs(folder),folder/'native'],folder/'native.log',timeout)
    write(folder/'native-process.json',result)
    if result['exit_code']!=0 or result['timed_out']:raise RuntimeError(('native incomplete',folder.name,result))
    raw=folder/'native';native=read(raw/'run.json');assert native['complete']
    phases=['geostatic']+[s['name'] for s in read(folder/'case.json')['stages']]
    groups=[native['geostatic']['attempts']]+[[a for a in native['attempts'] if a['stage']==name] for name in phases[1:]]
    for name,attempts in zip(phases,groups):
        assert all(a['accepted'] for a in attempts),(name,attempts)
        assert [a['progress_to'] for a in attempts]==[i/steps for i in range(1,steps+1)],('input grid differs',name)
        assert all(a['relative_force_residual']<=native['case']['settings']['relative_force_tolerance'] for a in attempts)
        assert all(s['true_relative_residual']<=native['case']['settings']['relative_linear_tolerance'] for a in attempts for s in a['linear_solves'])
    model=Model(*(read(p) for p in inputs(folder)))
    with (raw/'accepted.jsonl').open() as states,(raw/'constitutive.jsonl').open() as diagnostics,(raw/'structure-reports.jsonl').open() as reports,(folder/'native-fields.jsonl').open('x') as stream:
        sink=Sink(stream);count=0
        for a,b,c in itertools.zip_longest(states,diagnostics,reports):
            assert all(x is not None for x in (a,b,c)),'missing native observation'
            state=json.loads(a);d=json.loads(b);report=json.loads(c)
            phase=count//steps;progress=(count%steps+1)/steps
            assert state['stage_index']==phase
            # Equilibration completion resets the committed continuation
            # cursor to 0; the accepted diagnostic/event is still progress 1.
            expected_cursor=0. if phase==0 and progress==1. else progress
            assert state['stage_progress']==expected_cursor,(phase,progress,state['stage_progress'])
            assert d['stage']==phase and d['progress']==progress
            row=dict(stage=phases[phase],stage_index=phase,progress=progress,factor=state['boundary_factor'])
            sink(row,native_fields(model,state,d,report));count+=1
    assert count==len(phases)*steps
    # The observer's stop/seal/continue must preserve the whole direct solve.
    process=invoke([probe,*inputs(folder),folder/'whole','--whole'],folder/'whole.log',timeout)
    write(folder/'whole-process.json',process);assert process['exit_code']==0 and not process['timed_out']
    assert exact(read(folder/'whole/run.json'))==exact(native)
    write(folder/'native-verification.json',dict(complete=True,accepted_states=count,all_input_grid_points=True,
        direct_unobserved_run_exact_except_elapsed=True,solver_sha256=sha(probe),process_id=result['pid']))


def independent(folder,steps,user_pc=False):
    model=Model(*(read(p) for p in inputs(folder)),user_pc=user_pc)
    with (folder/'independent-fields.jsonl').open('x') as stream:
        result=run(model,steps,Sink(stream))
    write(folder/'independent-result.json',result)
    write(folder/'independent-process.json',dict(pid=os.getpid(),uses_native_results=False,inputs_sha256={p.name:sha(p) for p in inputs(folder)}))


def restart(folder,probe,steps,timeout):
    root=folder/'restarts';root.mkdir(exist_ok=False);raw=folder/'native';baseline=read(raw/'run.json')
    basepid=read(folder/'native-process.json')['pid'];rows=[]
    names=['geostatic']+[s['name'] for s in baseline['case']['stages']]
    for name in ('birth','open','rebirth'):
        phase=names.index(name);number=phase*steps+steps//2
        cp=raw/f'checkpoint-{number:04}.json';assert cp.is_file()
        saved=read(cp);assert saved['run']['last_committed']['stage_progress']==.5
        suffix=root/name;process=invoke([probe,*inputs(folder),suffix,'--resume',cp],root/(name+'.log'),timeout)
        write(root/(name+'-process.json'),process)
        assert process['exit_code']==0 and not process['timed_out'] and process['pid']!=basepid
        restored=read(suffix/'run.json');assert restored['complete'] and exact(restored)==exact(baseline)
        with (raw/'accepted.jsonl').open() as a,(suffix/'accepted.jsonl').open() as b:
            for _ in range(number):next(a)
            count=0
            for x,y in itertools.zip_longest(a,b):
                assert x is not None and y is not None and exact(json.loads(x))==exact(json.loads(y));count+=1
        assert count==len(names)*steps-number
        rows.append(dict(stage=name,accepted_number=number,saved_progress=.5,checkpoint_sha256=sha(cp),
                         baseline_pid=basepid,resume_pid=process['pid'],suffix_states_exact=count,whole_run_exact_except_elapsed=True))
    write(root/'verification.json',dict(passed=True,actual_files=True,separate_processes=True,checks=rows))


def refusals(out,solver,timeout):
    from verify_joint_ext_combos import tip_fixture
    root=out/'refusals';root.mkdir(exist_ok=False);rows=[]
    def check(name,values,expected):
        folder=root/name;folder.mkdir()
        for p,v in zip(inputs(folder),values):write(p,v)
        process=invoke([solver,'solve-stages',*inputs(folder),folder/'run.json'],folder/'cli.log',timeout)
        write(folder/'process.json',process);text=(folder/'cli.log').read_text()
        assert not process['timed_out'] and process['exit_code']==1 and expected in text,(name,text)
        assert not (folder/'run.json').exists(),(name,'refusal unexpectedly produced result')
        rows.append(dict(name=name,expected_message=expected,exit_code=1,result_file_exists=False,pid=process['pid']))
    for policy in ('transfer_force','stress_free'):
        for water in ('continuous','soil','prescribed'):
            mesh,g,c=tip_fixture(policy);g['initialization']['water_level_m']=1.5
            c['joint_system']['joints'][0]['water']=dict(kind=water)
            if water=='soil':
                c['joint_system']['joints'][0]['water']['soil_side']='a'
                c['joint_system']['structure_bodies']=['UpperClay']
            if water=='prescribed':c['joint_system']['joints'][0]['water']['level_m']=1.5
            check('wet-tip-'+policy+'-'+water,(mesh,g,c),'declared-tip JOINT BIRTH is limited to one dry interface')
        mesh,g,c=tip_fixture(policy);original=copy.deepcopy(mesh);n=len(mesh['nodes']);ne=len(mesh['tetrahedra'])
        mesh['nodes']+=(np.asarray(original['nodes'])+[4.,0.,0.]).tolist()
        for k in ('tetrahedra','tetrahedra10'):mesh[k]+=(np.asarray(original[k])+n).tolist()
        for k in ('material_ids','region_ids'):mesh[k]+=original[k]
        sides={'Seam:a','Seam:b'}
        for key in ('face_sets','node_sets'):
            for name,value in original[key].items():
                shifted=(np.asarray(value)+n).tolist()
                if name in sides:mesh[key][name.replace('Seam','Second')]=shifted
                else:mesh[key][name]+=shifted
        for r in mesh['regions']:r['bounds_m'][1][0]+=4.;r['volume_m3']*=2
        second=copy.deepcopy(c['joint_system']['joints'][0]);second.update(id='Second',side_a='Second:a',side_b='Second:b')
        tip=second['crack_tip'];tip['edges']=(np.asarray(tip['edges'])+n).tolist()
        for f in tip['faces']:
            for side in ('a','b'):f['face_'+side]=[x+n for x in f['face_'+side]];f['owner_'+side]+=ne
        c['joint_system']['joints'].append(second)
        for stage in c['stages']:
            if 'joints' in stage:stage['joints'].append('Second')
        check('multiple-tips-'+policy,(mesh,g,c),'declared-tip JOINT BIRTH is limited to one dry interface')
    for kind in ('beam','plate','anchor'):
        mesh,g,c=tip_fixture('stress_free');_,_,other=structure(kind)
        c['anchor_system' if kind=='anchor' else 'beam_system']=other['anchor_system' if kind=='anchor' else 'beam_system']
        if kind!='anchor':
            for node in c['beam_system']['nodes']:
                if node.get('soil_node_set'):
                    found=np.flatnonzero(np.linalg.norm(np.asarray(mesh['nodes'])-node['position_m'],axis=1)<1e-12).tolist()
                    preferred=[n for n in found if n in mesh['node_sets']['Seam:b']]
                    if preferred:found=preferred
                    assert found
                    mesh['node_sets'][node['soil_node_set']]=[found[0]]
        check('tip-'+kind,(mesh,g,c),'declared-tip JOINT BIRTH with structural members is unverified')
    mesh,g,c=tip_fixture('stress_free');c['joint_system']['joints'][0].pop('crack_tip')
    check('missing-tip',(mesh,g,c),'JOINT interpolation, declared-tip sharing or owner pair mismatch')
    mesh,g,c=tip_fixture('stress_free');c['joint_system']['joints'][0]['crack_tip']['edges'][0][2]=0
    check('inconsistent-tip',(mesh,g,c),'crack-tip')
    mesh,g,c=specimen();c['joint_system']['joints'][0].pop('birth_policy')
    check('missing-birth-policy',(mesh,g,c),'JOINT BIRTH/DEATH requires unique ids and an explicit birth_policy')
    mesh,g,c=specimen();c['joint_system']['joints'].append(copy.deepcopy(c['joint_system']['joints'][0]))
    check('duplicate-joint',(mesh,g,c),'JOINT')
    mesh,g,c=specimen();c['stages'][0]['joints'].append('Seam')
    check('duplicate-birth-id',(mesh,g,c),'JOINT BIRTH/DEATH requires unique ids and an explicit birth_policy')
    mesh,g,c=specimen();g['initialization']['water_level_m']=1.5
    c['stages'].insert(0,dict(kind='hydrostatic_water',name='same_explicit_water',target_level_m=1.5))
    check('missing-water-policy',(mesh,g,c),'staged water with JOINT needs an explicit hydraulic interface policy')
    mesh,g,c=specimen();g['initialization']['water_level_m']=1.5
    c['joint_system']['joints'][0]['water']=dict(kind='prescribed',level_m=1.4)
    check('inconsistent-prescribed-level',(mesh,g,c),'prescribed JOINT pressure conflicts with the adjacent hydrostatic soil')
    write(root/'verification.json',dict(passed=True,checks=rows,unsupported_are_not_accepted=True))


def all_cases(out,probe,solver,steps,timeout,quadratic=False,user_pc=False,selected=None):
    out.mkdir(parents=True,exist_ok=False)
    binaries=out/'_bin';binaries.mkdir()
    frozen_probe=binaries/('probe.exe' if probe.suffix=='.exe' else 'probe')
    frozen_solver=binaries/('solver.exe' if solver.suffix=='.exe' else 'solver')
    shutil.copy2(probe,frozen_probe);shutil.copy2(solver,frozen_solver)
    probe,solver=frozen_probe,frozen_solver
    source_paths=sorted(set(list((ROOT/'src').rglob('*.rs'))+list((ROOT/'tools').rglob('*.py'))+
        list(HERE.glob('*.py'))+[ROOT/'Cargo.toml',ROOT/'Cargo.lock',ROOT/'rust-toolchain.toml',
                               ROOT/'examples/joint_combo_probe.rs',ROOT/'tests/test_joint_combo_extension.py']))
    fingerprints={p.relative_to(ROOT).as_posix():sha(p) for p in source_paths}
    write(out/'source-fingerprints.json',fingerprints)
    with zipfile.ZipFile(out/'source-snapshot.zip','x',compression=zipfile.ZIP_DEFLATED) as z:
        for p in source_paths:z.write(p,p.relative_to(ROOT).as_posix())
    identity=dict(base=BASE,label='2026-11-23',executed_jst='2026-11-12',
        probe_sha256=sha(probe),solver_sha256=sha(solver),steps=steps,quadratic=quadratic,
        source_snapshot_sha256=sha(out/'source-snapshot.zip'),reference_inputs_only=True,
        production_source_changed=False,C01_C04='未達・今回対象外',A01='未達')
    write(out/'identity.json',identity);failures=[];comparisons={}
    for name,values in cases(steps,quadratic).items():
        if selected and name not in selected:continue
        folder=out/name;folder.mkdir()
        for p,v in zip(inputs(folder),values):write(p,v)
        print('case',name,flush=True)
        # Independent fresh solve always runs from inputs even if production
        # refuses/fails. Neither completed nor partial native JSON is an input.
        for mode in ('independent','native','compare','restart'):
            try:
                if mode=='independent':independent(folder,steps,user_pc)
                elif mode=='native':native(folder,probe,steps,timeout)
                elif mode=='compare':
                    result=compare(folder/'native-fields.jsonl',folder/'independent-fields.jsonl')
                    comparisons[name]=result;write(folder/'comparison.json',result)
                elif mode=='restart':restart(folder,probe,steps,timeout)
            except Exception as error:
                failure=dict(case=name,mode=mode,passed=False,error=repr(error),traceback=traceback.format_exc())
                write(folder/('failure-'+mode+'.json'),failure);failures.append(failure)
                print('failure',name,mode,repr(error),flush=True)
    try:refusals(out,solver,timeout)
    except Exception as error:
        failure=dict(mode='refusals',passed=False,error=repr(error),traceback=traceback.format_exc());write(out/'failure-refusals.json',failure);failures.append(failure)
    write(out/'measurements.json',dict(A01='未達',passed=False,requested_all_combinations_passed=False,
        comparisons=comparisons,execution_failures=failures,all_field_thresholds_not_invented=True,
        C01_C04='未達・対象外',NC01='未達・対象外',P6='未達'))
    write(out/'completion.json',dict(measurement_complete=not failures,passed=False,A01='未達',failures=len(failures)))
    assert fingerprints=={p.relative_to(ROOT).as_posix():sha(p) for p in source_paths},'source changed during measurements'
    return int(bool(failures))


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--steps',type=int,default=16)
    p.add_argument('--probe',type=Path,default=ROOT/'target/release/examples/joint_combo_probe')
    p.add_argument('--solver',type=Path,default=ROOT/'target/release/soil-fem-3d');p.add_argument('--timeout',type=int,default=900)
    p.add_argument('--quadratic',action='store_true');p.add_argument('--user-pc',action='store_true')
    p.add_argument('--case',action='append',choices=list(cases(2)),dest='selected')
    args=p.parse_args();return all_cases(args.output.resolve(),args.probe.resolve(),args.solver.resolve(),args.steps,args.timeout,args.quadratic,args.user_pc,args.selected)


if __name__=='__main__':sys.exit(main())
