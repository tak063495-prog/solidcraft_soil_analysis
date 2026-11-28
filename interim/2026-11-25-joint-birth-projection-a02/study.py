"""Fresh bounded solves and real-file process recovery; never replay old JSON as a solve."""
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
import zipfile

for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools'))
import numpy as np
from verify_joint_ext_combos import tip_fixture
from reference import Model, run
from fields import Sink, native_fields, compare, load_fields


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def write(path, value):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def exact(value):
    def strip(v):
        if isinstance(v, dict):
            return {k:strip(x) for k,x in v.items() if k != 'elapsed_seconds'}
        if isinstance(v, list):
            return [strip(x) for x in v]
        return v
    return json.dumps(strip(value), sort_keys=True, separators=(',', ':'), allow_nan=False)


def invoke(command, log, timeout=300):
    timeout=int(os.environ.get('A02_PROCESS_BUDGET_SECONDS',timeout))
    with Path(log).open('x', encoding='utf-8') as stream:
        stream.write(repr(list(map(str, command)))+'\n'); stream.flush()
        started = time.monotonic()
        p = subprocess.Popen(list(map(str, command)), stdout=stream, stderr=subprocess.STDOUT)
        timed_out = False
        try:
            code = p.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True; p.kill(); code = p.wait()
        result = dict(pid=p.pid, exit_code=code, timed_out=timed_out, elapsed_seconds=time.monotonic()-started)
        stream.write('\n'+json.dumps(result)+'\n')
    return result


def fixture(policy, steps, cycle=False):
    m,g,c = tip_fixture(policy)
    c['stages'] = [dict(kind='joint_birth',name='birth',joints=['Seam']),
        dict(kind='boundary',name='compress',target_factor=-.002),
        dict(kind='boundary',name='unload',target_factor=-.001),
        dict(kind='joint_death',name='death',joints=['Seam']),
        dict(kind='joint_birth',name='rebirth',joints=['Seam']),
        dict(kind='boundary',name='reload',target_factor=-.0015),
        dict(kind='boundary',name='return',target_factor=-.001)]
    if cycle:
        c['stages'] = [dict(kind='joint_birth',name='birth',joints=['Seam']),
            dict(kind='boundary',name='open',target_factor=1.),
            dict(kind='boundary',name='recontact',target_factor=-.001),
            dict(kind='joint_death',name='death',joints=['Seam']),
            dict(kind='joint_birth',name='rebirth',joints=['Seam']),
            dict(kind='boundary',name='return',target_factor=-.002)]
    c['settings']['initial_increment'] = c['settings']['minimum_increment'] = 1/steps
    for stage in c['stages']:
        stage['initial_increment'] = 1/steps
    return m,g,c


def input_paths(folder):
    return [folder/(name+'.json') for name in ('mesh','ground','case')]


def make_inputs(folder, values):
    folder.mkdir()
    for path,value in zip(input_paths(folder), values):
        write(path,value)


def native(folder, probe, steps):
    values = [read(p) for p in input_paths(folder)]
    process = invoke([probe,*input_paths(folder),folder/'native'], folder/'native.log')
    write(folder/'native-process.json',process)
    assert process['exit_code'] == 0 and not process['timed_out'], process
    result = read(folder/'native/run.json'); assert result['complete']
    groups = [result['geostatic']['attempts']]+[[a for a in result['attempts'] if a['stage']==s['name']] for s in values[2]['stages']]
    for attempts in groups:
        assert all(a['accepted'] for a in attempts)
        assert [a['progress_to'] for a in attempts] == [i/steps for i in range(1,steps+1)]
        assert all(a['relative_force_residual'] <= values[2]['settings']['relative_force_tolerance'] for a in attempts)
        assert all(x['true_relative_residual'] <= 1e-10 for a in attempts for x in a['linear_solves'])
    model=Model(*values); names=['geostatic']+[s['name'] for s in values[2]['stages']]
    with (folder/'native-fields.jsonl').open('x') as stream:
        sink=Sink(stream)
        with (folder/'native/accepted.jsonl').open() as a, (folder/'native/constitutive.jsonl').open() as b, (folder/'native/structure-reports.jsonl').open() as c:
            count=0
            for state,diagnostic,report in itertools.zip_longest(a,b,c):
                assert state is not None and diagnostic is not None and report is not None
                state,diagnostic,report=map(json.loads,(state,diagnostic,report))
                phase=count//steps; progress=(count%steps+1)/steps
                assert state['stage_index']==phase and diagnostic['stage']==phase and diagnostic['progress']==progress
                assert all(state['active_elements'])
                sink(dict(stage=names[phase],progress=progress),native_fields(model,state,diagnostic,report))
                count+=1
        assert count==len(names)*steps
    direct=invoke([probe,*input_paths(folder),folder/'whole','--whole'],folder/'whole.log')
    write(folder/'whole-process.json',direct)
    assert direct['exit_code']==0 and exact(read(folder/'whole/run.json'))==exact(result)
    # Every preparation starts history at exact zero. Archive values are never reintegrated.
    for audit in result['transition_audits']:
        if audit['stage'] in ('birth','rebirth'):
            for point in audit['prepared_state']['joints']['points']:
                for key in ('plastic_slip_m','accumulated_plastic_slip_m','plastic_dissipation_kpa_m','numerical_dissipation_kpa_m','endpoint_work_kpa_m'):
                    assert np.all(np.asarray(point[key])==0)
    death=next(i for i,s in enumerate(values[2]['stages']) if s['kind']=='joint_death')
    archive=result['last_committed']['joints']['lifecycle'][0]['archived'][0]['points']
    assert exact(archive)==exact(list(enumerate(result['frames'][death]['state']['joints']['points'])))
    return result


def restart(folder, probe, steps):
    baseline=read(folder/'native/run.json'); basepid=read(folder/'native-process.json')['pid']
    names=['geostatic']+[s['name'] for s in baseline['case']['stages']]
    out=folder/'restarts';out.mkdir();checks=[]
    for name in ('geostatic','birth','death','rebirth'):
        number=names.index(name)*steps+steps//2
        cp=folder/f'native/checkpoint-{number:04}.json';assert cp.is_file()
        process=invoke([probe,*input_paths(folder),out/name,'--resume',cp],out/(name+'.log'))
        write(out/(name+'-process.json'),process)
        assert process['exit_code']==0 and not process['timed_out'] and process['pid']!=basepid
        assert exact(read(out/name/'run.json'))==exact(baseline)
        with (folder/'native/accepted.jsonl').open() as a,(out/name/'accepted.jsonl').open() as b:
            for _ in range(number):next(a)
            count=0
            for x,y in itertools.zip_longest(a,b):
                assert x is not None and y is not None and exact(json.loads(x))==exact(json.loads(y));count+=1
        assert count==len(names)*steps-number
        checks.append(dict(stage=name,accepted_number=number,checkpoint_sha256=sha(cp),baseline_pid=basepid,
            resume_pid=process['pid'],suffix_states_exact=count,whole_run_exact_except_elapsed=True))
    return checks


def refuse_inputs(out, solver, steps):
    folder=out/'refusals';folder.mkdir();rows=[]
    def check(name,values,expected):
        f=folder/name;make_inputs(f,values)
        result=invoke([solver,'solve-stages',*input_paths(f),f/'run.json'],f/'cli.log')
        write(f/'process.json',result)
        text=(f/'cli.log').read_text()
        assert result['exit_code']==1 and not result['timed_out'] and expected in text,(name,text)
        assert not (f/'run.json').exists()
        rows.append(dict(name=name,expected=expected,process=result,result_file_exists=False))
    m,g,c=fixture('transfer_force',steps)
    bad=copy.deepcopy(c);bad['joint_system']['joints'][0]['crack_tip']['edges'][0][0]=0
    check('inconsistent-tip',(m,g,bad),'invalid crack-tip')
    bad=copy.deepcopy(c);bad['joint_system']['joints'][0].pop('crack_tip')
    check('hidden-shared-tip',(m,g,bad),'declared-tip sharing or owner pair mismatch')
    # Physical pre-birth forces from an equilibrated input, not requested JSON traction.
    for name,factor in (('tensile-transfer',1.),('outside-circle-transfer',-.001)):
        gg=copy.deepcopy(g);cc=copy.deepcopy(c);cc['stages']=cc['stages'][:1]
        gg['initialization']['geostatic']['initial_boundary_factor']=factor
        if name.startswith('outside'):
            cc['boundary']['constraints'][3]['value_m']=10.
        check(name,(m,gg,cc),'violates no-tension/Coulomb cone')
    for policy in ('transfer_force','stress_free'):
        mm,gg,cc=fixture(policy,steps)
        cc['initially_inactive_bodies']=['UpperClay']
        cc['stages']=cc['stages'][:1]
        check('inactive-owner-'+policy,(mm,gg,cc),'both active bodies')
        fresh=copy.deepcopy(gg['initialization']['geostatic']['seed']['states']['UpperClay'])
        fresh['effective_stress_compression_positive_kpa']=[0.]*6
        cc['stages'].insert(0,dict(kind='birth',name='body_fill',bodies=['UpperClay'],
            states={'UpperClay':fresh},placement='cad_new_nodes'))
        check('body-across-tip-'+policy,(mm,gg,cc),'body BIRTH across a declared crack tip is unverified')
    # Existing excavation-set input: a connected 3-Tet cut reaches the top.
    # It removes owner 13 while owner 12 of the other Tri6 face stays active.
    # No buried-cavity guard is bypassed and geometry/material are unchanged.
    mm,gg,cc=fixture('stress_free',steps)
    assert [f['owner_b'] for f in cc['joint_system']['joints'][0]['crack_tip']['faces']]==[13,12]
    mm['excavation_sets']={'PartialFace':[13,14,15]}
    gg['excavation_bodies']=['PartialFace']
    cc['stages']=[dict(kind='joint_birth',name='birth',joints=['Seam']),
        dict(kind='joint_death',name='death',joints=['Seam']),
        dict(kind='excavate',name='remove_one_owner',bodies=['PartialFace'],remaining_surface_pressure_kpa=0.),
        dict(kind='joint_birth',name='partial_rebirth',joints=['Seam'])]
    check('partial-face-rebirth',(mm,gg,cc),'both active bodies')
    partial=folder/'partial-face-rebirth';prefix=partial/'excavation-prefix.json'
    process=invoke([solver,'checkpoint-stages',*input_paths(partial),4*steps-1,prefix],partial/'prefix.log')
    assert process['exit_code']==0 and not process['timed_out']
    state=read(prefix)['run']['last_committed']
    assert not state['active_elements'][13] and state['active_elements'][12]
    life=state['joints']['lifecycle'][0]
    assert life['generation']==1 and not life['active'] and life['tied'] and len(life['archived'])==1
    rows[-1]['actual_partial_prefix']=dict(process=process,checkpoint_sha256=sha(prefix),
        inactive_owner=13,active_owner=12,generation=1,frozen_generations=1)
    return rows


def projection_requests(out, probe, steps):
    folder=out/'projection';folder.mkdir();m,g,c=fixture('transfer_force',steps)
    write(folder/'mesh.json',m);write(folder/'case.json',c)
    model=Model(m,g,c);s=model.system
    def force(q):return np.einsum('qag,a,q->g',s.c,q,s.weight)
    good=force([12.,5.,-100.]);tip=c['joint_system']['joints'][0]['crack_tip']['edges'][0]
    unrep=good.copy();direction=np.asarray(m['nodes'][tip[1]])-m['nodes'][tip[0]]
    unrep[3*tip[0]:3*tip[0]+3]+=direction;unrep[3*tip[1]:3*tip[1]+3]-=direction
    delta=(unrep-good).reshape(-1,3)
    assert np.max(np.abs(delta.sum(axis=0)))<1e-10
    assert np.max(np.abs(np.cross(m['nodes'],delta).sum(axis=0)))<1e-10
    requests=[('good',good,None,True,None),('tension',force([0.,0.,100.]),None,False,'violates no-tension/Coulomb cone'),
        ('outside-circle',force([40.,0.,-100.]),None,False,'violates no-tension/Coulomb cone'),
        ('balanced-unreproducible-tip',unrep,None,False,'every free nodal residual'),
        ('partial-rebirth',good,s.points[0]['owner_b'],False,'both active bodies'),
        ('body-across-tip',good,None,False,'body BIRTH across a declared crack tip')]
    rows=[]
    for name,rhs,owner,accepted,message in requests:
        operation={'partial-rebirth':'partial_rebirth','body-across-tip':'body_birth'}.get(name,'birth')
        request=dict(operation=operation,residual_kn=rhs.tolist(),inactive_owner=owner)
        write(folder/(name+'-request.json'),request)
        result=invoke([probe,folder/'mesh.json',folder/'case.json',folder/(name+'-request.json'),folder/(name+'-result.json')],folder/(name+'.log'))
        response=read(folder/(name+'-result.json'))
        assert response['accepted']==accepted and result['exit_code']==(0 if accepted else 1) and not result['timed_out']
        if message:assert message in response['error'],response['error']
        rows.append(dict(name=name,diagnostic_public_api=True,accepted=accepted,error=response['error'],process=result))
    # Independent weighted projection, using the same input force and no native field.
    matrix=np.einsum('qag,q->gqa',s.c,np.sqrt(s.weight)).reshape(model.nd,-1)
    projected=np.linalg.lstsq(matrix,good,rcond=1e-12)[0]
    assert np.max(np.abs(matrix@projected-good))<1e-10
    expected=projected.reshape(-1,3)/np.sqrt(s.weight)[:,None]
    actual=np.asarray([p['birth']['initial_traction_kpa'] for p in read(folder/'good-result.json')['result']['state']['points']])
    error=float(np.max(np.abs(actual-expected)))
    assert np.linalg.norm(actual-expected)/np.linalg.norm(expected)<1e-7
    return dict(requests=rows,initial_traction_maximum_absolute_kpa=error,balanced_tip_delta_force_moment_zero=True)


def tamper(folder,probe,rehash,steps):
    out=folder/'tamper';out.mkdir(); cp=read(folder/f'native/checkpoint-{steps+steps//2:04}.json');rows=[]
    for mode in ('tension','outside-circle','partial-active','generation','reference','frozen'):
        # Frozen tampering uses a post-DEATH checkpoint.
        base=read(folder/f'native/checkpoint-{4*steps+steps//2:04}.json') if mode=='frozen' else copy.deepcopy(cp)
        state=base['run']['last_committed'];j=state['joints']
        if mode=='tension':j['points'][0]['birth']['initial_traction_kpa']=[0.,0.,100.]
        elif mode=='outside-circle':j['points'][0]['birth']['initial_traction_kpa']=[40.,0.,-100.]
        elif mode=='partial-active':state['active_elements'][3]=False
        elif mode=='generation':j['lifecycle'][0]['generation']+=1
        elif mode=='reference':j['points'][0]['birth']['jump_m'][0]+=.001
        else:j['lifecycle'][0]['archived'][0]['points'][0][1]['traction_kpa'][0]+=1.
        raw=out/(mode+'-raw.json');sealed=out/(mode+'-rehashed.json');write(raw,base)
        result=invoke([rehash,raw,sealed],out/(mode+'-rehash.log'));assert result['exit_code']==0
        assert read(sealed)['identity']['solver_binary_sha256']==sha(probe)
        result=invoke([probe,*input_paths(folder),out/(mode+'-resume'),'--resume',sealed],out/(mode+'.log'))
        text=(out/(mode+'.log')).read_text()
        assert result['exit_code']==1 and not result['timed_out'] and not (out/(mode+'-resume/run.json')).exists()
        assert 'checkpoint accepted-path replay mismatch' in text,text
        assert 'input/solver hash mismatch' not in text,text
        rows.append(dict(mode=mode,outer_hash_recomputed=True,checkpoint_identity_matches_producer=True,
            producer_binary_sha256=sha(probe),process=result,expected_error='checkpoint accepted-path replay mismatch',error=text))
    return rows


def audits(folder):
    model=Model(*(read(p) for p in input_paths(folder)));rows={}
    for name in ('native','independent'):
        work_error=0.;force_residual=0.;count=0
        for row,f in load_fields(folder/(name+'-fields.jsonl')):
            balance=f['initial_energy']+f['endpoint_work']-f['elastic_energy']-f['plastic_dissipation']-f['numerical_dissipation']
            work_error=max(work_error,float(np.max(np.abs(balance))))
            t,_=model.spaces[bool(f['joint_tied'][0])]
            r=t.T@(f['internal']-f['applied_load']).ravel()
            force_residual=max(force_residual,float(np.linalg.norm(r)/max(np.linalg.norm(f['applied_load']),1.)))
            count+=1
        assert work_error<1e-10 and force_residual<1e-10
        rows[name]=dict(accepted_states=count,maximum_contact_work_identity_error_kpa_m=work_error,
            maximum_reassembled_free_force_relative_residual=force_residual,
            force_source='committed-state public production kernel reassembly' if name=='native' else 'input-only fresh equilibrium')
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--steps',type=int,default=4)
    p.add_argument('--independent',type=Path);args=p.parse_args()
    if args.independent:
        folder=args.independent
        with (folder/'independent-fields.jsonl').open('x') as stream:
            result=run(Model(*(read(x) for x in input_paths(folder))),args.steps,Sink(stream))
        write(folder/'independent-result.json',result)
        write(folder/'independent-process.json',dict(pid=os.getpid(),uses_native_results=False,input_sha256={x.name:sha(x) for x in input_paths(folder)}))
        return
    assert args.steps>=2 and args.steps%2==0
    out=args.output.resolve();out.mkdir(parents=True);binaries=out/'bin';binaries.mkdir()
    sources=[]
    for name in ('src','tests','tools','examples','.github','interim'):
        sources.extend(p for p in (ROOT/name).rglob('*') if p.is_file() and p.suffix in ('.py','.rs')
            and '__pycache__' not in p.parts)
    sources += [ROOT/name for name in ('Cargo.toml','Cargo.lock','rust-toolchain.toml','requirements-reference.lock.txt')]
    sources += [Path(__file__).parent/'SPEC.ja.md']
    fingerprint={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(set(sources))}
    write(out/'source-fingerprints.json',fingerprint)
    with zipfile.ZipFile(out/'source-snapshot.zip','x',compression=zipfile.ZIP_DEFLATED) as z:
        for name,digest in fingerprint.items():
            assert sha(ROOT/name)==digest
            z.writestr(name,(ROOT/name).read_bytes())
    exe='.exe' if os.name=='nt' else ''
    sources={'solver':ROOT/'target/release'/('soil-fem-3d'+exe),
        'probe':ROOT/'target/release/examples'/('joint_a02_history_probe'+exe),
        'projection':ROOT/'target/release/examples'/('joint_a02_projection_probe'+exe),
        'rehash':ROOT/'target/release/examples'/('joint_ext_checkpoint_rehash'+exe)}
    fixed={}
    for name,path in sources.items():
        fixed[name]=binaries/(name+exe);shutil.copy2(path,fixed[name])
    write(out/'binaries.json',{name:dict(original=str(sources[name]),sha256=sha(path)) for name,path in fixed.items()})
    summary=dict(passed=False,a02_complete=False,base='702a19a7ce189f1173caed727ee2a0059298a607',cases={})
    for name,policy,cycle in (('transfer_force','transfer_force',False),('stress_free','stress_free',False),('stress_free-cycle','stress_free',True)):
        folder=out/name;make_inputs(folder,fixture(policy,args.steps,cycle));native(folder,fixed['probe'],args.steps)
        result=invoke([sys.executable,Path(__file__).resolve(),'--output',out,'--steps',args.steps,'--independent',folder],folder/'independent.log')
        assert result['exit_code']==0 and not result['timed_out']
        comparison=compare(folder/'native-fields.jsonl',folder/'independent-fields.jsonl');write(folder/'comparison.json',comparison)
        fields=comparison['all_fields']
        assert all(fields[k]['relative_l2']<1e-7 for k in ('u','stress','traction'))
        assert all(fields[k]['maximum_absolute']<1e-10 for k in ('internal','applied_load'))
        assert sum(v.get('mismatches',0) for v in fields.values())==0
        summary['cases'][name]=dict(comparison=comparison,existing_bounded_gates_passed=True,
            audits=audits(folder),restarts=restart(folder,fixed['probe'],args.steps),tampering=tamper(folder,fixed['probe'],fixed['rehash'],args.steps))
        print('completed',name,flush=True)
    summary['refusals']=refuse_inputs(out,fixed['solver'],args.steps)
    summary['projection']=projection_requests(out,fixed['projection'],args.steps)
    summary['scope_judgments']={
        'dry_single_tip_all_owners_active_contact_death_rebirth':'existing limited acceptance reconfirmed; no scope expansion',
        'body_birth_across_tip':'unmet / unsupported; current refusal retained',
        'partially_active_face_rebirth':'unmet / unsupported; current refusal retained',
        'tensile_initial_force_projection':'unmet / unsupported; current refusal retained',
        'outside_friction_circle_initial_force_projection':'unmet / unsupported; current refusal retained',
        'unreproducible_initial_force_projection':'unmet / unsupported; current refusal retained'}
    summary['fresh_input_only_solve_pairs']=len(summary['cases'])
    summary['full_field_accepted_pairs']=sum(v['comparison']['accepted_pairs'] for v in summary['cases'].values())
    assert all(sha(ROOT/name)==digest for name,digest in fingerprint.items())
    summary['source_count']=len(fingerprint)
    summary['sources_unchanged_during_fresh_solves']=True
    write(out/'measurements.json',summary)


if __name__=='__main__':main()
