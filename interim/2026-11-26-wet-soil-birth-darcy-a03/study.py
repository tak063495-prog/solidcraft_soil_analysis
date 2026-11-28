"""Fresh 3D solves, actual CLI refusals and different-PID real-file recovery."""
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
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
import numpy as np
from verify_joint_ext_combos import wet_fixture
import fem_reference as fem
from reference import Model, solve
from fields import Sink, native_fields, compare, row


def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def write(path,value):
    with Path(path).open('x',encoding='utf-8') as f:json.dump(value,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def exact(value):
    def strip(v):
        if isinstance(v,dict):return {k:strip(x) for k,x in v.items() if k!='elapsed_seconds'}
        if isinstance(v,list):return [strip(x) for x in v]
        return v
    return json.dumps(strip(value),sort_keys=True,separators=(',',':'),allow_nan=False)
def invoke(command,log,timeout=300):
    timeout=int(os.environ.get('A03_PROCESS_BUDGET_SECONDS',timeout))
    with Path(log).open('x',encoding='utf-8') as f:
        f.write(repr(list(map(str,command)))+'\n');f.flush();start=time.monotonic()
        p=subprocess.Popen(list(map(str,command)),stdout=f,stderr=subprocess.STDOUT)
        timed_out=False
        try:code=p.wait(timeout=timeout)
        except subprocess.TimeoutExpired:timed_out=True;p.kill();code=p.wait()
        result=dict(pid=p.pid,exit_code=code,timed_out=timed_out,elapsed_seconds=time.monotonic()-start)
        f.write('\n'+json.dumps(result)+'\n')
    return result

def fixture(steps,level=1.5):
    m,g,c=wet_fixture();c['stages']=c['stages'][:3];c['verification_reject_once_stages']=[]
    c['settings']['kinematics']='side_rollers';c.pop('boundary',None)
    c['settings']['relative_force_tolerance']=1e-11;c['settings']['relative_linear_tolerance']=1e-12
    g['initialization']['water_level_m']=level;c['stages'][0]['hydraulic_initialization']['level_m']=level
    if level==2.:g['initialization']['geostatic']['hydrostatic_face_sets']=['exterior']
    c['settings']['initial_increment']=c['settings']['minimum_increment']=1/steps
    for s in c['stages']:s['initial_increment']=1/steps
    return m,g,c

def soil_role_fixture(steps):
    # One existing lower-soil Tet is assigned to a new, initially inactive body.
    # SOIL keeps lower old soil against the impermeable upper MC continuum.
    # No contact face owner is selected for the new body; no hydraulic role is faked.
    m,g,c=fixture(steps)
    region=copy.deepcopy(m['regions'][0]);region['id']=3;region['body']='FreshSoil'
    region['volume_m3']=1/12;region['source_sha256']=hashlib.sha256(json.dumps([m['nodes'][n] for n in m['tetrahedra'][0]]).encode()).hexdigest()
    m['regions'][0]['volume_m3']-=region['volume_m3'];m['regions'].append(region);m['region_ids'][0]=3
    g['bodies'].append(dict(body='FreshSoil',material_id=1))
    g['initialization']['geostatic']['seed']['states']['FreshSoil']=copy.deepcopy(g['initialization']['geostatic']['seed']['states']['LowerClay'])
    c['initially_inactive_bodies']=['FreshSoil'];birth=c['stages'][0];birth['bodies']=['FreshSoil'];birth['states']={'FreshSoil':birth['states']['UpperClay']}
    c['joint_system']['structure_bodies']=['UpperClay'];c['joint_system']['joints'][0]['water']=dict(kind='soil',soil_side='a')
    return m,g,c


def paths(folder):return [folder/(n+'.json') for n in ('mesh','ground','case')]
def inputs(folder,values):
    folder.mkdir()
    for p,v in zip(paths(folder),values):write(p,v)

def independent(folder,steps):
    m=Model(*[read(p) for p in paths(folder)])
    with (folder/'independent-fields.jsonl').open('x') as f,(folder/'independent-prepared.jsonl').open('x') as p:
        result=solve(m,steps,Sink(f),Sink(p,work=False))
    write(folder/'independent-result.json',result)

def geometry(m,run,folder):
    expected=[]
    for e,q in enumerate(m.sample):
        for i in range(len(q['weight'])):
            expected.append(dict(element_index=e,local_index=i,natural_coordinates=q['natural'][i].tolist(),position_m=q['position'][i].tolist(),weight_m3=float(q['weight'][i]),shape_functions=q['shape'][i].tolist()))
    assert len(expected)==len(run['integration_points'])
    errors={k:0. for k in ('position_m','weight_m3','shape_functions','natural_coordinates')}
    for a,b in zip(run['integration_points'],expected):
        assert (a['element_index'],a['local_index'])==(b['element_index'],b['local_index'])
        for k in errors:errors[k]=max(errors[k],float(np.max(np.abs(np.asarray(a[k])-b[k]))))
    assert errors['position_m']<1e-10
    contacts=m.system.points;assert len(contacts)==len(run['joint_integration_points'])
    contact_errors={k:0. for k in ('position_m','normal','weight_m2','shape')}
    for a,b in zip(run['joint_integration_points'],contacts):
        for k in ('joint_id','face_index','nodes_a','nodes_b','owner_a','owner_b'):assert a[k]==b[k]
        for k in contact_errors:contact_errors[k]=max(contact_errors[k],float(np.max(np.abs(np.asarray(a[k])-b[k]))))
    assert contact_errors['position_m']<1e-10
    write(folder/'independent-layout.json',dict(volume=expected,contact=contacts))
    return dict(volume_points=len(expected),contact_points=len(contacts),all_ids_and_owners_match=True,volume_maximum_absolute=errors,contact_maximum_absolute=contact_errors)


def native(folder,probe,solver,steps):
    proc=invoke([probe,*paths(folder),folder/'native'],folder/'native.log');write(folder/'native-process.json',proc)
    assert proc['exit_code']==0 and not proc['timed_out'],proc
    run=read(folder/'native/run.json');assert run['complete']
    assert run['provenance']['solver_binary_sha256']==sha(probe)
    assert not run.get('flow_audits')
    m=Model(*[read(p) for p in paths(folder)])
    layout=geometry(m,run,folder)
    birth=m.case['stages'][0];preparation=next(a for a in run['transition_audits'] if a['stage']==birth['name'])
    selected=set(preparation['selected_elements'])
    for point,state in zip(run['integration_points'],preparation['prepared_state']['material_states']):
        if point['element_index'] in selected:
            assert exact(state)==exact(birth['states']['UpperClay']),'fresh input state was modified during preparation'
    assert preparation['zero_release_relative_force_residual']<1e-10
    names=['geostatic']+[s['name'] for s in m.case['stages']]
    groups=[run['geostatic']['attempts']]+[[a for a in run['attempts'] if a['stage']==s['name']] for s in m.case['stages']]
    for attempts in groups:
        assert all(a['accepted'] for a in attempts)
        assert [a['progress_to'] for a in attempts]==[i/steps for i in range(1,steps+1)]
        assert all(a['relative_force_residual']<=m.case['settings']['relative_force_tolerance'] for a in attempts)
        assert all(x['true_relative_residual']<=1e-10 for a in attempts for x in a['linear_solves'])
    audits=[]
    with (folder/'native-fields.jsonl').open('x') as f:
        sink=Sink(f)
        with (folder/'native/accepted.jsonl').open() as a,(folder/'native/constitutive.jsonl').open() as b:
            count=0
            for state,diag in itertools.zip_longest(a,b):
                assert state is not None and diag is not None
                state,diag=map(json.loads,(state,diag));phase=count//steps;fraction=(count%steps+1)/steps
                assert state['stage_index']==phase
                if count==steps:
                    sink.placement(preparation['prepared_state']['displacement_m'])
                fields,audit=native_fields(m,state,diag,run)
                sink(row(names[phase],fraction),fields);audits.append(audit);count+=1
        assert count==len(names)*steps
    with (folder/'native-prepared.jsonl').open('x') as f:
        sink=Sink(f,work=False)
        for prep in read(folder/'native/preparations.json'):
            audit=next(a for a in run['transition_audits'] if a['stage']==prep['stage'])
            fields,check=native_fields(m,audit['prepared_state'],prep['diagnostic'],run,True)
            sink(row(prep['stage'],0.),fields)
            audits.append(dict(preparation=True,**check))
    # Same executable and full provenance; only measured time may differ.
    whole=invoke([probe,*paths(folder),folder/'whole','--whole'],folder/'whole.log');write(folder/'whole-process.json',whole)
    assert whole['exit_code']==0 and not whole['timed_out']
    assert exact(read(folder/'whole/run.json'))==exact(run)
    # Standard production CLI is an additional control, not the restart oracle.
    cli=invoke([solver,'solve-stages',*paths(folder),folder/'production-run.json'],folder/'production.log');write(folder/'production-process.json',cli)
    assert cli['exit_code']==0 and not cli['timed_out']
    standard=read(folder/'production-run.json');assert standard['complete']
    left=copy.deepcopy(run);right=copy.deepcopy(standard)
    for r in (left,right):r['provenance'].pop('solver_binary_sha256')
    assert exact(left)==exact(right)
    write(folder/'native-audits.json',audits)
    return dict(accepted=count,layout=layout,maximum_native_reassembled_relative_residual=max(a['relative_force_residual'] for a in audits),
                maximum_native_load_decomposition_error_kn=max(a['load_decomposition_error_kn'] for a in audits),maximum_constraint_error_m=max(a['constraint_error_m'] for a in audits),maximum_local_work_identity_error_kpa_m=max(a['local_work_identity_error_kpa_m'] for a in audits),whole_exact_except_elapsed=True,
                standard_cli_exact_except_elapsed_and_binary_sha=True)

def restart(folder,probe,steps):
    run=read(folder/'native/run.json');baseline_pid=read(folder/'native-process.json')['pid']
    names=['geostatic']+[s['name'] for s in run['case']['stages']];out=folder/'restarts';out.mkdir();rows=[]
    for phase,name in enumerate(names):
        number=phase*steps+steps//2;cp=folder/f'native/checkpoint-{number:04}.json';assert cp.is_file()
        proc=invoke([probe,*paths(folder),out/name,'--resume',cp],out/(name+'.log'));write(out/(name+'-process.json'),proc)
        assert proc['exit_code']==0 and not proc['timed_out'] and proc['pid']!=baseline_pid
        assert exact(read(out/name/'run.json'))==exact(run)
        count=0
        with (folder/'native/accepted.jsonl').open() as a,(out/name/'accepted.jsonl').open() as b:
            for _ in range(number):next(a)
            for x,y in itertools.zip_longest(a,b):
                assert x is not None and y is not None and exact(json.loads(x))==exact(json.loads(y));count+=1
        assert count==len(names)*steps-number
        rows.append(dict(stage=name,accepted_number=number,checkpoint_sha256=sha(cp),baseline_pid=baseline_pid,resume_pid=proc['pid'],
                         suffix_states_exact=count,whole_run_exact_except_elapsed=True))
    return rows

def refusals(out,solver,steps):
    folder=out/'refusals';folder.mkdir();rows=[]
    def check(name,values,expected):
        f=folder/name;inputs(f,values)
        proc=invoke([solver,'solve-stages',*paths(f),f/'forbidden.json'],f/'cli.log');write(f/'process.json',proc)
        text=(f/'cli.log').read_text()
        assert proc['exit_code']==1 and not proc['timed_out'] and expected in text,(name,text)
        assert not (f/'forbidden.json').exists()
        row=dict(name=name,expected=expected,process=proc,result_file_exists=False);rows.append(row)
        return f,row
    guard='wet BIRTH requires one fully submerged fresh soil body'
    for mode in ('mismatched-level','partial-submerged','partial-disguised','missing-hydraulic','aged-plastic-strain','aged-equivalent','aged-work','aged-return','active-body','placement','mode-A','soil-policy','prescribed-policy','prescribed-conflict','undefined-level'):
        m,g,c=fixture(steps);birth=c['stages'][0];state=birth['states']['UpperClay'];joint=c['joint_system']['joints'][0];expected=guard
        if mode=='mismatched-level':birth['hydraulic_initialization']['level_m']=1.4
        elif mode=='partial-submerged':g['initialization']['water_level_m']=birth['hydraulic_initialization']['level_m']=.75
        elif mode=='partial-disguised':g['initialization']['water_level_m']=.75
        elif mode=='missing-hydraulic':birth.pop('hydraulic_initialization');expected='stress-free BIRTH requires dry, zero-stress, fresh material states'
        elif mode=='aged-plastic-strain':state['plastic_strain_compression_positive'][0]=.001;expected='stress-free BIRTH requires dry, zero-stress, fresh material states'
        elif mode=='aged-equivalent':state['equivalent_plastic_strain']=.1;expected='BIRTH cannot inherit prior plastic history'
        elif mode=='aged-work':state['plastic_work_kpa']=1.;expected='BIRTH cannot inherit prior plastic history'
        elif mode=='aged-return':state['last_return']='shear_face';expected='BIRTH cannot inherit prior plastic history'
        elif mode=='active-body':c['initially_inactive_bodies']=[];expected='BIRTH requires a wholly inactive body'
        elif mode=='placement':birth['placement']='current_displacement'
        elif mode=='mode-A':joint['birth_policy']='transfer_force'
        elif mode=='soil-policy':joint['water']=dict(kind='soil',soil_side='a');expected='SOIL requires the declared soil side and one impermeable structure side'
        elif mode=='prescribed-policy':joint['water']=dict(kind='prescribed',level_m=1.5)
        elif mode=='prescribed-conflict':joint['water']=dict(kind='prescribed',level_m=1.4)
        elif mode=='undefined-level':g['initialization']['water_level_m']=None
        check(mode,(m,g,c),expected)
    m,g,c=fixture(steps)
    mg=read(ROOT/'examples/ground-joint-water-mcc.json')
    upper=next(m for m in mg['materials'] if m['id']==2)
    g['materials']=[copy.deepcopy(upper) if mat['id']==2 else mat for mat in g['materials']]
    fresh=copy.deepcopy(mg['initialization']['geostatic']['seed']['states']['UpperClay'])
    g['initialization']['geostatic']['seed']['states']['UpperClay']=fresh
    c['stages'][0]['states']['UpperClay']=copy.deepcopy(fresh)
    f,row=check('fresh-MCC',(m,g,c),'MCC has zero bulk stiffness at p=0')
    control=copy.deepcopy(c);control['stages']=[];write(f/'control-case.json',control)
    proc=invoke([solver,'solve-stages',f/'mesh.json',f/'ground.json',f/'control-case.json',f/'control-run.json'],f/'control.log')
    write(f/'control-process.json',proc)
    assert proc['exit_code']==0 and not proc['timed_out'] and read(f/'control-run.json')['complete']
    row['valid_MCC_seed_control']=dict(process=proc,complete=True,upper_inactive=True)
    # Invalid BIRTH initialization tags are schema refusals, NOT physical tests.
    for kind in ('soil','continuous','prescribed'):
        m,g,c=fixture(steps);c['stages'][0]['hydraulic_initialization']['kind']=kind
        check('initialization-enum-'+kind,(m,g,c),'unknown variant')
    m,g,c=fixture(steps);c['initially_inactive_bodies']=[];c['stages']=[]
    c['joint_system']['joints'][0]['water']=dict(kind='prescribed',level_m=1.4)
    check('prescribed-pressure-conflict-control',(m,g,c),'prescribed JOINT pressure conflicts with the adjacent hydrostatic soil')
    check('soil-valid-roles-hydrostatic',soil_role_fixture(steps),guard)
    # Flow policies are valid inputs. Require a completed flow-only checkpoint
    # with a real pressure field before interpreting the following BIRTH refusal.
    for policy in ('soil','continuous','prescribed'):
        m,g,c=soil_role_fixture(steps) if policy=='soil' else fixture(steps)
        c['joint_system']['joints'][0]['water']={'kind':policy,**({'soil_side':'a'} if policy=='soil' else {'level_m':1.5} if policy=='prescribed' else {})}
        flow=read(ROOT/'examples/stages-joint-water-flow.json')['stages'][0]['flow']
        if policy=='soil':
            topology=fem.topology(m)
            m['face_sets']['old_soil_bottom']=[f for f in m['face_sets']['bottom'] if all(e!=0 for e in topology[tuple(sorted(f))])]
            assert m['face_sets']['old_soil_bottom']
        flow['heads']=[dict(face_set='old_soil_bottom' if policy=='soil' else 'bottom',head_m=1.5,gradient=[.02,0.,0.],reference_point_m=[0.,0.,0.])]
        c['stages'].insert(0,dict(kind='steady_seepage',name='flow',flow=flow,exterior_water_level_m=1.5,initial_increment=1/steps))
        f,row=check('darcy-before-birth-'+policy,(m,g,c),'BIRTH after seepage needs an explicit initial hydraulic state')
        prefix=copy.deepcopy(c);prefix['stages']=prefix['stages'][:1];write(f/'prefix-case.json',prefix)
        proc=invoke([solver,'checkpoint-stages',*paths(f),2*steps,f/'flow-complete-checkpoint.json'],f/'prefix.log');write(f/'prefix-process.json',proc)
        assert proc['exit_code']==0 and not proc['timed_out'],(policy,proc)
        cp=read(f/'flow-complete-checkpoint.json');state=cp['run']['last_committed']
        assert state['pore_pressure_kpa'] is not None and len(state['pore_pressure_kpa'])>0
        assert not all(state['active_elements']) and cp['run']['flow_audits']
        solution=cp['run']['flow_audits'][-1]['solution']
        maximum_flux=max(float(np.linalg.norm(p['darcy_flux_m_s'])) for p in solution['points'])
        assert solution['complete'] and maximum_flux>1e-10
        row['completed_flow_prefix']=dict(process=proc,checkpoint_sha256=sha(f/'flow-complete-checkpoint.json'),maximum_darcy_flux_m_s=maximum_flux,balances=solution['balances'],
              pore_pressure_count=len(state['pore_pressure_kpa']),flow_audits=len(cp['run']['flow_audits']),birth_has_not_occurred=True)
    return rows

def tamper(folder,probe,rehash,steps):
    out=folder/'tamper';out.mkdir();rows=[]
    original_path=folder/f'native/checkpoint-{2*steps+steps//2:04}.json';original=read(original_path)
    for mode in ('water-pressure','water-load','birth-initial-effective-stress','birth-holding','active-element','aged-material'):
        bad=copy.deepcopy(original);state=bad['run']['last_committed'];audit=bad['run']['transition_audits'][0]
        if mode=='water-pressure':state['joint_water']['physical_pressure_kpa'][0]+=1.
        elif mode=='water-load':state['joint_water']['nodal_load_kn'][0]+=.1
        elif mode=='birth-initial-effective-stress':audit['prepared_state']['material_states'][-1]['effective_stress_compression_positive_kpa'][0]+=1.
        elif mode=='birth-holding':audit['initial_holding_load_kn'][0]+=.1
        elif mode=='active-element':state['active_elements'][-1]=False
        else:state['material_states'][-1]['equivalent_plastic_strain']=.1
        raw=out/(mode+'-raw.json');sealed=out/(mode+'-sealed.json');write(raw,bad)
        proc=invoke([rehash,raw,sealed],out/(mode+'-rehash.log'))
        assert proc['exit_code']==0 and not proc['timed_out']
        assert read(sealed)['run_sha256']!=original['run_sha256']
        proc=invoke([probe,*paths(folder),out/(mode+'-resume'),'--resume',sealed],out/(mode+'.log'));text=(out/(mode+'.log')).read_text()
        assert proc['exit_code']!=0 and not proc['timed_out'] and 'checkpoint accepted-path replay mismatch' in text,(mode,text)
        rows.append(dict(mutation=mode,outer_hash_recomputed=True,same_binary_sha256=sha(probe),process=proc,reason='checkpoint accepted-path replay mismatch'))
    return rows

def snapshot(out):
    files=subprocess.check_output(['git','ls-files','src','tools','tests','examples','Cargo.toml','Cargo.lock','rust-toolchain.toml'],cwd=ROOT,text=True).splitlines()
    files+=['examples/wet_birth_a03_probe.rs']
    files+=[str(p.relative_to(ROOT)) for p in Path(__file__).parent.iterdir() if p.suffix=='.py' or p.name=='SPEC.ja.md']
    files=sorted(set(files));fingerprints={f:sha(ROOT/f) for f in files}
    write(out/'source-fingerprints.json',fingerprints)
    with zipfile.ZipFile(out/'source-snapshot.zip','x',zipfile.ZIP_DEFLATED) as z:
        for f in files:z.write(ROOT/f,f)
    return fingerprints

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--steps',type=int,default=4)
    p.add_argument('--independent',type=Path);a=p.parse_args()
    if a.independent:independent(a.independent,a.steps);return
    assert a.steps>=2 and a.steps%2==0
    out=a.output.resolve();out.mkdir(parents=True);assert not any(out.iterdir())
    fingerprints=snapshot(out);binpath=out/'bin';binpath.mkdir();suffix='.exe' if os.name=='nt' else ''
    binaries={}
    for name,source in [('solver',ROOT/('target/release/soil-fem-3d'+suffix)),('probe',ROOT/('target/release/examples/wet_birth_a03_probe'+suffix)),('rehash',ROOT/('target/release/examples/joint_ext_checkpoint_rehash'+suffix))]:
        dest=binpath/(name+suffix);shutil.copy2(source,dest);binaries[name]=dict(path=str(dest),sha256=sha(dest))
    write(out/'binaries.json',binaries);probe=Path(binaries['probe']['path']);solver=Path(binaries['solver']['path']);rehash=Path(binaries['rehash']['path'])
    report=dict(passed=False,a03_complete=False,scope='existing hydrostatic fully submerged fresh MC / single CONTINUOUS B only',cases={},refusals=[])
    for name,level in [('hydrostatic-h1_5',1.5),('hydrostatic-h2_0',2.)]:
        folder=out/name;inputs(folder,fixture(a.steps,level))
        observed=native(folder,probe,solver,a.steps)
        proc=invoke([sys.executable,Path(__file__),'--output',out,'--steps',a.steps,'--independent',folder],folder/'independent.log');write(folder/'independent-process.json',proc)
        assert proc['exit_code']==0 and not proc['timed_out'],proc
        comparison=compare(folder/'native-fields.jsonl',folder/'independent-fields.jsonl');write(folder/'comparison.json',comparison)
        prepared=compare(folder/'native-prepared.jsonl',folder/'independent-prepared.jsonl');write(folder/'prepared-comparison.json',prepared)
        for key in ('internal','applied_load','effective_internal','pore_internal','holding_load'):
            assert prepared['all_fields'][key]['maximum_absolute']<1e-10,(name,'prepared',key,prepared['all_fields'][key])
        for key in ('u','stress','traction'):
            assert comparison['all_fields'][key]['relative_l2']<1e-7,(name,key,comparison['all_fields'][key])
        for key in ('internal','applied_load','effective_internal','pore_internal','contact_internal','gravity_load','surcharge_load','exterior_water_load','joint_water_load','holding_load'):
            assert comparison['all_fields'][key]['maximum_absolute']<1e-10,(name,key,comparison['all_fields'][key])
        assert all(v['mismatches']==0 for v in comparison['all_fields'].values() if 'mismatches' in v)
        report['cases'][name]=dict(native=observed,comparison=comparison,preparation=prepared,independent=read(folder/'independent-result.json'),
                                  restart=restart(folder,probe,a.steps),tampering=tamper(folder,probe,rehash,a.steps))
    report['refusals']=refusals(out,solver,a.steps)
    assert fingerprints=={f:sha(ROOT/f) for f in fingerprints}
    report['fixed_sources_unchanged_after_solves']=True
    report['limited_hydrostatic_scope_gates_passed']=True
    report['decisions']={
        'existing_fully_submerged_fresh_MC_cad_new_nodes_matching_level_single_CONTINUOUS_B':'limited acceptance reconfirmed',
        'steady_Darcy_SOIL_BIRTH':'not achieved / current refusal retained',
        'steady_Darcy_CONTINUOUS_BIRTH':'not achieved / current refusal retained',
        'steady_Darcy_PRESCRIBED_BIRTH':'not achieved / current refusal retained',
        'partial_submergence':'not achieved / current refusal retained',
        'fresh_MCC_BIRTH':'not achieved / current refusal retained',
        'aged_material_history_transfer':'not achieved / current refusal retained',
        'undefined_initial_pore_pressure':'not achieved / current refusal retained'}
    write(out/'measurements.json',report)
    print(json.dumps(dict(output=str(out),a03_complete=False,cases=len(report['cases']),refusals=len(report['refusals'])),ensure_ascii=False))
if __name__=='__main__':main()
