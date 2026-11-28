"""Reproducible small extension-(a) evidence; failures are never accepted."""
import argparse
import copy
import hashlib
import json
import os
import shutil
from pathlib import Path
import subprocess
import numpy as np
from joint_ext_reference import audit_tip_elastic, audit_reassembly, audit_wet_birth_elastic, audit_wet_birth_holding, audit_clean_path

ROOT=Path(__file__).resolve().parents[1]
def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def dump(path,value):Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def without_clock(v):
    if isinstance(v,dict):return {k:without_clock(x) for k,x in v.items() if k!='elapsed_seconds'}
    if isinstance(v,list):return [without_clock(x) for x in v]
    return v
def exact(v):return json.dumps(without_clock(v),sort_keys=True,separators=(',',':'),allow_nan=False)

def tip_fixture(policy='transfer_force',lifecycle=False):
    mesh=read(ROOT/'examples/joint-crack-tip-mesh.json')
    ground=read(ROOT/'examples/ground-joint-crack-tip.json')
    case=read(ROOT/'examples/stages-joint-crack-tip.json')
    case['joint_system']['joints'][0]['birth_policy']=policy
    ground['initialization']['geostatic']['initial_boundary_factor']=-.001
    case['stages']=[dict(kind='joint_birth',name='contact',joints=['Seam']),dict(kind='pressure',name='reload',target_pressure_kpa=20.),dict(kind='joint_death',name='death',joints=['Seam']),dict(kind='joint_birth',name='rebirth',joints=['Seam'])]
    case['verification_reject_once_stages']=[]
    if lifecycle:
        # Supported top displacements provide a deterministic slip path before
        # removal; the elastic comparison uses the separate force-driven case.
        case=read(ROOT/'examples/stages-joint-crack-tip.json')
        ground=read(ROOT/'examples/ground-joint-crack-tip.json')
        case['joint_system']['joints'][0]['birth_policy']=policy
        case['boundary']['constraints'][-1]['value_m']=0.
        case['stages']=[dict(kind='joint_birth',name='contact',joints=['Seam']),dict(kind='boundary',name='slide',target_factor=1.),dict(kind='boundary',name='unload',target_factor=0.),dict(kind='joint_death',name='death',joints=['Seam']),dict(kind='joint_birth',name='rebirth',joints=['Seam'])]
        case['verification_reject_once_stages']=['contact','slide','unload','death','rebirth']
    return mesh,ground,case

def wet_fixture():
    mesh=read(ROOT/'examples/joint-water-column-tet10-mesh.json'); ground=read(ROOT/'examples/ground-joint-water-mc.json')
    ground['initialization']['water_level_m']=1.5
    case=read(ROOT/'examples/stages-joint-water-column.json')
    case['joint_system']['joints'][0]['birth_policy']='stress_free'
    case['settings']['relative_force_tolerance']=1e-11
    case['settings']['kinematics']='explicit'
    case['boundary']=read(ROOT/'examples/stages-joint-column.json')['boundary']
    case['initially_inactive_bodies']=['UpperClay']
    seed=dict(effective_stress_compression_positive_kpa=[5.,5.,5.,0.,0.,0.],plastic_strain_compression_positive=[0.]*6,void_ratio=.7,equivalent_plastic_strain=0.,plastic_work_kpa=0.,last_return='elastic')
    birth=dict(kind='birth',name='fill',bodies=['UpperClay'],states={'UpperClay':seed},placement='cad_new_nodes',hydraulic_initialization=dict(kind='hydrostatic',level_m=1.5))
    case['stages']=[birth,dict(kind='joint_birth',name='contact',joints=['Seam']),dict(kind='pressure',name='reload',target_pressure_kpa=20.),dict(kind='boundary',name='slide',target_factor=-.1),dict(kind='joint_death',name='death',joints=['Seam']),dict(kind='joint_birth',name='rebirth',joints=['Seam'])]
    again=copy.deepcopy(birth);again['name']='fill_again'
    case['stages'] += [dict(kind='excavate',name='remove_soil',bodies=['UpperClay'],remaining_surface_pressure_kpa=0.),again,dict(kind='joint_birth',name='contact_again',joints=['Seam'])]
    case['verification_reject_once_stages']=['fill','contact','slide','death','rebirth','remove_soil','fill_again','contact_again']
    return mesh,ground,case

def invoke(solver,args,log,success=True):
    result=subprocess.run([str(solver),*map(str,args)],capture_output=True,text=True,encoding='utf-8')
    Path(log).write_text(result.stdout+result.stderr+f'\nexit_code={result.returncode}\n',encoding='utf-8')
    assert (result.returncode==0)==success,(str(log),result.returncode,result.stderr)
    return result

def solve(solver,folder,mesh,ground,case):
    folder.mkdir(parents=True,exist_ok=False)
    for name,data in [('mesh',mesh),('ground',ground),('case',case)]:dump(folder/f'{name}.json',data)
    inputs=[folder/f'{n}.json' for n in ('mesh','ground','case')]
    invoke(solver,['solve-stages',*inputs,folder/'run.json'],folder/'solve.log')
    run=read(folder/'run.json');assert run['complete'],run.get('failure_reason')
    return run,inputs

def lifecycle(solver,folder,mesh,ground,case,run,inputs,rehash=None):
    for audit in run['transition_audits']:
        if next(s['kind'] for s in case['stages'] if s['name']==audit['stage'])=='joint_birth':
            for p in audit['prepared_state']['joints']['points']:
                for key in ('plastic_slip_m','accumulated_plastic_slip_m','plastic_dissipation_kpa_m','numerical_dissipation_kpa_m','endpoint_work_kpa_m'):
                    assert np.all(np.asarray(p[key])==0.),key
    death=next(i for i,s in enumerate(case['stages']) if s['kind']=='joint_death')
    old=run['frames'][death]['state']['joints']['points']
    final=run['last_committed']['joints']['lifecycle'][0]
    assert exact(final['archived'][0]['points'])==exact(list(enumerate(old)))
    assert any(p['accumulated_plastic_slip_m']>0. and p['plastic_dissipation_kpa_m']>0. for p in old)
    rejected=[a for a in run['attempts'] if not a['accepted']]
    assert rejected and all(a['committed_state_unchanged_on_rejection'] for a in rejected)
    assert all(a['reason']=='verification_injected_rejection_after_convergence' for a in rejected)
    clean=copy.deepcopy(case);clean['verification_reject_once_stages']=[]
    for s in clean['stages']:
        if s['name'] in case['verification_reject_once_stages']:s['initial_increment']=s.get('initial_increment',case['settings']['initial_increment'])/2
    clean_run,_=solve(solver,folder/'clean',mesh,ground,clean)
    for key in ('last_committed','frames','transition_audits','flow_audits','reference_displacement_m','joint_integration_points','integration_points'):
        assert exact(run.get(key))==exact(clean_run.get(key)),key
    accepted=lambda r:[a for a in r['attempts'] if a['accepted']]
    assert exact(accepted(run))==exact(accepted(clean_run))
    complete_comparison=audit_clean_path(run,clean_run)
    count=sum(a['accepted'] for a in run['geostatic']['attempts']); locations=[count];seen=set()
    for a in accepted(run):
        count+=1
        if a['stage'] not in seen or a['progress_to']==1.:
            locations.append(count);seen.add(a['stage'])
    checks=[]
    for count in sorted(set(locations)):
        cp=folder/f'checkpoint-{count}.json';rp=folder/f'resumed-{count}.json'
        invoke(solver,['checkpoint-stages',*inputs,count,cp],folder/f'checkpoint-{count}.log')
        saved=read(cp);state=saved['run']['last_committed'];si=state['stage_index']-1
        if 0<=si<len(case['stages']) and case['stages'][si]['kind']=='joint_birth' and case['joint_system']['joints'][0]['birth_policy']=='transfer_force':
            assert exact(state['applied_nodal_load_kn'])==exact(run['frames'][si]['state']['applied_nodal_load_kn'])
            assert all(v==0. for v in state['transition_release_load_kn'])
            if saved['run'].get('continuation'):assert all(v==0. for v in saved['run']['continuation']['holding_load'])
        invoke(solver,['resume-stages',*inputs,cp,rp],folder/f'resume-{count}.log')
        assert exact(read(rp))==exact(run)
        checks.append(dict(accepted_increments=count,exact_except_elapsed_seconds=True))
    tampering=[]
    if rehash is not None:
        original=read(cp)
        for mode in range(6):
            bad=copy.deepcopy(original);state=bad['run']['last_committed'];joint=state['joints']
            if mode==0:joint['points'][0]['birth']['initial_traction_kpa'][2]-=1.
            elif mode==1:joint['points'][0]['birth']['jump_m'][0]+=.001
            elif mode==2:joint['lifecycle'][0]['generation']+=1
            elif mode==3:
                archived=joint['lifecycle'][0]['archived'][0]['points'][0][1]
                # Preserve its local energy identity: authentication must also
                # reject a self-consistent forgery of a frozen generation.
                archived['plastic_dissipation_kpa_m']+=1.;archived['endpoint_work_kpa_m']+=1.
            elif mode==4:
                if state.get('joint_water'):state['joint_water']['physical_pressure_kpa'][0]+=1.
                else:bad['run']['joint_integration_points'][0]['position_m'][0]+=.01
            else:bad['run']['transition_audits'][0]['prepared_state']['material_states'][0]['effective_stress_compression_positive_kpa'][0]+=1.
            changed=folder/f'tampered-{mode}.json';sealed=folder/f'rehashed-{mode}.json';dump(changed,bad)
            invoke(rehash,[changed,sealed],folder/f'rehash-{mode}.log')
            assert read(sealed)['run_sha256']!=original['run_sha256']
            result=invoke(solver,['resume-stages',*inputs,sealed,folder/f'forbidden-{mode}.json'],folder/f'tamper-{mode}.log',False)
            tampering.append(dict(mutation=mode,outer_hash_recomputed=True,rejected=True,reason=result.stderr.strip()))
    return dict(passed=True,restarts=checks,injected_rejections=len(rejected),same_accepted_attempts_exact=True,states_frames_audits_exact=True,complete_accepted_run=complete_comparison,rehashed_tampering=tampering)

def scope_rejections(solver,out):
    out.mkdir();mesh,ground,case=wet_fixture();checks=[]
    for name in ('missing','mismatched_level','aged_history','active_body','placement','mode_a','steady_continuous','soil_initialization','prescribed_initialization'):
        folder=out/name;folder.mkdir();bad=copy.deepcopy(case)
        if name=='missing':bad['stages'][0].pop('hydraulic_initialization')
        elif name=='mismatched_level':bad['stages'][0]['hydraulic_initialization']['level_m']=1.4
        elif name=='aged_history':bad['stages'][0]['states']['UpperClay']['equivalent_plastic_strain']=.1
        elif name=='active_body':bad['initially_inactive_bodies']=[]
        elif name=='placement':bad['stages'][0]['placement']='current_displacement'
        elif name=='mode_a':bad['joint_system']['joints'][0]['birth_policy']='transfer_force'
        elif name=='steady_continuous':
            flow=read(ROOT/'examples/stages-joint-water-flow.json')['stages'][0]['flow'];flow['heads']=[dict(face_set='bottom',head_m=1.5)]
            bad['stages'].insert(0,dict(kind='steady_seepage',name='flow',flow=flow,exterior_water_level_m=1.5))
        else:bad['stages'][0]['hydraulic_initialization']['kind']=name.removesuffix('_initialization')
        for n,d in [('mesh',mesh),('ground',ground),('case',bad)]:dump(folder/f'{n}.json',d)
        result=invoke(solver,['solve-stages',*[folder/f'{n}.json' for n in ('mesh','ground','case')],folder/'forbidden.json'],folder/'refused.log',False)
        assert not (folder/'forbidden.json').exists()
        checks.append(dict(name=name,rejected=True,reason=result.stderr.strip()))
    richards=[]
    for name,head in [('infiltration',2.),('drainage',-1.),('two_way',0.)]:
        folder=out/f'richards-{name}';folder.mkdir();input=read(ROOT/'examples/one-way-mc.json')
        mechanics=copy.deepcopy(case);mechanics['stages']=[];mechanics['initially_inactive_bodies']=[];mechanics['verification_reject_once_stages']=[]
        input['mechanics']=mechanics;input['flow']['initial_head']['head_m']=head
        if name=='two_way':input['coupling']='two_way'
        for n,d in [('mesh',mesh),('ground',read(ROOT/'examples/ground-one-way-mc.json')),('case',input)]:dump(folder/f'{n}.json',d)
        result=invoke(solver,['solve-one-way',*[folder/f'{n}.json' for n in ('mesh','ground','case')],folder/'forbidden.json'],folder/'refused.log',False)
        assert not (folder/'forbidden.json').exists()
        richards.append(dict(name=name,coupling='two_way' if name=='two_way' else 'ONE_WAY',accepted=False,rejected=True,reason=result.stderr.strip(),coupled_mass_balance_checked=False,time_refinement_checked=False))
    return dict(passed=True,wet_birth=checks,richards=richards)

def main():
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--rehash',type=Path)
    a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    (out/'solver').mkdir();solver=out/'solver'/a.solver.name;shutil.copy2(a.solver.resolve(),solver)
    if a.rehash:
        helper=out/'solver'/a.rehash.name;shutil.copy2(a.rehash.resolve(),helper);a.rehash=helper
    report=dict(schema='soil-fem-joint-ext-combos/1',passed=False,requested_all_combinations_passed=False,solver_sha256=hashlib.sha256(solver.read_bytes()).hexdigest(),cases={},remaining=['steady SOIL/CONTINUOUS/PRESCRIBED wet soil BIRTH','Richards wet JOINT ONE_WAY/two-way','tip body BIRTH/partially active interface rebirth','P6 low-preload NC full removal and wet sliding full-history convergence'])
    try:
        for policy in ('transfer_force','stress_free'):
            mesh,ground,case=tip_fixture(policy);run,_=solve(solver,out/f'tip-elastic-{policy}',mesh,ground,case)
            report['cases'][f'tip-elastic-{policy}']=dict(elastic=audit_tip_elastic(mesh,ground,run),reassembly=audit_reassembly(mesh,ground,run))
            mesh,ground,case=tip_fixture(policy,True);folder=out/f'tip-history-{policy}';run,inputs=solve(solver,folder,mesh,ground,case)
            assert any(p['accumulated_plastic_slip_m']>0. for p in run['frames'][2]['state']['joints']['points'])
            report['cases'][f'tip-history-{policy}']=dict(reassembly=audit_reassembly(mesh,ground,run),lifecycle=lifecycle(solver,folder,mesh,ground,case,run,inputs,a.rehash))
        mesh,ground,case=wet_fixture();case['stages']=case['stages'][:3];case['verification_reject_once_stages']=[]
        case['settings']['kinematics']='side_rollers';case.pop('boundary',None)
        run,_=solve(solver,out/'wet-elastic',mesh,ground,case)
        report['cases']['wet-elastic']=dict(elastic=audit_wet_birth_elastic(mesh,ground,run),reassembly=audit_reassembly(mesh,ground,run),holding=audit_wet_birth_holding(mesh,ground,run))
        mesh,ground,case=wet_fixture();folder=out/'wet-hydrostatic';run,inputs=solve(solver,folder,mesh,ground,case)
        report['cases']['wet-hydrostatic']=dict(reassembly=audit_reassembly(mesh,ground,run),holding=audit_wet_birth_holding(mesh,ground,run),lifecycle=lifecycle(solver,folder,mesh,ground,case,run,inputs,a.rehash))
        report['scope_rejections']=scope_rejections(solver,out/'scope-rejections')
        report['passed']=True
    except Exception as error:
        report['failure']=repr(error);raise
    finally:dump(out/'measurements.json',report)
    print(json.dumps(report,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
