"""NEW offset scope: deterministic soil + plate + beam lifecycle and direct-system study."""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import numpy as np
import fem_reference as ref
import plate_surfaces
from verify_plate6_offset import verify
from verify_high_order import prepare,internal,physical
from beam_reference import coupled_forces
from plate6_offset_reference import System
import plate_io
from beam_reference import verify_state

ROOT=Path(__file__).resolve().parents[1]

def inputs(e=0.05):
    values=[json.loads((ROOT/('examples/plate6-offset-'+name+'.json')).read_text()) for name in ('mesh','ground','stages')]
    for plate in values[2]['beam_system']['plates']:plate['attachment_offset_m']=e
    return *values,dict(schema='soil-fem-offset-direct-input/1',attachment='structural positions/translations at soil attachment; neutral plane x_s-e*n',offset_m=e)

def write(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def command(binary,args):
    p=subprocess.run([str(binary),*map(str,args)],capture_output=True,text=True,check=False)
    if p.returncode:raise RuntimeError(f'CLI failed {p.returncode}: {p.stdout}\n{p.stderr}')
    return p

def semantic(run):
    # Runtime is the sole nonphysical field excluded from restart comparison.
    if isinstance(run,dict):return {k:semantic(v) for k,v in run.items() if k!='elapsed_seconds'}
    if isinstance(run,list):return [semantic(v) for v in run]
    return run

def zero_equivalence(output,binary=None):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    binary=Path(binary or ROOT/'target/release'/('soil-fem-3d.exe' if os.name=='nt' else 'soil-fem-3d'))
    mesh,ground,case,_=inputs(0.)
    for name,value in [('mesh',mesh),('ground',ground),('offset-stages',case)]:write(output/(name+'.json'),value)
    base=[output/'mesh.json',output/'ground.json']
    command(binary,['solve-stages',*base,output/'offset-stages.json',output/'offset.json'])
    for plate in case['beam_system']['plates']:plate['formulation']='mindlin_tri6';plate.pop('attachment_offset_m')
    write(output/'planar-stages.json',case)
    command(binary,['solve-stages',*base,output/'planar-stages.json',output/'planar.json'])
    a,b=(json.loads((output/n).read_text()) for n in ('offset.json','planar.json'))
    assert a['complete'] and b['complete']
    keys=('frames','last_committed','attempts','geostatic','transition_audits')
    for key in keys:assert semantic(a[key])==semantic(b[key]),key
    result=dict(passed=True,roundoff_or_better='exact equality',fields=list(keys),excluded_fields=['elapsed_seconds'])
    write(output/'zero-equivalence.json',result);return result

def study(output,binary=None,e=0.05):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    binary=Path(binary or ROOT/'target/release'/('soil-fem-3d.exe' if os.name=='nt' else 'soil-fem-3d'))
    mesh,ground,case,binding=inputs(e)
    (output/'ground').mkdir();(output/'solver').mkdir()
    solver=output/'solver/soil-fem-3d.exe';shutil.copy2(binary,solver)
    for name,value in [('ground/mesh.json',mesh),('ground/ground.json',ground),('stages.json',case),('plate-surface-binding.json',binding)]:write(output/name,value)
    digest=hashlib.sha256(solver.read_bytes()).hexdigest();write(output/'solver-binding.json',dict(solver_sha256=digest))
    base=[output/'ground/mesh.json',output/'ground/ground.json',output/'stages.json']
    command(solver,['solve-stages',*base,output/'stage-results.json'])
    run=json.loads((output/'stage-results.json').read_text());assert run['complete'],run['failure_reason']
    verify(output,output/'independent.json')
    reference=json.loads((output/'independent.json').read_text())
    increment=np.asarray(run['frames'][2]['state']['displacement_m'])-run['frames'][1]['state']['displacement_m']
    component_maxima=np.max(np.abs(increment),axis=0)
    assert np.all(component_maxima>1e-10)
    assert max(p['membrane_energy_kn_m'] for f in run['frames'] for p in f['beam_results'].get('plates',[]))>0.
    # Reassemble accepted saved stresses and structural forces, independently.
    normal_ground=ref.normalize_ground(copy.deepcopy(ground));model=System(mesh,run['case']['beam_system'])
    quadrature,offset=prepare(mesh,run['integration_points']);reassembled=[]
    for frame in run['frames']:
        state=frame['state'];active=np.asarray(state['active_elements'])
        fi,_=internal(mesh,normal_ground,state['material_states'],active,quadrature,offset)
        fe,_=physical(mesh,normal_ground,active,state['surface_pressure_kpa'],case.get('pressure_on_original_surface_only',False))
        fe+=state['boundary_factor']*ref.conditions_load(mesh,case['boundary'],active)
        fi,fe,_=coupled_forces(model,state,fi,fe,frame['beam_results'])
        basis,_=ref.staged_constraint_space(mesh,run['case'],state)
        residual=np.concatenate([b.T@(fi[n]-fe[n]) for n,b in basis.items()])
        relative=float(np.linalg.norm(residual)/max(np.linalg.norm(fe),1.))
        assert relative<1e-10,relative
        reassembled.append(relative)
    attempts=run['geostatic']['attempts']+run['attempts']
    assert all(a['relative_force_residual']<=1e-8 for a in attempts if a['accepted'])
    linear=[l['true_relative_residual'] for a in attempts if a['accepted'] for l in a['linear_solves']]
    assert all(l<=1e-10 for l in linear)
    assert all(a['committed_state_unchanged_on_rejection'] for a in attempts if not a['accepted'])
    assert max(p['transverse_shear_energy_kn_m'] for f in run['frames'] for p in f['beam_results'].get('plates',[]))>0.
    # Verify saved report equations as well as the independent directly coupled solve.
    for frame in run['frames']:verify_state(model,frame['state']['beams'],frame['beam_results'],frame['state']['displacement_m'],frame['state']['boundary_factor'])
    trace_error=0.
    for frame in run['frames']:
        history=frame['state']['beams'];d=np.c_[history['node_displacement_m'],history['node_rotation_rad']]
        for i,(nodes,_,kernel) in enumerate(model.plate_data):
            if not history['active_plates'][i]:continue
            soil=[model.attach[n] for n in nodes];assert all(n is not None for n in soil)
            for l in ([.2,.3,.5],[.1,.8,.1],[0.,.37,.63],[.23,0.,.77],[.61,.39,0.]):
                shape=ref.basis(np.array([[l[1],l[2],0.]]),True)[0][0][[0,1,2,4,5,6]]
                expected=shape@np.asarray(frame['state']['displacement_m'])[soil]
                trace_error=max(trace_error,float(np.max(np.abs(expected-kernel.attachment_translation(l,d[nodes].ravel())))))
    assert trace_error<1e-10
    export=output/'plate-output';export.mkdir();plate_io.export(run,export)
    # Bounded compilation, including a zero-byte cache, must reproduce the default.
    bounded=[]
    for cache in (0,64):
        variant=copy.deepcopy(case);variant['beam_system']['storage_budget']=dict(maximum_nodes=500,maximum_compilation_mib=64,maximum_cached_mib=cache)
        cp=output/f'bounded-{cache}.json';rp=output/f'bounded-{cache}-run.json';write(cp,variant)
        command(solver,['solve-stages',*base[:2],cp,rp]);rr=json.loads(rp.read_text());assert rr['complete']
        for a,b in zip(rr['frames'],run['frames'],strict=True):
            np.testing.assert_allclose(a['state']['displacement_m'],b['state']['displacement_m'],atol=1e-12,rtol=1e-10)
            np.testing.assert_allclose(a['state']['beams']['node_rotation_rad'],b['state']['beams']['node_rotation_rad'],atol=1e-12,rtol=1e-10)
        bounded.append(cache)
    no_reject=copy.deepcopy(case);no_reject['verification_reject_once_stages']=[];no_reject['settings']['initial_increment']=0.125
    for stage in no_reject['stages']:
        first=next(a for a in run['attempts'] if a['stage']==stage['name'] and a['accepted'])
        stage['initial_increment']=first['progress_to']-first['progress_from']
    cp=output/'no-rejection-case.json';rp=output/'no-rejection-run.json';write(cp,no_reject)
    ng=copy.deepcopy(ground);ng['initialization']['geostatic']['verification_reject_once']=False
    ngp=output/'no-rejection-ground.json';write(ngp,ng)
    command(solver,['solve-stages',base[0],ngp,cp,rp]);clean=json.loads(rp.read_text());assert clean['complete']
    assert all(a['accepted'] for a in clean['geostatic']['attempts']+clean['attempts'])
    accepted=lambda seq:[semantic(a) for a in seq if a['accepted']]
    assert accepted(clean['attempts'])==accepted(run['attempts']), 'accepted construction path differs'
    assert accepted(clean['geostatic']['attempts'])==accepted(run['geostatic']['attempts']), 'accepted geostatic path differs'
    assert semantic(clean['frames'])==semantic(run['frames']), 'rejection-free frames differ'
    assert semantic(clean['last_committed'])==semantic(run['last_committed']), 'rejection-free final state differs'
    def accepted_run(value):
        value=semantic(copy.deepcopy(value))
        value['attempts']=[a for a in value['attempts'] if a['accepted']]
        value['geostatic']['attempts']=[a for a in value['geostatic']['attempts'] if a['accepted']]
        value['geostatic'].pop('verification_rejection_done',None)
        value['case'].pop('verification_reject_once_stages',None)
        value['case']['settings'].pop('initial_increment')
        for stage in value['case']['stages']:stage.pop('initial_increment',None)
        value['ground']['initialization']['geostatic'].pop('verification_reject_once',None)
        for key in ('ground_sha256','case_sha256'):value['provenance'].pop(key)
        return value
    assert accepted_run(clean)==accepted_run(run), 'accepted Run differs beyond rejection/step-control inputs and their hashes'
    # The accepted-increment CLI checkpoint option is discovered from the repo's
    # stable CLI contract. Every solve and resume is a separate process.
    resumes=[]
    for stop in (1,9,17,25,29):
        target=output/f'stopped-{stop}.json';checkpoint=output/f'checkpoint-{stop}.json'
        command(solver,['checkpoint-stages',*base,stop,checkpoint])
        sealed=json.loads(checkpoint.read_text())
        result=output/f'resumed-{stop}.json'
        command(solver,['resume-stages',*base,checkpoint,result])
        assert semantic(json.loads(result.read_text()))==semantic(run)
        resumes.append(stop)
    # Corrupt history and change formulation: checksum/input identity must reject.
    sealed=json.loads((output/'checkpoint-17.json').read_text())
    bad=copy.deepcopy(sealed)
    state=bad['run']['last_committed']['beams']
    state['plate_reference_global_kinematics'][0][0]+=.001
    corrupt=output/'corrupt-reference.json';write(corrupt,bad)
    p=subprocess.run([str(solver),'resume-stages',*map(str,base),str(corrupt),str(output/'must-not-exist.json')],capture_output=True)
    assert p.returncode!=0
    variant=copy.deepcopy(case);variant['beam_system']['plates'][0]['formulation']='cst_dkt'
    changed=output/'changed-type.json';write(changed,variant)
    p=subprocess.run([str(solver),'resume-stages',*map(str,base[:2]),str(changed),str(output/'checkpoint-17.json'),str(output/'must-not-exist-type.json')],capture_output=True)
    assert p.returncode!=0
    result=dict(schema='soil-fem-plate6-offset-acceptance/1',passed=True,source_head='cecbb39edf7359f02ca3968998532ebd8612d544',attachment_offset_m=e,
        reassembled_relative_residuals=reassembled,maximum_attachment_trace_error_m=trace_error,maximum_accepted_true_linear_residual=max(linear,default=0.),no_rejection_comparison=True,accepted_run_exact_comparison=True,
        loaded_soil_increment_max_abs_by_component_m=component_maxima.tolist(),frames=len(run['frames']),accepted_increments=sum(a['accepted'] for a in run['attempts']),rejected_attempts=sum(not a['accepted'] for a in attempts),maximum_accepted_force_residual=max(a['relative_force_residual'] for a in attempts if a['accepted']),
        separate_process_restart_stops=resumes,bounded_cache_mib=bounded,reference=reference,quadratic_vtu_roundtrip=True,
        checkpoint_reference_and_type_tampering_rejected=True,p6_complete=False,
        scope='NEW: straight planar offset Tri6, |e|<=t/2, dry elastic soil + plate + beam; no JOINT, water, soil BIRTH/DEATH, folds, curved shell, large model or GUI/GPU acceptance.')
    write(output/'acceptance.json',result)
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('output',type=Path);parser.add_argument('--binary',type=Path);parser.add_argument('--offset',type=float,default=0.05)
    args=parser.parse_args();print(json.dumps(study(args.output,args.binary,args.offset),indent=2))
