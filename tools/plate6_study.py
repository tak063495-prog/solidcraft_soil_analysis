"""Small deterministic Plate6/soil/beam acceptance study; no GUI or large jobs."""
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
from verify_beam_coupling import verify
from verify_high_order import prepare,internal,physical
from beam_reference import System,coupled_forces
import plate_io
from solidcraft_structure_fields import fields

ROOT=Path(__file__).resolve().parents[1]

def inputs():
    mesh=json.loads((ROOT/'examples/hardening-column-tet10-mesh.json').read_text())
    mesh['face_sets']={'top':[list(f) for f,owners in ref.topology(mesh).items() if len(owners)==1 and all(mesh['nodes'][i][2]==1. for i in f)]}
    ground=json.loads((ROOT/'examples/ground-plate-beam-column.json').read_text())
    case=json.loads((ROOT/'examples/stages-plate-beam-column.json').read_text())
    # Keep the bottom fixed; retain all three soil translation components above it.
    for constraint in case['boundary']['constraints']:
        if constraint['node_set']=='all':constraint['node_set']='bottom'
    old=case.pop('beam_system');old['nodes']=[];old['plates']=[];old['rotation_constraints']=[];old['nodal_loads']=[];old['members']=[]
    case['beam_system']=old
    for stage in case['stages']:
        if 'members' in stage:stage['members']=['wall','strut']
    owner={r['id']:r['body'] for r in mesh['regions']}
    body=owner[mesh['region_ids'][ref.topology(mesh)[tuple(sorted(mesh['face_sets']['top'][0]))][0]]]
    spec=dict(schema='soil-fem-plate-surfaces/1',sections=[],surfaces=[dict(id='wall',face_set='top',owner_body=body,section_id='wall_section',initially_active=False,formulation='mindlin_tri6',boundary_traction_kpa=[0.2,-0.1,-0.1])])
    mesh,case,binding=plate_surfaces.attach(mesh,case,spec);cfg=case['beam_system']
    corner=next(n['id'] for n in cfg['nodes'] if n['position_m']==[1.,0.,1.])
    cfg['nodes'].append(dict(id='tip',position_m=[2.,0.,1.],match_soil_node=False))
    cfg['members']=[dict(id='strut',nodes=[corner,'tip'],section_id='strut_section',local_y=[0.,1.,0.],initially_active=False,boundary_line_load_kn_m=[0.,0.,-.2])]
    cfg['nodal_loads']=[dict(node_id='tip',boundary_force_kn=[0.,0.,-5.])]
    return mesh,ground,case,binding

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

def study(output,binary=None):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    binary=Path(binary or ROOT/'target/release'/('soil-fem-3d.exe' if os.name=='nt' else 'soil-fem-3d'))
    mesh,ground,case,binding=inputs()
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
    for i in range(len(run['frames'])):fields(mesh,run,i)
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
    cp=output/'no-rejection-case.json';rp=output/'no-rejection-run.json';write(cp,no_reject)
    ng=copy.deepcopy(ground);ng['initialization']['geostatic']['verification_reject_once']=False
    ngp=output/'no-rejection-ground.json';write(ngp,ng)
    command(solver,['solve-stages',base[0],ngp,cp,rp]);clean=json.loads(rp.read_text());assert clean['complete']
    assert all(a['accepted'] for a in clean['geostatic']['attempts']+clean['attempts'])
    for a,b in zip(clean['frames'],run['frames'],strict=True):
        np.testing.assert_allclose(a['state']['displacement_m'],b['state']['displacement_m'],atol=1e-9,rtol=0.)
        for key in ('node_displacement_m','node_rotation_rad','plate_reference_global_kinematics'):
            np.testing.assert_allclose(a['state']['beams'][key],b['state']['beams'][key],atol=1e-9,rtol=0.)
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
    result=dict(schema='soil-fem-plate6-mindlin-acceptance/1',passed=True,source_head='86d2e3c3972b9bc5ae17e77972911df755b58366',
        reassembled_relative_residuals=reassembled,maximum_accepted_true_linear_residual=max(linear,default=0.),no_rejection_comparison=True,
        loaded_soil_increment_max_abs_by_component_m=component_maxima.tolist(),frames=len(run['frames']),accepted_increments=sum(a['accepted'] for a in run['attempts']),rejected_attempts=sum(not a['accepted'] for a in attempts),maximum_accepted_force_residual=max(a['relative_force_residual'] for a in attempts if a['accepted']),
        separate_process_restart_stops=resumes,bounded_cache_mib=bounded,reference=reference,quadratic_vtu_roundtrip=True,
        checkpoint_reference_and_type_tampering_rejected=True,p6_complete=False,
        scope='Straight planar same-reference-surface Tri6, small elastic soil + plate + beam. No offset, JOINT, water, soil birth, curved shell, large model or GUI/GPU acceptance.')
    write(output/'acceptance.json',result)
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('output',type=Path);parser.add_argument('--binary',type=Path)
    args=parser.parse_args();print(json.dumps(study(args.output,args.binary),indent=2))
