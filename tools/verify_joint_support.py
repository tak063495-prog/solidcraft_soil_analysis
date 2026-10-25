"""A zero-pressure contact must not certify a free vertical rigid mode."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import numpy as np
import fem_reference as ref
from joint_reference import System
from verify_checkpoint import invoke
from verify_rigid_mode_theory import verify as mode_theory


def verify(solver,output):
    root=Path(__file__).resolve().parents[1];output.mkdir(parents=True,exist_ok=False)
    theory=mode_theory()
    mesh=json.loads((root/'examples/joint-column-tet10-mesh.json').read_text())
    case=json.loads((root/'examples/stages-joint-explicit-normal-column.json').read_text())
    xyz=np.asarray(mesh['nodes']);mesh['face_sets']['free_top']=[list(face) for face,owners in ref.topology(mesh).items()
        if len(owners)==1 and np.all(abs(xyz[list(face),2]-1)<1e-10)]
    assert mesh['face_sets']['free_top']
    case['boundary']['surface_loads']=[{'kind':'traction','face_set':'free_top','traction_kpa':[0.,0.,110.]}]
    case['stages']=[{'kind':'boundary','name':'unseat','target_factor':1.,'initial_increment':1.}]
    case['verification_reject_once_stages']=[]
    for name,value in (('mesh.json',mesh),('stages.json',case)):
        (output/name).write_text(json.dumps(value,indent=2),encoding='utf-8')
    shutil.copyfile(root/'examples/ground-joint-mc.json',output/'ground.json')
    shutil.copyfile(solver,output/'solver.exe')
    inputs=[str(output/name) for name in ('mesh.json','ground.json','stages.json')]
    # Upward 110 kN balances 100 kN top pressure and 10 kN upper self-weight.
    # With zero tensile capacity, a contact-free vertical translation remains.
    invoke([str(output/'solver.exe'),'solve-stages',*inputs,str(output/'result.json')],output/'solve.log',expected=1)
    data=json.loads((output/'result.json').read_text());assert not data['complete']
    rejected=[a for a in data['attempts'] if not a['accepted']]
    assert rejected and all(a['committed_state_unchanged_on_rejection'] for a in rejected)
    assert any('JOINT rigid mechanism' in a['reason'] for a in rejected)
    assert 0<data['last_committed']['stage_progress']<1
    system=System(mesh,case['joint_system'],data['joint_integration_points'])
    retained=system.rigid_support(case,data['last_committed']);assert retained['unrestrained_modes']==0
    assert len(data['frames'])==1,'unsupported final stage must not have a result frame'
    # Independently replay only the accepted physical boundary factors in a
    # fresh process, omitting every failed endpoint trial. Constitutive work,
    # displacement, load and active histories must agree to roundoff. Splitting
    # one stage into several changes floating-point evaluation order.
    clean=json.loads(json.dumps(case))
    accepted=[a['progress_to'] for a in data['attempts'] if a['accepted']]
    clean['stages']=[{'kind':'boundary','name':f'accepted_{i}','target_factor':factor,'initial_increment':1.}
                     for i,factor in enumerate(accepted)]
    (output/'clean-stages.json').write_text(json.dumps(clean,indent=2),encoding='utf-8')
    invoke([str(output/'solver.exe'),'solve-stages',*inputs[:2],str(output/'clean-stages.json'),str(output/'clean-result.json')],output/'clean.log')
    restored=json.loads((output/'clean-result.json').read_text());assert restored['complete']
    def physical_state(state):return {k:v for k,v in state.items() if k not in ('stage_index','stage_progress')}
    initial=data['frames'][0]['state']
    stress_scale=max(1.,max(abs(v) for s in initial['material_states'] for v in s['effective_stress_compression_positive_kpa']),
                     max(abs(v) for s in initial['joints']['points'] for v in s['traction_kpa']))
    maximum_error=0.
    def compare(a,b,field=''):
        nonlocal maximum_error
        if isinstance(a,dict):
            assert a.keys()==b.keys()
            for key in a:compare(a[key],b[key],key)
        elif isinstance(a,list):
            assert len(a)==len(b)
            for x,y in zip(a,b):compare(x,y,field)
        elif isinstance(a,float):
            # Unloading pressure is a cancellation of the initial stress.
            # Scale stress components by that physical load, not final zero.
            scale=stress_scale if field.endswith('_kpa') else 1.
            error=abs(a-b)/max(scale,abs(a),abs(b));maximum_error=max(maximum_error,error)
            assert error<=128*np.finfo(float).eps,(a,b,error)
        else:assert a==b
    compare(physical_state(restored['last_committed']),physical_state(data['last_committed']))
    report={'schema':'soil-fem-joint-unseating-verification/1','passed':True,
        'rigid_mode_theory':theory,
        'zero_pressure_vertical_mechanism_rejected':True,'all_rejected_trials_preserved_committed_state':True,
        'retained_state_agrees_with_clean_to_roundoff':True,'clean_maximum_normalized_component_error':maximum_error,
        'clean_stress_cancellation_scale_kpa':stress_scale,
        'last_committed_progress':data['last_committed']['stage_progress'],'retained_rigid_support':retained,
        'solver_sha256':hashlib.sha256((output/'solver.exe').read_bytes()).hexdigest(),
        'sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in output.iterdir() if p.is_file()}}
    with (output/'verification.json').open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2)
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();verify(a.solver.resolve(),a.output_dir.resolve())
