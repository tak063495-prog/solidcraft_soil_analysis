"""Bind actual slope loads and propose a bounded plane-strain work witness."""
from pathlib import Path
import copy,datetime,hashlib,json,math,os,subprocess,sys,time
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import ijson,numpy as np
p=Path(__file__).resolve().parent
study=p/'plane-strain-slope-comparison-v1';case=study/'xy1-ny1'
out=p/'srm-general-virtual-work-v1/slope1200-LP-v1';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def active():
    if (p/'pause-requested').exists():raise SystemExit(2)
active()
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),
    created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
assert read(case/'solve-v1/terminal.json')['passed']
assert read(case/'field-comparison-v1/full-independent-audit-terminal.json')['passed']
original_sha=read(case/'solve-v1/terminal.json')['status']['artifact_sha256']
result=case/'solver-job/result.json';assert sha(result)==original_sha
mesh=read(case/'mesh.json');ground=read(case/'ground.json');stages=read(case/'stages.json')
assert len(mesh['tetrahedra'])==1200 and ground['materials'][0]['unit_weight_kn_m3']==20.
assert len(ground['materials'])==1 and ground['initialization']['water_level_m'] is None
assert ground['initialization']['surface_surcharge_kpa']==0.
assert set(stages['boundary'])=={'constraints'}
assert stages['stages']==[{'kind':'boundary','name':'retain_full_gravity_state','target_factor':0.,'initial_increment':1.}]
reference=read(study/'reference-2d-v2/terminal.json');Fss=reference['fixed_factor_for_future3d_match']
assert reference['passed'] and abs(Fss-1.35)<1e-14
original=ground['materials'][0]['mohr_coulomb']
assert abs(original['cohesion_kpa']*Fss-10.)<1e-12
assert abs(math.degrees(math.atan(math.tan(math.radians(original['friction_deg']))*Fss))-20.)<1e-12
force=[];active_mask=[];pressure=None;factor=None
with result.open('rb') as f:
    for prefix,event,value in ijson.parse(f,use_float=True):
        if prefix=='last_committed.applied_nodal_load_kn.item':force.append(float(value))
        elif prefix=='last_committed.active_elements.item':active_mask.append(value)
        elif prefix=='last_committed.surface_pressure_kpa':pressure=value
        elif prefix=='last_committed.boundary_factor':factor=value
        elif prefix=='last_committed.stage_progress' and event=='number':
            assert value==1.;break
assert len(force)==3*len(mesh['nodes']) and len(active_mask)==1200 and all(active_mask)
assert pressure==0. and factor==0. and all(math.isfinite(x) for x in force)
# Independent exact integrals of straight Tet10 shape functions: vertices -V/20,
# edge nodes V/5. Compare every force component against the real accepted Run.
xyz=np.asarray(mesh['nodes']);calculated=np.zeros(len(force));volume=0.
for cell,high in zip(mesh['tetrahedra'],mesh['tetrahedra10']):
    x=xyz[cell];V=float(np.linalg.det((x[1:]-x[0]).T))/6.;assert V>0.;volume+=V
    for n in high[:4]:calculated[3*n+2]+=V
    for n in high[4:]:calculated[3*n+2]-=4.*V
error=float(np.max(np.abs(calculated-force)))/max(1.,float(np.max(np.abs(force))))
assert error<1e-11 and abs(volume-220.)<1e-10 and abs(sum(force[2::3])+4400.)<1e-9
assert sha(result)==original_sha
dump(out/'physical-load-proof.json',dict(passed=True,elements=1200,nodes=len(xyz),volume_m3=volume,
    actual_result_sha256=original_sha,every_accepted_load_component_matched_independent_Tet10_shape_integrals=True,
    maximum_scaled_force_error=error,total_body_force_z_kn=sum(force[2::3]),
    same_geometry_boundary_and_body_weight_at_every_future_strength_factor=True,
    original_c0_kpa=10.,original_phi0_deg=20.,archived_fixed_Fss=Fss,
    no_previous_numerical_Fs_bracket_adopted_as_a_physical_upper=True,p6_complete=False))
with (out/'mesh.json').open('x',encoding='utf-8') as f:json.dump(mesh,f,allow_nan=False)
dump(out/'source-binding.json',dict(prepared=True,driver_sha256=sha(Path(__file__)),
    generator_sha256=sha(p/'generate-plane-strain-mc-virtual-work.py'),
    input_files_sha256={str(q):sha(q) for q in [case/'mesh.json',case/'ground.json',case/'stages.json',result]},
    original5percent_displacement_acceptance_remains_FAILED=True,
    only_static_virtual_field_proposals_not_new_nonlinear_solutions=True,p6_complete=False))
rows=[]
for strength_factor in [1.5,2.,3.,4.]:
    active()
    folder=out/('F'+str(strength_factor).replace('.','p'));folder.mkdir()
    parameters=copy.deepcopy(original);parameters['cohesion_kpa']=10./strength_factor
    parameters['friction_deg']=math.degrees(math.atan(math.tan(math.radians(20.))/strength_factor))
    problem=dict(schema='soil-fem-mc-virtual-work-bound/1',water_policy='dry',
        materials=[dict(id=1,parameters=parameters)],constraints=stages['boundary']['constraints'],
        physical_nodal_force_kn=force,active_elements=active_mask,virtual_displacement_m=[[0.,0.,0.] for _ in xyz])
    dump(folder/'problem-without-witness.json',problem)
    command=[sys.executable,str(p/'generate-plane-strain-mc-virtual-work.py'),str(out/'mesh.json'),
        str(folder/'problem-without-witness.json'),str(folder/'candidate.json'),'--sides','16','--time-limit','45']
    start=time.monotonic()
    with (folder/'generator.log').open('x',encoding='utf-8') as log:
        child=subprocess.Popen(command,cwd=p,stdout=log,stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        dump(folder/'generator-child.json',dict(pid=child.pid,parent_pid=os.getpid(),executable=sys._base_executable,
            arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code=child.wait()
    row=dict(factor=strength_factor,exit_code=code,elapsed_seconds=time.monotonic()-start,
        proposed=(folder/'candidate.json').exists(),folder=str(folder))
    rows.append(row);dump(folder/'process-terminal.json',row);print(json.dumps(row),flush=True)
    assert code in [0,2]
    if row['proposed']:break
dump(out/'terminal.json',dict(passed=True,preparation_and_bounded_proposal_attempts_completed=True,
    any_candidate_proposed=any(r['proposed'] for r in rows),rows=rows,
    actual_Rust_all_GP_verification_and_independent_full_GP_comparison_required=True,
    no_physical_slope_upper_or_Fs_or_new_FEM_or_convergence_acceptance_from_LP_status=True,
    original_inputs_and_result_unchanged=sha(result)==original_sha,p6_complete=False))
print('Bounded slope LP proposals finished; full Rust GP verification remains',flush=True)
