"""A certified P1 field must embed in the quadratic proposal's feasible set."""
from pathlib import Path
import hashlib,importlib.util,json,math,os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np
p=Path(__file__).resolve().parent;study=p/'slope1200-durable-SRM-v1';out=study/'quadratic-known-witness-v1';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
source=p/'generate-plane-strain-quadratic-work-v3.py';known=study/'tighter-static-bounds-v1/F1p75/candidate.json'
actual=study/'tighter-static-bounds-v1/F1p75/actual-Rust-all-GP-result.json';meshfile=p/'srm-general-virtual-work-v1/slope1200-LP-v1/mesh.json'
assert read(study/'tighter-static-bounds-v1/F1p75/independent-PASS.json')['passed']
binding={str(q):sha(q) for q in [source,known,actual,meshfile,Path(__file__)]}
spec=importlib.util.spec_from_file_location('quadratic_embedding',source);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
mesh=read(meshfile);problem=read(known);xyz=np.asarray(mesh['nodes']);velocity=np.asarray(problem['virtual_displacement_m']);factor=1/max(1.,float(np.max(np.abs(velocity))))
velocity*=factor;used=set(np.asarray(mesh['tetrahedra10']).reshape(-1).tolist());groups={};group_of={}
for n in sorted(used):
    key=tuple((0. if xyz[n,i]==0 else float(xyz[n,i])).hex() for i in [0,2])
    if key not in groups:groups[key]=len(groups)
    group_of[n]=groups[key]
fixed=set()
for c in problem['constraints']:
    axis=next(i for i,x in enumerate(c['direction']) if x!=0)
    assert all(velocity[n,axis]==0 for n in mesh['node_sets'][c['node_set']] if n in used)
    if axis!=1:fixed.update((group_of[n],axis) for n in mesh['node_sets'][c['node_set']] if n in used)
variables={item:i for i,item in enumerate((g,a) for g in range(len(groups)) for a in [0,2] if (g,a) not in fixed)}
values=np.zeros(len(variables));seen_values={}
for n in sorted(used):
    for axis in [0,2]:
        key=(group_of[n],axis)
        if key in seen_values:assert seen_values[key]==float(velocity[n,axis]),'P1 field is not in the depth-tied quadratic space'
        seen_values[key]=float(velocity[n,axis])
        if key in variables:values[variables[key]]=velocity[n,axis]
assert max(abs(values))<=1.
phi=math.radians(problem['materials'][0]['parameters']['friction_deg']);apex=problem['materials'][0]['parameters']['cohesion_kpa']/math.tan(phi)
angles=np.arange(64)*2*math.pi/64;inflation=(1+1e-4)/math.cos(math.pi/64)
maximum_violation=-math.inf;capacity=[];triangles=set()
for corners,nodes in zip(mesh['tetrahedra'],mesh['tetrahedra10']):
    signature,indices,T,D,S,mapping=m.projected_rows(xyz,np.asarray(nodes),np.asarray(corners),group_of,variables);triangles.add(signature)
    trace=T@values[indices];difference=D@values[indices];shear=S@values[indices]
    lhs=math.sin(phi)*inflation*(difference[:,None]*np.cos(angles)+shear[:,None]*np.sin(angles))-trace[:,None]
    maximum_violation=max(maximum_violation,float(np.max(lhs)));assert np.max(lhs)<1e-10
    x=xyz[corners];J=(x[1:]-x[0]).T;V=float(np.linalg.det(J))/6
    grad=m.derivatives([.25,.25,.25])@np.linalg.inv(J)
    capacity.append(V*apex*float(np.trace((velocity[nodes]-velocity[nodes[0]]).T@grad)))
work=float(np.asarray(problem['physical_nodal_force_kn'])@velocity.reshape(-1));cap=math.fsum(capacity)
certificate=read(actual)['assessment']['certificate']
for got,key in [(work,'external_virtual_work_kn_m'),(cap,'maximum_admissible_internal_virtual_work_kn_m')]:
    expected=factor*certificate[key];assert abs(got-expected)<=2e-12*max(1.,abs(expected)),(got,expected)
assert work-cap>1e-5*max(abs(work),abs(cap))
assert all(sha(Path(q))==h for q,h in binding.items())
proof=dict(passed=True,source_files_sha256=binding,quadratic_variables=len(variables),condensed_triangles=len(triangles),
    all307200_uncondensed_vertex_inequalities_covered=True,
    known_independently_certified_F1p75_P1_field_embeds_in_same_quadratic_LP=True,
    homogeneity_scale_to_respect_unit_variable_box=factor,maximum_inner_polygon_inequality_violation=maximum_violation,
    actual_work_and_capacity_match_existing_certificate_scaled_at_original2e_minus12=True,
    scaled_work=work,scaled_capacity=cap,scaled_excess=work-cap,
    bounded_optimizer_timeout_is_not_absence_of_a_feasible_witness=True,p6_complete=False)
with (out/'terminal.json').open('x',encoding='utf-8') as f:json.dump(proof,f,indent=2,allow_nan=False)
print(json.dumps({k:v for k,v in proof.items() if k!='source_files_sha256'}),flush=True)
