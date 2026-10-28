"""Supplement the unmodified acceptance checker with separate force/moment units."""
from pathlib import Path
import hashlib,json,os,sys
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
import numpy as np
p=Path(__file__).resolve().parent;root=p.parents[1];sys.path.insert(0,str(root/'tools'))
import fem_reference as ref
import beam_reference
from verify_mcc_structural_study import verify_reference_bindings,physical_case
read=lambda q:json.loads(q.read_text(encoding='utf-8-sig'));sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
accepted=read(p/'independent-verification.json');assert accepted['independent_all_history_and_member_comparison']=='passed'
assert all(sha(Path(q))==h for q,h in accepted['input_sha256'].items())
paths=[Path(__file__),p/'independent-verification.json',root/'tools/beam_reference.py',root/'tools/plate_reference.py',root/'tools/fem_reference.py',root/'tools/verify_mcc_structural_study.py']
native=[read(p/str(n)/'analysis/stage-results.json') for n in [16,32]]
references=[read(p/'reference-results'/('reference'+str(n)+'.json')) for n in [100,200]]
for r in references:verify_reference_bindings(r)
paths.extend(Path(q) for q in accepted['input_sha256'])
for r in references:paths.extend(Path(q) for q in r['input_and_source_sha256'])
binding={str(q.resolve()):sha(q) for q in paths}
case=read(p/'16/analysis/stages.json');assert physical_case(json.loads(json.dumps(case)))==physical_case(read(p/'32/analysis/stages.json'))
cfg=ref.normalize_stage_case(case);mesh=read(p/'16/analysis/ground/mesh.json');system=beam_reference.System(mesh,cfg['beam_system'])
assert not system.members and len(system.plates)==128
def relative(a,b,weights=None,floor=1e-12):
 a,b=np.asarray(a,dtype=float),np.asarray(b,dtype=float);assert a.shape==b.shape
 if not a.size:return dict(relative=0.,absolute_l2=0.,reference_l2=0.)
 if weights is None:weights=np.ones(len(a))
 weights=np.asarray(weights);assert len(weights)==len(a) and np.all(weights>0)
 d=(a-b).reshape(len(a),-1);v=b.reshape(len(b),-1)
 dn=float(np.sqrt(np.sum(weights[:,None]*d*d)));vn=float(np.sqrt(np.sum(weights[:,None]*v*v)))
 return dict(relative=dn/max(vn,floor),absolute_l2=dn,reference_l2=vn,normalization_floor=floor)
def kinematics(a,b):
 assert a['active_members']==b['active_members'] and a['active_plates']==b['active_plates']
 result={k:relative(a[k],b[k]) for k in ['node_displacement_m','node_rotation_rad']}
 for k in ['reference_global_kinematics','plate_reference_global_kinematics']:
  aa,bb=np.asarray(a[k]).reshape(-1,6),np.asarray(b[k]).reshape(-1,6)
  result[k+'_translation_m']=relative(aa[:,:3],bb[:,:3]);result[k+'_rotation_rad']=relative(aa[:,3:],bb[:,3:])
 return result
def sections(a,b):
 assert set(a)==set(b);result={}
 for key in ['membrane_force_kn_m','bending_moment_kn']:
  x=[];y=[];weights=[]
  for identifier in sorted(a):
   left,right=a[identifier],b[identifier]
   assert abs(left['area_m2']-right['area_m2'])<1e-12
   assert np.max(abs(np.asarray(left['local_axes'])-np.asarray(right['local_axes'])))<1e-12
   assert len(left['samples'])==len(right['samples'])==3
   for la,lb in zip(left['samples'],right['samples'],strict=True):
    assert np.max(abs(np.asarray(la['barycentric'])-np.asarray(lb['barycentric'])))<1e-14
    x.append(np.asarray(la[key])*np.sqrt([1,1,2]));y.append(np.asarray(lb[key])*np.sqrt([1,1,2]));weights.append(left['area_m2']/3)
  result[key]=relative(x,y,weights,floor=1e-6)
 return result
rows=[]
for index,frames in enumerate(zip(*[r['frames'] for r in native+references],strict=True)):
 c,f,r100,r200=frames;assert len({v['name'] for v in frames})==1
 factor=c['state']['boundary_factor'];assert factor==f['state']['boundary_factor']
 response=[system.solve_interface(r['beams'],np.asarray(r['displacement_m']),factor).get('plates',{}) for r in [r100,r200]]
 actual=[{v['id']:v for v in r['beam_results'].get('plates',[])} for r in [c,f]]
 for r,v in zip([c,f],actual):assert len(v)==sum(r['state']['beams']['active_plates'])
 refined=dict(kinematics=kinematics(r100['beams'],r200['beams']),sections=sections(*response))
 candidates={tag:dict(kinematics=kinematics(r['state']['beams'],r200['beams']),sections=sections(v,response[1])) for tag,r,v in zip(['coarse16','fine32'],[c,f],actual)}
 for group in refined.values():assert all(value['relative']<=1e-4 for value in group.values()),(c['name'],'reference member refinement',refined)
 for tag,groups in candidates.items():
  for group in groups.values():assert all(value['relative']<=.01 for value in group.values()),(c['name'],tag,groups)
 rows.append(dict(stage=c['name'],reference100_200=refined,candidates=candidates))
assert all(sha(Path(q))==h for q,h in binding.items())
output=dict(passed=True,original_acceptance_checker_preserved=True,independent_member_reference_refinement_gate=1e-4,candidate_member_gate=.01,force_and_moment_units_separate=True,translation_and_rotation_units_separate=True,reference_offsets_units_separate=True,section_norm='Area-weighted Frobenius tensor norm; three degree2-exact triangular quadrature samples on each unchanged linear plate',section_force_source='Read-only independent3D plate recovery from independently propagated structural and soil states',frames=rows,sha256=binding,spatial_convergence_accepted=False,p6_complete=False)
with (p/'member-units-verification.json').open('x',encoding='utf-8') as stream:json.dump(output,stream,indent=2,allow_nan=False)
print(json.dumps(dict(passed=True,stages=len(rows),max_reference_relative=max(v['relative'] for row in rows for group in row['reference100_200'].values() for v in group.values()),max_candidate_relative=max(v['relative'] for row in rows for groups in row['candidates'].values() for group in groups.values() for v in group.values()))),flush=True)
