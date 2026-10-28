"""Companion independent volume norms; original checker remains authoritative."""
from pathlib import Path
import hashlib,json,os,sys
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
import numpy as np
p=Path(__file__).resolve().parent;root=p.parents[1];sys.path.insert(0,str(root/'tools'))
import fem_reference as ref
from verify_mcc_structural_study import verify_reference_bindings
read=lambda q:json.loads(q.read_text(encoding='utf-8-sig'));sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
gate=read(p/'independent-verification.json');assert gate['independent_all_history_and_member_comparison']=='passed'
assert all(sha(Path(q))==h for q,h in gate['input_sha256'].items())
native=[read(p/str(n)/'analysis/stage-results.json') for n in [16,32]];references=[read(p/'reference-results'/('reference'+str(n)+'.json')) for n in [100,200]]
for r in references:verify_reference_bindings(r)
paths=[Path(__file__),Path(ref.__file__),root/'tools/verify_mcc_structural_study.py',p/'independent-verification.json',*[Path(q) for q in gate['input_sha256']]]
for r in references:paths.extend(Path(q) for q in r['input_and_source_sha256'])
binding={str(q.resolve()):sha(q) for q in paths}
mesh=read(p/'16/analysis/ground/mesh.json');cells=ref.cells(mesh);geometry=ref.Geometry(mesh)
assert len(cells)==274 and native[0]['integration_points']==native[1]['integration_points']
quad=[geometry.volume(e,independent=True,independent_order=5) for e in range(len(cells))]
moments=np.zeros(4)
for e,q in enumerate(quad):
 assert len(q['weight'])==125 and np.all(q['weight']>0)
 xyz=q['shape']@np.asarray(mesh['nodes'])[cells[e]]
 moments+=np.r_[q['weight'].sum(),q['weight']@(xyz**2)]
expected=np.asarray([240.,8000.,1280.,2880.]);assert np.max(abs(moments-expected)/expected)<1e-12
owners=np.asarray([q['element_index'] for q in native[0]['integration_points']]);weights=np.asarray([q['weight_m3'] for q in native[0]['integration_points']])
assert len(owners)==7398 and np.all(weights>0) and np.all(np.bincount(owners)==27)
assert abs(weights.sum()-240)<1e-9
def norm(a,b,w,floor=1e-18):
 a,b=np.asarray(a),np.asarray(b);w=np.asarray(w);assert a.shape==b.shape and len(a)==len(w)
 d=(a-b).reshape(len(a),-1);v=b.reshape(len(b),-1)
 dn=float(np.sqrt(np.sum(w[:,None]*d*d)));vn=float(np.sqrt(np.sum(w[:,None]*v*v)))
 return dict(relative_l2=dn/max(vn,floor),l2_difference=dn,reference_l2=vn,rms_difference=dn/float(np.sqrt(w.sum())),maximum_absolute_difference=float(abs(d).max()))
def displacement(u,active):
 x=[];w=[]
 for e,q in enumerate(quad):
  if active[e]:x.extend(q['shape']@u[cells[e]]);w.extend(q['weight'])
 return np.asarray(x),np.asarray(w)
def state(frame,is_native):
 if is_native:
  f=frame['state'];points=f['material_states']
  material=np.asarray([[*v['effective_stress_compression_positive_kpa'],np.log(v['preconsolidation_pressure_kpa']),np.log1p(v['void_ratio']),*v['plastic_strain_compression_positive']] for v in points])
  return np.asarray(f['displacement_m']),np.asarray(f['active_elements']),material
 return np.asarray(frame['displacement_m']),np.asarray(frame['active_elements']),np.asarray(frame['material_state']).reshape(-1,14)
def fields(material):
 return dict(stress_frobenius_kpa=material[:,:6]*np.sqrt([1,1,1,2,2,2]),pc_kpa=np.exp(material[:,6]),void_ratio=np.expm1(material[:,7]),plastic_strain_frobenius=material[:,8:]*np.sqrt([1,1,1,.5,.5,.5]))
histories=native+references;previous=[np.zeros((len(mesh['nodes']),3)) for _ in histories];rows=[]
for index,frames in enumerate(zip(*[r['frames'] for r in histories],strict=True)):
 assert len({f['name'] for f in frames})==1
 states=[state(f,i<2) for i,f in enumerate(frames)];active=states[0][1]
 assert all(np.array_equal(v[1],active) for v in states)
 materials=[fields(s[2]) for s in states];keep=active[owners];all_gp=np.ones(len(owners),bool)
 evaluations=[displacement(s[0],active) for s in states];increments=[displacement(s[0]-prior,active)[0] for s,prior in zip(states,previous)]
 w=evaluations[0][1];groups={}
 for label,i,j in [('reference100_200',2,3),('coarse16_reference200',0,3),('fine32_reference200',1,3)]:
  group=dict(displacement=norm(evaluations[i][0],evaluations[j][0],w),stage_increment=norm(increments[i],increments[j],w))
  for scope,mask in [('active',keep),('all_history',all_gp)]:
   for name in materials[i]:group[scope+'_'+name]=norm(materials[i][name][mask],materials[j][name][mask],weights[mask])
  groups[label]=group
 rows.append(dict(stage=frames[0]['name'],active_volume_m3=float(w.sum()),comparisons=groups));previous=[s[0] for s in states]
failures=[]
for row in rows:
 for label,group in row['comparisons'].items():
  tolerance=1e-4 if label=='reference100_200' else .01
  for name,value in group.items():
   if value['relative_l2']>tolerance:failures.append(dict(stage=row['stage'],comparison=label,field=name,value=value['relative_l2'],tolerance=tolerance))
maxima={tag:{k:max(row['comparisons'][tag][k]['relative_l2'] for row in rows) for k in rows[0]['comparisons'][tag]} for tag in rows[0]['comparisons']}
for k,coarse in maxima['coarse16_reference200'].items():
 fine=maxima['fine32_reference200'][k]
 if coarse>1e-8 and fine>coarse/2:failures.append(dict(field=k,coarse=coarse,fine=fine,cause='maximum error did not halve'))
assert all(sha(Path(q))==h for q,h in binding.items())
report=dict(passed=not failures,all5_stages=True,reference_refinement_gate=1e-4,candidate_gate=.01,error_halving_required=True,displacement_norm='Positive tensor-Legendre5 Duffy quadrature; squared quadratic Tet10 fields integrated exactly on each affine cell',displacement_quadrature_points=34250,stress_and_plastic_norms='Physical3D Frobenius tensors with original positive material-point volume weights',analytic_domain_moments=dict(computed=moments.tolist(),expected=expected.tolist(),relative_error_max=float(np.max(abs(moments-expected)/expected))),active_and_removed_material_histories_compared_separately=True,frames=rows,maxima=maxima,failures=failures,sha256=binding,spatial_convergence_accepted=False,p6_complete=False)
with (p/'volume-fields-verification.json').open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2,allow_nan=False)
print(json.dumps(dict(passed=not failures,failures=failures,maxima=maxima)),flush=True)
assert not failures
