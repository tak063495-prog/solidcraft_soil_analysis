"""Exact overlap of quadratic soil fields with the spatial wall discretization fixed."""
from pathlib import Path
import datetime,hashlib,json,os,subprocess,sys,time
for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
import numpy as np
from scipy.spatial import cKDTree
p=Path(__file__).resolve().parent;root=p.parents[1]
left=p/'mesh5144'
right=p/'mesh8061'
read=lambda q:json.loads(q.read_text(encoding='utf-8-sig'));sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
assert not (p/'comparison.json').exists()
for case in [left,right]:
 while not (case/'study-terminal.json').exists():
  if (p/'pause-requested').exists():raise SystemExit(2)
  time.sleep(5)
sys.path.insert(0,str(left/'runtime/tools'))
from tet10_overlap import Overlay
for name,case in [('left',left),('right',right)]:
 output=p/(name+'-preflight.json')
 with (p/(name+'-preflight.log')).open('x',encoding='utf-8') as log:
  command=[sys.executable,str(p/'preflight.py'),str(case),str(output)]
  child=subprocess.Popen(command,cwd=root,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
  with (p/(name+'-child.json')).open('x') as stream:json.dump(dict(pid=child.pid,parent_pid=os.getpid(),executable=sys.executable,arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),stream)
  assert child.wait()==0
 assert read(output)['passed']
proof_path=Path(read(right/'preparation.json')['original_preparation'])/'wall-spatial-topology-proof.json';proof=read(proof_path)
assert proof['wall_mesh_equivalent_within_tolerance'] and proof['structural_node_position_maximum_matching_distance_m']==0
assert all(sha(Path(q))==h for q,h in proof['sha256'].items())
assert sha(left/'runtime/tools/tet10_overlap.py')==sha(right/'runtime/tools/tet10_overlap.py')
paths=[Path(__file__),p/'preflight.py',proof_path,p/'left-preflight.json',p/'right-preflight.json',left/'runtime/tools/tet10_overlap.py',left/'runtime/tools/fem_reference.py']
for name in ['left','right']:paths.extend(Path(q) for q in read(p/(name+'-preflight.json'))['sha256'])
binding={str(q.resolve()):sha(q) for q in paths}
analyses=[c/'workflow/analysis' for c in [left,right]]
meshes=[read(a/'ground/mesh.json') for a in analyses];results=[read(a/'stage-results.json') for a in analyses]
assert [len(m['tetrahedra']) for m in meshes]==[5144,8061]
assert read(analyses[0]/'ground/ground.json')==read(analyses[1]/'ground/ground.json')
assert read(Path(read(left/'preparation.json')['original_preparation'])/'workflow/conformal/stages.json')==read(Path(read(right/'preparation.json')['original_preparation'])/'workflow/conformal/stages.json')
assert meshes[0]['regions']==meshes[1]['regions']
assert len({r['provenance']['solver_binary_sha256'] for r in results})==1
configs=[r['case']['beam_system'] for r in results]
xyz=[np.asarray([n['position_m'] for n in cfg['nodes']]) for cfg in configs]
distance,node_map=cKDTree(xyz[1]).query(xyz[0]);assert distance.max()==0 and len(set(node_map))==len(xyz[0])==440
def plate_keys(cfg):
 nodes={n['id']:tuple(n['position_m']) for n in cfg['nodes']}
 answer={}
 for plate in cfg['plates']:
  corners=[nodes[n] for n in plate['nodes']]
  key=min(tuple(corners[k:]+corners[:k]) for k in range(3))
  assert key not in answer;answer[key]=plate
 return answer
keyed=[plate_keys(c) for c in configs];assert set(keyed[0])==set(keyed[1]) and len(keyed[0])==832
def metric(a,b,w):
 a,b=np.asarray(a),np.asarray(b);w=np.asarray(w);assert a.shape==b.shape and len(w)==len(a)
 d=(a-b).reshape(len(a),-1);v=b.reshape(len(b),-1)
 dn=float(np.sqrt(np.sum(w[:,None]*d*d)));vn=float(np.sqrt(np.sum(w[:,None]*v*v)))
 return dict(relative_l2=dn/vn if vn>1e-20 else None,weighted_rms_difference=dn/float(np.sqrt(w.sum())),maximum_absolute_difference=float(abs(d).max()),reference_l2=vn)
def tensor(values,axes):
 nxx,nyy,nxy=values;axes=np.asarray(axes);assert np.max(abs(axes@axes.T-np.eye(3)))<1e-12
 return axes[:2].T@np.asarray([[nxx,nxy],[nxy,nyy]])@axes[:2]
# Independent closed-form Frobenius and rotated-axis checks for the comparison operator.
angle=.37;c,s=np.cos(angle),np.sin(angle);axes=np.asarray([[c,s,0],[-s,c,0],[0,0,1]])
value=tensor([2,3,4],axes);assert abs(float(np.sum(value*value))-(2**2+3**2+2*4**2))<1e-12
assert abs(metric([[1.,2.]],[[2.,4.]],[3.])['relative_l2']-.5)<1e-14
overlay=Overlay(*meshes,progress=lambda e,n,pairs:print(f'Overlay {e}/{n}, intersections {pairs}',flush=True))
origins=[np.asarray(r['reference_displacement_m']) for r in results]
previous=[np.zeros_like(v) for v in origins];rows=[]
for fa,fb in zip(results[0]['frames'],results[1]['frames'],strict=True):
 assert fa['name']==fb['name'];a,b=fa['state'],fb['state'];u,v=np.asarray(a['displacement_m']),np.asarray(b['displacement_m'])
 active_a,active_b=a['active_elements'],b['active_elements']
 row=dict(stage=fa['name'],from_seed=overlay.compare(u,v,active_a,active_b),from_construction_origin=overlay.compare(u-origins[0],v-origins[1],active_a,active_b),stage_increment=overlay.compare(u-previous[0],v-previous[1],active_a,active_b))
 ba,bb=a['beams'],b['beams'];assert not any(ba['active_members']) and not any(bb['active_members'])
 row['wall_node_coefficients']={k:metric(ba[k],np.asarray(bb[k])[node_map],np.ones(440)) for k in ['node_displacement_m','node_rotation_rad']}
 plate_results=[{r['id']:r for r in frame['beam_results'].get('plates',[])} for frame in [fa,fb]]
 assert all(len(pr)==sum(state['beams']['active_plates']) for pr,state in zip(plate_results,[a,b]))
 samples={k:([],[],[]) for k in ['membrane_force_kn_m','bending_moment_kn']}
 for key in sorted(keyed[0]):
  pa,pb=[keyed[i][key] for i in [0,1]];ra=plate_results[0].get(pa['id']);rb=plate_results[1].get(pb['id'])
  assert (ra is None)==(rb is None)
  if ra is None:continue
  assert abs(ra['area_m2']-rb['area_m2'])<1e-13
  coords=[]
  for cfg,plate,res in [(configs[0],pa,ra),(configs[1],pb,rb)]:
   nodes={n['id']:np.asarray(n['position_m']) for n in cfg['nodes']};tri=np.asarray([nodes[n] for n in plate['nodes']])
   coords.append(np.asarray([np.asarray(s['barycentric'])@tri for s in res['samples']]))
  dist,mapping=cKDTree(coords[1]).query(coords[0]);assert dist.max()<1e-12 and len(set(mapping))==len(coords[0])==3
  assert np.dot(ra['local_axes'][2],rb['local_axes'][2])>1-1e-12
  for i,j in enumerate(mapping):
   for name,(x,y,w) in samples.items():
    x.append(tensor(ra['samples'][i][name],ra['local_axes']));y.append(tensor(rb['samples'][j][name],rb['local_axes']));w.append(ra['area_m2']/3)
 row['wall_area_weighted_global_tensors']={k:metric(x,y,w) for k,(x,y,w) in samples.items() if x}
 rows.append(row);previous=[u,v]
assert all(sha(Path(q))==h for q,h in binding.items())
report=dict(comparison_recorded=True,left_elements=2865,right_elements=5144,wall_nodes=440,wall_plates=832,wall_spatial_topology_and_sections_and_loads_preserved=True,each_physical_case_bound_to_original_preparation_and_independent_frozen_audit=True,frames=rows,geometry=overlay.metadata,wall_translations_m_and_rotations_rad_compared_separately=True,wall_sample_coordinates_matched=True,wall_resultants_are_global_area_weighted_frobenius_norms=True,wall_node_norms_are_coefficients_not_surface_integrals=True,independent_solver_comparison=False,spatial_convergence_accepted=False,p6_complete=False,sha256=binding)
with (p/'comparison.json').open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2,allow_nan=False)
print(json.dumps([dict(stage=r['stage'],increment_relative_l2=r['stage_increment']['relative_l2']) for r in rows]),flush=True)
