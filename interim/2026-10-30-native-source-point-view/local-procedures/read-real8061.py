"""Complete the already-compiled public Native probe using independently read VTUs."""
from pathlib import Path
import datetime,hashlib,json,os,subprocess,sys,time
import meshio,numpy as np
p=Path(__file__).resolve().parent;root=p.parents[1];out=p/'public8061-reader-v1'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
 with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
 with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2)
assert read(p/'actual8061-v2/terminal.json')['passed']
assert read(p/'actual8061-v2/verification.json')['passed']
binding=read(out/'compiled-probe-binding.json');binary=out/'public-source-point-reader.exe'
assert sha(binary)==binding['binary_sha256'] and sha(out/'source-point-reader-probe.rs')==binding['probe_sha256']
sources=read(p/'native-quality-v2/source73.json')['files_sha256'];native=root/'runs/native-typed-source-integration-20261029-v1/native'
assert all(sha(native/name)==h for name,h in sources.items())
workflow=p/'actual8061-v2/workflow';manifest=workflow/'view/view-manifest.json';document=workflow/'conformal/source.solidcraft.json'
original={str(q):sha(q) for q in [manifest,document,workflow/'analysis/stage-results.json',workflow/'view-verification.json',binary]}
results=workflow/'analysis/results';plan=[];inventory=read(manifest)
assert inventory['schema']=='soil-fem-solidcraft-view/4' and len(inventory['frames'])==5
for i,entry in enumerate(inventory['frames']):
 raw=meshio.read(results/f'integration-points-{i:03}.vtu');vtk=meshio.read(results/f'stage-{i:03}.vtu')
 ids=np.asarray(raw.point_data['integration_point_id'],dtype=int);materials=np.asarray(raw.point_data['material_id'],dtype=int)
 chosen=set()
 for material in np.unique(materials):
  candidates=ids[materials==material];chosen.update(map(int,candidates[[0,len(candidates)//2,-1]]))
 active=set(map(int,ids));inactive=next((v for v in range(217647) if v not in active),None)
 nodes=np.asarray(vtk.point_data['original_node_id'],dtype=int);frame=read(workflow/'view'/entry['frame'])
 plates=[k for k,value in enumerate(frame['structures']['plates']) if value['active']]
 plan.append(dict(index=i,point_ids=sorted(chosen),node_ids=sorted(set(map(int,nodes[[0,len(nodes)//2,-1]]))),
  inactive_point=inactive,plate_ids=sorted(set([plates[0],plates[-1]])) if plates else []))
dump(out/'query-plan.json',dict(frames=plan,selection='first/middle/last active GP per material and active node; first/last active plate; inactive GP refusal'))
dump(out/'data-driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
command=[str(binary),str(document),str(manifest),str(out/'query-plan.json'),str(out/'reader-receipt.json')];start=time.monotonic()
with (out/'actual-reader.log').open('x',encoding='utf-8') as log:
 child=subprocess.Popen(command,cwd=root,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
 dump(out/'actual-reader-child.json',dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()));code=child.wait()
dump(out/'actual-reader-terminal.json',dict(exit_code=code,elapsed_seconds=time.monotonic()-start));assert code==0
receipt=read(out/'reader-receipt.json');assert len(receipt['frames'])==5
raw_keys=['effective_stress_compression_positive_kpa','plastic_strain_compression_positive','void_ratio','preconsolidation_pressure_kpa']
checks=0;derived_checks=0;plate_checks=0;inactive_checks=0
for i,(actual,request) in enumerate(zip(receipt['frames'],plan,strict=True)):
 status=actual['configured'];assert status['raw_point_storage']=='authoritative_source'
 assert status['source_storage']['representation']=='typed_projected_MCC_selected_frame'
 assert status['source_storage']['material_history_count']==status['source_storage']['integration_layout_count']==217647
 assert status['triangles']>0 and actual['section']['section_triangles']>0 and actual['section']['reference_section_area_m2']>0
 raw=meshio.read(results/f'integration-points-{i:03}.vtu');vtk=meshio.read(results/f'stage-{i:03}.vtu')
 lookup={int(v):row for row,v in enumerate(raw.point_data['integration_point_id'])}
 for entry in actual['points']:
  value=entry['value'];row=lookup[entry['id']]
  assert value['raw_point_storage']=='authoritative_source' and value['element']==int(raw.point_data['original_element_id'][row])
  assert np.array_equal(value['position_m'],raw.points[row])
  assert value['weight_m3']==raw.point_data['integration_weight_m3'][row]
  assert value['material_id']==int(raw.point_data['material_id'][row]) and value['stored_model_code']==0
  for key in raw_keys:assert np.array_equal(value[key],raw.point_data[key][row]),(i,entry['id'],key);checks+=1
  for key,item in value['mixed_fields'].items():
   assert key in raw.point_data and np.allclose(item,raw.point_data[key][row],rtol=1e-11,atol=1e-11),(i,entry['id'],key)
   derived_checks+=1
 first=actual['points'][0]['value'];after=actual['raw_after_settings']
 for key in [*raw_keys,'position_m','element','weight_m3','material_id','stored_model_code','mixed_fields']:assert after[key]==first[key],(i,key)
 if request['inactive_point'] is not None:assert actual['inactive_point_rejected'] is True;inactive_checks+=1
 node_lookup={int(v):row for row,v in enumerate(vtk.point_data['original_node_id'])}
 for entry in actual['nodes']:
  assert np.array_equal(entry['value']['stored_node_fields']['total_displacement_from_seed_m'],vtk.point_data['total_displacement_from_seed_m'][node_lookup[entry['id']]])
  checks+=1
 if actual['plates']:
  plate_raw=meshio.read(results/f'plate-samples-{i:03}.vtu')
  for entry in actual['plates']:
   member=entry['value']['member'];assert member['active']
   for sample,data in enumerate(member['result']['samples']):
    rows=np.flatnonzero((plate_raw.point_data['plate_index']==entry['id'])&(plate_raw.point_data['sample_index']==sample));assert len(rows)==1
    for key in ['membrane_strain','curvature_per_m','membrane_force_kn_m','bending_moment_kn','drilling_mismatch_rad']:
     assert np.array_equal(data[key],plate_raw.point_data[key][rows[0]]),(i,entry['id'],sample,key);plate_checks+=1
assert original=={q:sha(Path(q)) for q in original} and all(sha(native/name)==h for name,h in sources.items())
dump(out/'terminal.json',dict(passed=True,source_files=73,frames=5,actual_public_Native_API_and_original_CAD_identity_used=True,
 original1459184313byte_source_unchanged=True,representative_raw_and_node_checks=checks,representative_projected_GP_checks=derived_checks,
 representative_raw_plate_checks=plate_checks,inactive_point_refusal_checks=inactive_checks,
 original_cell_source_and_projected_checks_retained=True,raw_values_unchanged_by_section_deformation_and_quantity=True,
 GUI_full_CI_and_general_mesh_acceptance_separate=True,p6_complete=False))
print(json.dumps(dict(passed=True,raw_and_node_checks=checks,derived_checks=derived_checks,plate_checks=plate_checks)),flush=True)
