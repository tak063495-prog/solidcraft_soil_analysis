"""Independently audit the actual four-body seed, loads and accepted histories."""
from pathlib import Path
import hashlib,json,os,sys
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import numpy as np
p=Path(__file__).resolve().parent;root=p.parents[2];sys.path.insert(0,str(p/'runtime/tools'))
import fem_reference as ref
from mcc_rate_fem import physical
read=lambda path:json.loads(path.read_text(encoding='utf-8-sig'))
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
analysis=p/'workflow/analysis';files=[analysis/'ground/mesh.json',analysis/'ground/ground.json',analysis/'stages.json',analysis/'stage-results.json',p/'preparation.json',p/'mechanical-acceptance.json']
binding={str(path):sha(path) for path in files};mesh,ground,case,result,prepared,view=[read(path) for path in files]
assert result['complete'] and len(result['frames'])==5 and view['passed']
regions={r['id']:r for r in mesh['regions']};fields={f['body']:f for f in ground['initialization']['geostatic']['seed']['fields']};materials={m['id']:m['mcc'] for m in ground['materials']}
layout=result['integration_points'];owners=np.array([q['element_index'] for q in layout]);weights=np.array([q['weight_m3'] for q in layout]);material=np.asarray(mesh['material_ids'])[owners]
assert len(layout)==prepared['points'] and set(material)=={1,2}
geometry=ref.Geometry(mesh);position_error=0.;seed_error=0.;cursor=0
for e in range(len(mesh['tetrahedra'])):
 q=geometry.volume(e);count=len(q['weight']);actual=np.asarray([x['position_m'] for x in layout[cursor:cursor+count]])
 position_error=max(position_error,float(abs(q['position']-actual).max()));cursor+=count
for q,state in zip(layout,result['geostatic']['seed']['material_states'],strict=True):
 region=regions[mesh['region_ids'][q['element_index']]];f=fields[region['body']];m=materials[region['material_id']]
 stress=np.asarray(f['effective_stress_at_reference_compression_positive_kpa'])+np.asarray(f['effective_stress_gradient_kpa_per_m'])@(np.asarray(q['position_m'])-f['reference_position_m'])
 mean=stress[:3].mean();pc=f['pc_over_mean_effective']*mean
 void=m['specific_volume_ncl_at_reference']-m['lambda']*np.log(pc/m['reference_pressure_kpa'])+m['kappa']*np.log(pc/mean)-1
 expected=[*stress,pc,void,*([0.]*6)];actual=[*state['effective_stress_compression_positive_kpa'],state['preconsolidation_pressure_kpa'],state['void_ratio'],*state['plastic_strain_compression_positive']]
 seed_error=max(seed_error,float(abs(np.asarray(actual)-expected).max()))
 assert abs(stress[0]/stress[2]-.6)<1e-12 and abs(stress[1]/stress[2]-.6)<1e-12
assert position_error<1e-10 and seed_error<1e-10
seed_pc=np.array([s['preconsolidation_pressure_kpa'] for s in result['geostatic']['seed']['material_states']]);rows=[]
for frame in result['frames']:
 active=np.asarray(frame['state']['active_elements'])[owners];states=frame['state']['material_states'];pc=np.array([s['preconsolidation_pressure_kpa'] for s in states]);void=np.array([s['void_ratio'] for s in states]);plastic=np.array([s['plastic_strain_compression_positive'] for s in states]);stats={}
 for mid in (1,2):
  mask=active&(material==mid);w=weights[mask]
  stats[str(mid)]=dict(active_volume_m3=float(w.sum()),active_points=int(mask.sum()),plastic_history_points=int(np.count_nonzero(np.max(abs(plastic[mask]),axis=1)>1e-12)),maximum_pc_over_seed=float((pc[mask]/seed_pc[mask]).max()),minimum_pc_over_seed=float((pc[mask]/seed_pc[mask]).min()),volume_mean_void_ratio=float(w@void[mask]/w.sum()))
 rows.append(dict(stage=frame['name'],materials=stats))
rejected=[a for a in result['attempts'] if not a['accepted']];assert all(a['committed_state_unchanged_on_rejection'] for a in rejected)
last=result['frames'][-1]['state'];prior=result['frames'][-2]['state'];frozen=0
for q,current,old in zip(layout,last['material_states'],prior['material_states'],strict=True):
 if not last['active_elements'][q['element_index']]:assert current==old;frozen+=1
assert frozen==27*len(mesh['excavation_sets']['Pit'])
active=np.array(last['active_elements']);topology=ref.topology(mesh);floor_nodes=set();surface_rows={}
for label,body,centroid in [('floor','Lower11',[5.,2.,3.]),('surrounding_surface','UpperSoil',[5.,2.,6.])]:
 weights_surface=np.zeros(len(mesh['nodes']));area=0.;moment=np.zeros(3);rid=next(r['id'] for r in mesh['regions'] if r['body']==body)
 for face in mesh['face_sets'][f'body:{body}:zmax']:
  candidates=[e for e in topology[tuple(sorted(face))] if mesh['region_ids'][e]==rid and active[e]]
  if not candidates:continue
  e=candidates[0];q=ref.face(mesh,e,face);da=np.linalg.norm(q['area_vector'],axis=1)
  np.add.at(weights_surface,q['nodes'],da@q['shape']);area+=da.sum();moment+=da@q['position']
  if label=='floor':floor_nodes.update(q['nodes'])
 assert abs(area-(8. if label=='floor' else 32.))<1e-10 and max(abs(moment/area-centroid))<1e-10
 u=np.asarray(last['displacement_m']);before=np.asarray(prior['displacement_m']);origin=np.asarray(result['reference_displacement_m'])
 surface_rows[label]=dict(area_m2=float(area),excavation_increment_mm=(weights_surface@(u-before)*1000/area).tolist(),cumulative_display_mm=(weights_surface@(u-origin)*1000/area).tolist())
pressure=physical(mesh,ground,active,100.,True)[0];zero=physical(mesh,ground,active,0.,True)[0];loaded_floor=physical(mesh,ground,active,100.,False)[0]
assert np.max(abs((pressure-zero)[sorted(floor_nodes)]))<1e-10
np.testing.assert_allclose((pressure-zero).sum(axis=0),[0.,0.,-3200.],rtol=0,atol=1e-8)
np.testing.assert_allclose((loaded_floor-pressure).sum(axis=0),[0.,0.,-800.],rtol=0,atol=1e-8)
assert binding=={n:sha(Path(n)) for n in binding};binding.update({str(path):sha(path) for path in (Path(__file__),Path(ref.__file__))})
with (p/'field-verification.json').open('x') as f:json.dump(dict(seed_load_and_recorded_cutback_gates_passed=True,all_seed_states_independently_evaluated=True,seed_state_count=len(layout),maximum_gp_position_difference_m=position_error,maximum_seed_state_difference=seed_error,materials=2,material_bodies=10,excavation_selector_bodies=1,frames=rows,removed_gp_histories_exact=frozen,rejected_trials_restored=len(rejected),surface_measures=surface_rows,maximum_retained_nodal_y_displacement_m=float(np.max(abs(np.asarray(last['displacement_m'])[sorted(set(np.asarray(mesh['tetrahedra10'])[active].ravel())),1]))),floor_prescribed_pressure_kpa=0.,retained_top_pressure_force_kn=[0.,0.,-3200.],loaded_floor_control_additional_force_kn=[0.,0.,-800.],sha256=binding,independent_nonlinear_solution_comparison_complete=False,spatial_convergence_accepted=False,p6_complete=False),f,indent=2,allow_nan=False)
print(json.dumps(dict(passed=True,maximum_seed_error=seed_error,surface_measures=surface_rows,material_frames=rows)),flush=True)
