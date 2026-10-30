"""Full unchanged mechanics/VTU audit, then the predeclared2D/3D field gates."""
from pathlib import Path
import contextlib,datetime,hashlib,json,math,os,shutil,subprocess,sys,time
import numpy as np
p=Path(__file__).resolve().parent
study=p/'plane-strain-slope-comparison-v1'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
assert len(sys.argv)==2 and sys.argv[1] in ['xy1-ny1-global64','xy1-ny1-local-half']
folder=study/sys.argv[1]
out=folder/'field-comparison-v1';out.mkdir()
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,
    script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
terminal=folder/'solve-v1/terminal.json'
while not terminal.exists():
    if (p/'pause-requested').exists():raise SystemExit(2)
    time.sleep(5)
solved=read(terminal)
if not solved['passed']:
    dump(out/'terminal.json',dict(passed=False,reason='actual native solve not complete; no field comparison or physical Fs upper inferred',
        actual_solver_terminal=solved,p6_complete=False))
    raise SystemExit(1)
if folder.name=='xy2-ny2':
    # Avoid overlapping the largest additional DOM/VTU audit with the original
    # 16k complete production/export audit. This does not change field criteria.
    for gate in [p/'geometry-audit-real8061-v1/terminal.json',p/'global-volume-sensitivity-study-v1/terminal.json']:
        while not gate.exists():
            if (p/'pause-requested').exists():raise SystemExit(2)
            time.sleep(5)
        assert read(gate)['passed'],str(gate)
required=(16 if folder.name=='xy2-ny2' else 8)*(1<<30)
while True:
    if (p/'pause-requested').exists():raise SystemExit(2)
    memory=json.loads(subprocess.check_output(['pwsh','-NoProfile','-Command',
        'Get-CimInstance Win32_OperatingSystem | Select-Object FreePhysicalMemory,FreeVirtualMemory | ConvertTo-Json -Compress'],
        encoding='utf-8-sig',creationflags=subprocess.CREATE_NO_WINDOW))
    if min(memory['FreePhysicalMemory'],memory['FreeVirtualMemory'])*1024>=required and shutil.disk_usage(folder).free>=16*(1<<30):break
    time.sleep(5)
dump(out/'resource-preflight.json',dict(passed=True,free_physical_bytes=memory['FreePhysicalMemory']*1024,
    free_virtual_bytes=memory['FreeVirtualMemory']*1024,minimum_each_free_bytes=required))
runtime=p/'source635-offline-setup-v2/runtime'
assert read(p/'source635-offline-setup-v2/quality-v1/terminal.json')['passed']
sys.path.insert(0,str(runtime/'tools'))
from freeze_sources import manifest as runtime_manifest
runtime_sources=read(p/'source635-offline-setup-v2/quality-v1/source-snapshot.json')
assert runtime_manifest()==runtime_sources
reference=study/'reference-2d-v2'
source_binding=read(reference/'source-binding.json')['files_sha256']
assert all(sha(Path(q))==h for q,h in source_binding.items())
verified=read(reference/'terminal.json');assert verified['passed']
protocol=read(study/'protocol.json');limits=protocol['predeclared_acceptance']
assert protocol['reference_fixed_factor']==verified['fixed_factor_for_future3d_match']
prepared=read(folder/'prepared.json')
binding={str(folder/name):h for name,h in prepared['input_sha256'].items()}
binding.update({str(q):sha(q) for q in [Path(__file__),study/'protocol.json',reference/'terminal.json',
    reference/'source-binding.json',reference/'published-2d-result.json',reference/'reference-points.json',
    folder/'stage-results.json',folder/'solver-job/job.json']})
assert sha(folder/'stage-results.json')==solved['result_sha256']
assert all(sha(Path(q))==h for q,h in binding.items())
dump(out/'binding.json',binding)
dump(out/'audit-runtime635-source-snapshot.json',runtime_sources)
from export_stages import export
from verify_stage_package import verify
mp=folder/'mesh.json';gp=folder/'ground.json';rp=folder/'stage-results.json'
start=time.monotonic()
with (out/'independent-full-audit.log').open('x',encoding='utf-8') as log,contextlib.redirect_stdout(log):
    export(mp,rp,folder/'results')
    verify(mp,gp,rp,folder/'results',solver=folder/'solver-job/solver/soil-fem-3d.exe')
audit=read(folder/'stage-verification.json')
assert audit['schema']=='soil-fem-stage-verification/2' and audit['raw_integration_point_roundtrip']=='passed'
dump(out/'full-independent-audit-terminal.json',dict(passed=True,elapsed_seconds=time.monotonic()-start,
    original_audit_not_modified=True,frames=len(audit['frames']),all_raw_VTU_and_mechanical_geostatic_checks_passed=True))
mesh=read(mp);run=read(rp);mapping=read(folder/'reference-map.json')
target=read(reference/'published-2d-result.json')
assert run['complete'] and run['geostatic']['complete'] and run['geostatic']['progress']==1.
assert len(run['frames'])==len(audit['frames'])==2
state=run['geostatic']['last_committed']
u=np.asarray(state['displacement_m'],float)
reference_u=np.asarray(target['displacements_flat'],float).reshape(-1,3)[:,1:]
assert len(reference_u)==661 and np.isfinite(u).all()
norm_u=float(np.linalg.norm(reference_u));maximum_u=float(np.linalg.norm(reference_u,axis=1).max())
assert norm_u>0 and maximum_u>0
displacements={key:u[np.asarray(ids,dtype=int)][:,[0,2]] for key,ids in mapping['matching_3d_nodes_front_middle_back'].items()}
relative_u={key:float(np.linalg.norm(value-reference_u)/norm_u) for key,value in displacements.items()}
maximum_u_error=max(float(np.linalg.norm(value-reference_u,axis=1).max()/maximum_u) for value in displacements.values())
depth_variation=max(float(np.linalg.norm(a-b,axis=1).max()/maximum_u)
    for a in displacements.values() for b in displacements.values())
material_states=state['material_states']
layout=run['integration_points'];assert len(layout)==len(material_states)==prepared['integration_points']
sigma=np.asarray([s['effective_stress_compression_positive_kpa'] for s in material_states])
assert sigma.shape==(len(layout),6) and np.isfinite(sigma).all()
parents=np.asarray([mapping['original_quad_per_tet'][q['element_index']]-1 for q in layout],int)
weights=np.asarray([q['weight_m3'] for q in layout])
assert np.all(weights>0) and parents.min()==0 and parents.max()==199
volumes=np.bincount(parents,weights=weights,minlength=200)
area=np.asarray(mapping['source_quad_area_m2'])
assert np.allclose(volumes,area,rtol=0,atol=1e-10)
mean_native=np.column_stack([np.bincount(parents,weights=weights*sigma[:,j],minlength=200)/volumes for j in range(6)])
points=read(reference/'reference-points.json');assert len(points)==800
mean_reference=np.zeros((200,6))
for point in points:
    sx,sz,tau=point['stress_x_z_shear_kpa']
    mean_reference[point['element']-1]+=np.asarray([sx,point['out_of_plane_stress_kpa'],sz,0.,0.,tau])*point['weight_m2']
mean_reference/=area[:,None]
stress_relative=float(np.linalg.norm((mean_native-mean_reference)*np.sqrt(area[:,None]))/np.linalg.norm(mean_reference*np.sqrt(area[:,None])))
transverse_shear=float(np.linalg.norm(sigma[:,[3,4]]*np.sqrt(weights[:,None]))/np.linalg.norm(sigma*np.sqrt(weights[:,None])))
material=verified['material'];phi=math.radians(material['friction_deg']);cohesion=material['cohesion_kpa']
tensors=np.zeros((len(layout),3,3));tensors[:,0,0]=sigma[:,0];tensors[:,1,1]=sigma[:,1];tensors[:,2,2]=sigma[:,2]
for component,(a,b) in enumerate([(0,1),(1,2),(0,2)],3):tensors[:,a,b]=tensors[:,b,a]=sigma[:,component]
eigen=np.linalg.eigvalsh(tensors)
yield_value=eigen[:,-1]-eigen[:,0]-(eigen[:,-1]+eigen[:,0])*math.sin(phi)-2*cohesion*math.cos(phi)
yield_normalized=float(np.max(yield_value/np.maximum(1.,np.maximum(np.max(np.abs(eigen),axis=1),cohesion))))
frame=run['frames'][0]
force_error=float(np.linalg.norm(np.asarray(frame['support_reaction_kn'])+[0,0,-4400])/4400)
geostatic=run['geostatic']
rejected=[a for a in geostatic['attempts'] if not a['accepted']]
gates=dict(full_independent_audit=True,
    native_relative_physical_force=geostatic['final_relative_physical_force_residual']<=limits['relative_physical_force'],
    global_force=force_error<=limits['relative_global_force'],
    all_3d_MC_yield=yield_normalized<=limits['maximum_normalized_3d_yield'],
    all_661_displacements_each_depth=max(relative_u.values())<=limits['displacement_sample_relative_L2'],
    maximum_displacement_error=maximum_u_error<=limits['maximum_displacement_error_over_reference_max'],
    all_200_prism_mean_stresses=stress_relative<=limits['prism_mean_stress_relative_L2'],
    transverse_shear=transverse_shear<=limits['transverse_shear_over_stress_L2'],
    depth_variation=depth_variation<=limits['depth_displacement_variation_over_max'],
    all_transverse_displacements_fixed=np.max(np.abs(u[:,1]))==0.,
    strain_applicability=geostatic['maximum_absolute_engineering_strain_from_seed']<=limits['maximum_engineering_strain_from_seed'],
    actual_rejected_history_unchanged=all(a['committed_state_unchanged_on_rejection'] for a in rejected))
gates={key:bool(value) for key,value in gates.items()}
metrics=dict(displacement_sample_relative_L2_per_depth=relative_u,maximum_displacement_error_over_reference_max=maximum_u_error,
    prism_mean_stress_relative_L2=stress_relative,depth_displacement_variation_over_max=depth_variation,
    transverse_shear_over_stress_L2=transverse_shear,maximum_normalized_3d_yield=yield_normalized,
    relative_global_force_error=force_error,actual_geostatic_rejections=len(rejected),
    maximum_engineering_strain_from_seed=geostatic['maximum_absolute_engineering_strain_from_seed'])
dump(out/'field-samples.json',dict(source2d_u_xz_m=reference_u.tolist(),native3d_u_xz_m={k:v.tolist() for k,v in displacements.items()},
    source2d_prism_mean_stress_xyz_xy_yz_xz_kpa=mean_reference.tolist(),native3d_prism_mean_stress_xyz_xy_yz_xz_kpa=mean_native.tolist(),
    original200_quad_prism_volumes_m3=volumes.tolist()))
assert all(sha(Path(q))==h for q,h in binding.items()) and runtime_manifest()==runtime_sources
assert all(sha(Path(q))==h for q,h in source_binding.items())
report=dict(passed=all(gates.values()),case=folder.name,metrics=metrics,gates=gates,original_predeclared_limits=limits,
    protocol_sha256=sha(study/'protocol.json'),fixed_factor=protocol['reference_fixed_factor'],
    all_661_reference_nodes_at3depths_and_all200_prism_means_compared=True,
    actual_full_native3d_solve_not_synthetic_and_published_actual2d_reference_not_new_Excel_run=True,
    independent_source_and_raw_result_and_original_audit_unchanged=True,
    numerical_error_not_reinterpreted_as_physical_Fs_upper=True,
    general_mesh_convergence_finite_width_3D_IDH_SRM_and_allP6_still_separate=True,p6_complete=False)
dump(out/'assessment.json',report);dump(out/'terminal.json',report)
print(json.dumps(dict(passed=report['passed'],case=folder.name,metrics=metrics,gates=gates)),flush=True)
if not report['passed']:raise SystemExit(1)
