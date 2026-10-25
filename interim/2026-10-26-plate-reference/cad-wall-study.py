"""Run both CAD meshes and declared drilling variants with independent audits."""
import copy
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np

folder=Path(__file__).resolve().parent;root=folder.parents[1]
sys.path.insert(0,str(root/'tools'))
import freeze_sources
import fem_reference as ref

def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,v):
    with p.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)

freeze_sources.freeze(folder);sources=freeze_sources.manifest()
solver_source=root/'target/release/soil-fem-3d.exe';digest=sha(solver_source)
assert digest=='8a3eda114a0d92a46f5854b84e08c010928c513f7b2a0ccf0fc74dfd86da2e7d'
for label,alpha in [('coarse-drill-low',1e-6),('coarse-drill-high',.01)]:
    run=folder/label;run.mkdir();(run/'ground').mkdir()
    for name in ('mesh.json','ground.json'):shutil.copyfile(folder/'coarse/ground'/name,run/'ground'/name)
    case=read(folder/'coarse/stages.json');case['beam_system']['plate_sections'][0]['drilling_factor']=alpha
    write(run/'stages.json',case)
    write(run/'variant-binding.json',dict(base='coarse',base_stages_sha256=sha(folder/'coarse/stages.json'),drilling_factor=alpha,
        all_other_physical_inputs_unchanged=True))

def invoke(run,name,command):
    assert freeze_sources.manifest()==sources
    code=subprocess.call([sys.executable,'-X','utf8',str(root/'tools/run_logged.py'),str(run/(name+'.log')),*map(str,command)],cwd=root)
    if code:raise SystemExit(code)
    assert freeze_sources.manifest()==sources

for label in ('coarse','fine','coarse-drill-low','coarse-drill-high'):
    run=folder/label;(run/'solver').mkdir();solver=run/'solver/soil-fem-3d.exe';shutil.copyfile(solver_source,solver)
    write(run/'solver-binding.json',dict(solver_sha256=digest))
    invoke(run,'solve',[solver,'solve-stages',run/'ground/mesh.json',run/'ground/ground.json',run/'stages.json',run/'stage-results.json'])
    invoke(run,'export',[sys.executable,'-X','utf8',root/'tools/export_stages.py',run/'ground/mesh.json',run/'stage-results.json',run/'exports'])
    invoke(run,'force-vtu-audit',[sys.executable,'-X','utf8',root/'tools/verify_stage_package.py',run/'ground/mesh.json',run/'ground/ground.json',run/'stage-results.json',run/'exports','--solver',solver])
    invoke(run,'independent',[sys.executable,'-X','utf8',root/'tools/verify_beam_coupling.py',run,run/'independent.json'])
    result=read(run/'stage-results.json')
    assert result['complete'] and result['schema']=='soil-fem-stage-result/14'
    assert all(s['equivalent_plastic_strain']==0. and s['plastic_work_kpa']==0. and
        all(v==0. for v in s['plastic_strain_compression_positive']) for frame in result['frames'] for s in frame['state']['material_states'])
    print(json.dumps(dict(case=label,all_five_frames_independent_passed=True)),flush=True)

run=folder/'coarse';solver=run/'solver/soil-fem-3d.exe'
invoke(run,'history',[sys.executable,'-X','utf8',root/'tools/verify_beam_history.py','--solver',solver,'--run',run,'--output-dir',run/'history'])
invoke(run,'checkpoints',[sys.executable,'-X','utf8',root/'tools/verify_construction_checkpoint.py','--solver',solver,'--run',run,'--output-dir',run/'checkpoints'])

def quantities(label):
    run=folder/label;data=read(run/'stage-results.json');mesh=read(run/'ground/mesh.json');cfg=data['case']['beam_system']
    ids={n['id']:i for i,n in enumerate(cfg['nodes'])}
    before,after=data['frames'][1:3]
    delta=np.asarray(after['state']['beams']['node_displacement_m'])-before['state']['beams']['node_displacement_m']
    weight=np.zeros(len(ids));area=0.
    xyz=np.array([n['position_m'] for n in cfg['nodes']])
    for p in cfg['plates']:
        cell=[ids[n] for n in p['nodes']];a=np.linalg.norm(np.cross(xyz[cell[1]]-xyz[cell[0]],xyz[cell[2]]-xyz[cell[0]]))/2
        weight[cell]+=a/3;area+=a
    assert abs(area-80.)<1e-8
    plate_mean=float(weight@delta[:,0]/area)
    soil_delta=np.asarray(after['state']['displacement_m'])-before['state']['displacement_m'];value=0.;soil_area=0.
    topology=ref.topology(mesh)
    for face in mesh['face_sets']['body:Soil:xmax']:
        q=ref.face(mesh,topology[tuple(sorted(face))][0],face);w=np.linalg.norm(q['area_vector'],axis=1)
        value+=float(w@(q['shape']@soil_delta[q['nodes'],0]));soil_area+=float(w.sum())
    return dict(mesh=label,nodes=len(mesh['nodes']),tet10=len(mesh['tetrahedra']),plates=len(cfg['plates']),
        shared_wall_nodes=len(cfg['nodes']),area_m2=area,wall_plate_mean_loading_ux_m=plate_mean,
        soil_quadratic_trace_mean_loading_ux_m=value/soil_area,
        interpolation_mean_difference_m=plate_mean-value/soil_area)

metrics=[quantities(label) for label in ('coarse','fine')]
mesh_difference=abs(metrics[0]['wall_plate_mean_loading_ux_m']-metrics[1]['wall_plate_mean_loading_ux_m'])/abs(metrics[1]['wall_plate_mean_loading_ux_m'])
baseline=read(folder/'coarse/stage-results.json')
def fields(data):
    a,b=data['frames'][1:3]
    displacement=np.asarray(b['state']['displacement_m'])-a['state']['displacement_m']
    stress=np.array([s['effective_stress_compression_positive_kpa'] for s in b['state']['material_states']])-np.array([s['effective_stress_compression_positive_kpa'] for s in a['state']['material_states']])
    moments=np.array([s['bending_moment_kn'] for p in b['beam_results']['plates'] for s in p['samples']])-np.array([s['bending_moment_kn'] for p in a['beam_results']['plates'] for s in p['samples']])
    return displacement,stress,moments
base_fields=fields(baseline);sensitivity=[]
for label in ('coarse-drill-low','coarse-drill-high'):
    result=read(folder/label/'stage-results.json');values=fields(result)
    errors={name:float(np.linalg.norm(other-base)/max(np.linalg.norm(base),1e-12)) for name,base,other in zip(('displacement','stress','bending_moment'),base_fields,values,strict=True)}
    sensitivity.append(dict(case=label,relative_l2_errors=errors,one_percent_gate_satisfied=max(errors.values())<.01))
report=dict(schema='soil-fem-cad-wall-scoped-study/1',passed=True,solver_sha256=digest,source_files=len(sources),
    reference='Each variant independently solves original-input 3D soil+plate, retaining all plate rotations; no Rust solved field generates the reference',
    mesh_metrics=metrics,mesh_pair_relative_wall_loading_displacement_difference=mesh_difference,
    two_percent_mesh_pair_gate_satisfied=mesh_difference<.02,drilling_sensitivity=sensitivity,
    separate_process_checkpoints_on_coarse_only=5,exact_coarse_cutback_clean_path=True,
    exact_shared_nodes=True,pointwise_soil_plate_trace_identity=False,spatial_convergence_accepted=False,
    native_gui_wall_accepted=False,excavation_retaining_wall_accepted=False,goal_complete=False,p6_complete=False,
    evidence_sha256={str(p.relative_to(folder)):sha(p) for p in folder.rglob('*.json') if p.name in ('independent.json','verification.json','stage-verification.json','plate-surface-binding.json','source-snapshot.json')})
write(folder/'verification.json',report)
print(json.dumps(report),flush=True)
