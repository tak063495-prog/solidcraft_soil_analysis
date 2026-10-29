"""Extend failed rain validation without altering its inputs or thresholds."""
from pathlib import Path
from collections import Counter
import hashlib,json,shutil,sys,zipfile
p=Path(__file__).resolve().parent;root=p.parents[1];old=root/'runs/richards-local-rain-3d-20261029-v1'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
assert read(old/'raw-terminal.json')['passed'] and not read(old/'study-terminal.json')['passed']
original=read(old/'comparison.json');assert sum(not c['passed'] for c in original['checks'])==5
for name in ['source-snapshot.json','source-snapshot.zip','reference.py','reference-case.json']:
    shutil.copy2(old/name,p/name);assert sha(old/name)==sha(p/name)
runtime=p/'runtime';runtime.mkdir()
with zipfile.ZipFile(p/'source-snapshot.zip') as z:z.extractall(runtime)
assert len(read(p/'source-snapshot.json'))==611
(p/'solver').mkdir();shutil.copy2(old/'solver/soil-fem-3d.exe',p/'solver/soil-fem-3d.exe')
assert sha(p/'solver/soil-fem-3d.exe')=='d55ffe42c3b61e0891fec513489c9f5881830733ac9e003a01a0cf69f495905d'
sys.path.insert(0,str(runtime/'tools'));from richards_dam_mesh import make
n=24;mesh=make(n);faces=Counter(tuple(sorted(t[i] for i in f)) for t in mesh['tetrahedra'] for f in [(1,2,3),(0,3,2),(0,1,3),(0,2,1)])
mesh['face_sets']=dict(bottom=[],rain=[],dry=[])
for face,count in faces.items():
    if count!=1:continue
    xyz=[mesh['nodes'][i] for i in face]
    if all(q[2]==0 for q in xyz):mesh['face_sets']['bottom'].append(list(face))
    elif all(q[2]==1 for q in xyz):
        patch=all(.25<=q[0]<=.75 and .25<=q[1]<=.75 for q in xyz)
        mesh['face_sets']['rain' if patch else 'dry'].append(list(face))
assert len(mesh['tetrahedra'])==82944 and len(mesh['nodes'])==15625 and len(faces)==169344
case=read(old/'rt0-16-dt100/case.json');assert case['settings']==dict(initial_time_step_s=100.,maximum_time_step_s=100.,absolute_flow_tolerance_m3_s=1e-17,store_integration_points=False)
folder=p/'rt0-24-dt100';folder.mkdir()
for name,value in [('mesh.json',mesh),('case.json',case)]:
    with (folder/name).open('x') as f:json.dump(value,f,allow_nan=False)
assert case==read(old/'rt0-16-dt100/case.json')
reused=['comparison.json','study-terminal.json','scope.json','reference-case.json','fv-16-dt100.json','fv-32-dt100.json','rt0-8-dt100/result.json','rt0-16-dt100/result.json','rt0-16-dt100/audit.json']
# Hash large inputs in streaming mode to avoid a second in-memory field copy.
def stream_sha(q):
    h=hashlib.sha256()
    with q.open('rb') as f:
        for b in iter(lambda:f.read(4<<20),b''):h.update(b)
    return h.hexdigest()
with (p/'reused-parent-files.json').open('x') as f:json.dump({name:stream_sha(old/name) for name in reused},f,indent=2)
protocol=dict(schema='soil-fem-local-rain-refinement/1',parent='richards-local-rain-3d-20261029-v1',parent_failed_checks=[c for c in original['checks'] if not c['passed']],
    physics_periods_and_materials_and_boundaries_unchanged=True,solver_sha256=sha(p/'solver/soil-fem-3d.exe'),original_failures_and_thresholds_preserved=True,
    native_space_window=[8,16,24],reference_space_window=[16,32,64],dt_s=100,
    planned_native_time_subdivisions=16,planned_reference_time_subdivisions=32,planned_time_step_s=[100,50,25],
    gates=read(old/'scope.json')['predeclared_gates'],additional_finest_time_difference_rms_gate_m=.0005,
    phase1=['Native24 first accepted checkpoint then fresh-process remaining transient solve','Independent FV64/dt100 and FV32/dt50,25'],
    phase2=['Independent complete discrete audit; any faster audit implementation must first match the original on retained controls','Native16 time refinements','Full space/time and independent PDE comparisons under unchanged original thresholds'],
    native_elements=82944,native_nodes=15625,total_mesh_faces=169344,existing_default_maximum_free_faces=300000,
    raw_output_estimate_bytes=int((old/'rt0-16-dt100/result.json').stat().st_size*24**3/16**3),raw_audit_memory_and_runtime_require_separate_review=True,
    phase1_raw_completion_is_not_physical_acceptance=True,mechanics_or_stability_certification=False,p6_complete=False)
with (p/'protocol.json').open('x') as f:json.dump(protocol,f,indent=2)
print('Prepared82944 Tet4/15625 nodes, same d55 solver/physics/tolerances; original5 failed gates retained',flush=True)
