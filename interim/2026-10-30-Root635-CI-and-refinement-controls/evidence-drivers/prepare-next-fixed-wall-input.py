"""Prepare a finer soil mesh while keeping the same wall and all physics."""
from pathlib import Path
import ast,datetime,hashlib,json,os,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1]
out=p/'fixed-wall-finer-input-v1'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
if (p/'pause-requested').exists():raise SystemExit(2)
old=p/'global-volume-sensitivity-input-v1'
assert read(old/'wall-spatial-topology-proof.json')['wall_mesh_equivalent_within_tolerance']
assert read(p/'mesh14909/study-terminal.json')['passed']
source=(old/'run.py').read_text('utf-8')
replacements={'inside_size_m=.2':'inside_size_m=.15',
    '"Mesh.MeshSizeMin",200.':'"Mesh.MeshSizeMin",150.',
    "\\'VIn\\':200.":"\\'VIn\\':150.",
    'local_box_size_m=.2':'local_box_size_m=.15'}
for a,b in replacements.items():
    assert source.count(a)==1,(a,source.count(a))
    source=source.replace(a,b,1)
ast.parse(source)
out.mkdir()
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,
    script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
dump(out/'protocol.json',dict(prepared=True,original_preparer=str(old/'run.py'),
    original_preparer_sha256=sha(old/'run.py'),inside_target_m=.15,outside_target_m=.5,
    fixed_surface_target_m=1.,transition_m=.5,
    original_private600_geometry_preparation_used_for_like_for_like_inputs=True,
    all_material_boundary_gravity_load_and_stage_settings_retained=True,
    fixed_wall_all440_nodes_and832_oriented_plates_must_match_exactly=True,
    no_new_solver_job_or_nonlinear_solve_by_this_preparer=True,
    this_input_is_not_a_convergence_result=True,p6_complete=False))
(out/'run.py').write_text(source,'utf-8',newline='\n')
shutil.copyfile(old/'run.py',out/'original-preparer.py')
proof=(old/'verify-wall.py').read_text('utf-8')
old_label="soil8061->next refinement"
assert proof.count(old_label)==1
proof=proof.replace(old_label,'same14909 reference -> finer prepared soil mesh')
ast.parse(proof)
(out/'verify-wall.py').write_text(proof,'utf-8',newline='\n')
records=[]
for label,script in [('prepare','run.py'),('fixed-wall-proof','verify-wall.py')]:
    if (p/'pause-requested').exists():raise SystemExit(2)
    while True:
        if (p/'pause-requested').exists():raise SystemExit(2)
        memory=json.loads(subprocess.check_output(['pwsh','-NoProfile','-Command',
            'Get-CimInstance Win32_OperatingSystem | Select-Object FreePhysicalMemory,FreeVirtualMemory | ConvertTo-Json -Compress'],
            encoding='utf-8-sig',creationflags=subprocess.CREATE_NO_WINDOW))
        if min(memory['FreePhysicalMemory'],memory['FreeVirtualMemory'])*1024>=12*(1<<30) and shutil.disk_usage(out).free>=32*(1<<30):break
        time.sleep(5)
    command=[sys.executable,'-I','-X','utf8',str(out/script)];start=time.monotonic()
    with (out/(label+'-driver.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(command,cwd=root,stdout=log,stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        dump(out/(label+'-driver-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],
            arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code=child.wait()
    row=dict(step=label,exit_code=code,elapsed_seconds=time.monotonic()-start)
    records.append(row);dump(out/(label+'-driver-terminal.json'),row)
    print(json.dumps(row),flush=True)
    assert code==0,label
prep=read(out/'preparation.json');proof=read(out/'wall-spatial-topology-proof.json')
assert prep['elements']>16088 and proof['wall_mesh_equivalent_within_tolerance']
assert read(out/'preparation-only.json')['solver_not_started']
assert proof['structural_node_position_maximum_matching_distance_m']==0.
assert proof['left_plates']==proof['right_plates']==832
assert proof['left_structural_nodes']==proof['right_structural_nodes']==440
dump(out/'terminal.json',dict(passed=True,records=records,elements=prep['elements'],nodes=prep['nodes'],
    integration_points=prep['points'],all440_wall_nodes_and832_oriented_plates_exact=True,
    all_original_physics_and_case_settings_retained=True,
    solver_not_started=True,next_exact_job_registration_and_resource_gate_still_required=True,
    general_space_convergence_not_established_by_input_preparation=True,p6_complete=False))
print(json.dumps(dict(prepared=True,elements=prep['elements'],wall_exact=True,solver_not_started=True)),flush=True)
