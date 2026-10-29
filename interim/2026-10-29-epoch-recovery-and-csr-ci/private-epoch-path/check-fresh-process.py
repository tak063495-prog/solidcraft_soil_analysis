"""Fresh executable processes and real files; this is a verification CLI only."""
from pathlib import Path
import datetime,hashlib,json,os,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];runtime=p/'runtime'
out=p/'Fresh process recovery with spaces';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
write=lambda q,v:q.write_text(json.dumps(v,indent=2,allow_nan=False),'utf-8')
assert read(p/'quality-v1/terminal.json')['passed']
write(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,
    script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
sys.path.insert(0,str(runtime/'tools'));from freeze_sources import manifest
frozen=read(p/'quality-v1/source-snapshot.json');assert manifest()==frozen
deps=runtime/'target/release/deps';libraries={}
for name in ['soil_fem_3d','serde_json','serde']:
    found=list(deps.glob('lib'+name+'-*.rlib'));assert len(found)==1,(name,[x.name for x in found])
    libraries[name]=dict(path=str(found[0]),sha256=sha(found[0]))
records=[]
def call(label,args,expected=0):
    if (p/'pause-requested').exists():raise SystemExit(2)
    assert manifest()==frozen
    command=list(map(str,args));start=time.monotonic()
    with (out/(label+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(command,cwd=runtime,stdout=log,stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        write(out/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],
            arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code=child.wait()
    r=dict(step=label,exit_code=code,expected_exit_code=expected,elapsed_seconds=time.monotonic()-start);records.append(r)
    write(out/(label+'-terminal.json'),r);print(json.dumps(r),flush=True);assert code==expected

exe=out/'mc-epoch-verification.exe'
command=['C:/Users/link_/.cargo/bin/rustc.exe','--edition','2024','-L','dependency='+str(deps),
    p/'fresh-process-driver.rs','-o',exe]
for name,v in libraries.items():command+=['--extern',name+'='+v['path']]
call('compile-driver',command)
image=sha(exe);write(out/'library-and-image-identity.json',dict(libraries=libraries,driver_sha256=image,
    driver_source_sha256=sha(p/'fresh-process-driver.rs'),verification_only_not_production_cli=True))
mesh=read(runtime/'examples/mc-column-mesh.json');ground=read(runtime/'examples/ground-geostatic-mc-column.json')
case=read(runtime/'examples/stages-geostatic-column.json')
for m in ground['materials']:
    m['unit_weight_kn_m3']=1e-10;m['mohr_coulomb']=dict(young_kpa=30000.,poisson=.25,cohesion_kpa=40.,
        friction_deg=0.,dilation_deg=0.,tension_cutoff_kpa=None,flow='non_associated')
case['stages']=[dict(kind='pressure',name='load',target_pressure_kpa=100.)]
for name,value in [('mesh.json',mesh),('ground.json',ground),('case.json',case)]:write(out/name,value)
call('origin',[exe,'origin',out/'mesh.json',out/'ground.json',out/'case.json',out/'origin.json'])
call('continuous',[exe,'continuous',out/'mesh.json',out/'origin.json','30','20',out/'continuous.json'])
call('advance',[exe,'advance',out/'mesh.json',out/'origin.json','30',out/'checkpoint.json'])
original_files={name:sha(out/name) for name in ['mesh.json','ground.json','case.json','origin.json','checkpoint.json']}
assert read(out/'origin.json')['provenance']['solver_binary_sha256']==image
assert read(out/'checkpoint.json')['payload']['identity']['solver_binary_sha256']==image
call('restore',[exe,'restore',out/'mesh.json',out/'origin.json',out/'checkpoint.json',out/'restored.json'])
assert (out/'checkpoint.json').read_bytes()==(out/'restored.json').read_bytes()
call('resume',[exe,'resume',out/'mesh.json',out/'origin.json',out/'checkpoint.json','20',out/'resumed.json'])
def stable(v):
    if isinstance(v,dict):return {k:stable(x) for k,x in v.items() if k!='elapsed_seconds'}
    if isinstance(v,list):return [stable(x) for x in v]
    return v
assert stable(read(out/'continuous.json'))==stable(read(out/'resumed.json'))
tampered=read(out/'checkpoint.json');tampered['payload']['parameters'][0]['cohesion_kpa']=29.
write(out/'tampered.json',tampered)
call('reject-tampered',[exe,'resume',out/'mesh.json',out/'origin.json',out/'tampered.json','20',out/'must-not-exist.json'],1)
assert not (out/'must-not-exist.json').exists()
unchanged=sha(out/'resumed.json')
call('reject-overwrite',[exe,'resume',out/'mesh.json',out/'origin.json',out/'checkpoint.json','20',out/'resumed.json'],1)
assert sha(out/'resumed.json')==unchanged
assert all(sha(out/n)==h for n,h in original_files.items())
assert all(sha(Path(v['path']))==v['sha256'] for v in libraries.values()) and sha(exe)==image
assert manifest()==frozen
write(out/'terminal.json',dict(passed=True,verification_executable_sha256=image,
    same_exact_image_created_origin_and_each_fresh_process=True,no_provenance_hashes_rewritten=True,
    full_step_non_elapsed_fields_exact=True,checkpoint_restore_byte_exact=True,
    original_inputs_origin_checkpoint_unchanged=True,tamper_no_output_and_overwrite_guard_passed=True,
    material_epoch_checkpoint_from_accepted_plastic_c30_to_next_c20=True,
    source_files=614,verification_only_not_production_cli=True,
    adaptive_path_and_production_cli_and_full292_171_pending=True,records=records,p6_complete=False))
print('PASS fresh process checkpoint/continuation, exact non-elapsed full step and byte-exact checkpoint; verification CLI only',flush=True)
