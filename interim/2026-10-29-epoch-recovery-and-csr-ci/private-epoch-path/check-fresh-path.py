"""New process tests for an accepted and a rejected numerical path prefix."""
from pathlib import Path
import datetime,hashlib,json,math,os,subprocess,sys,time
p=Path(__file__).resolve().parent;runtime=p/'runtime';out=p/'Fresh adaptive path with spaces';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
write=lambda q,v:q.write_text(json.dumps(v,indent=2,allow_nan=False),'utf-8')
assert read(p/'quality-v2/terminal.json')['passed']
write(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,
    script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
sys.path.insert(0,str(runtime/'tools'));from freeze_sources import manifest
frozen=read(p/'quality-v2/source-snapshot.json');assert manifest()==frozen
deps=runtime/'target/release/deps';libraries={};records=[]
for name in ['soil_fem_3d','serde_json','serde']:
    matches=list(deps.glob('lib'+name+'-*.rlib'));assert len(matches)==1
    libraries[name]=dict(path=str(matches[0]),sha256=sha(matches[0]))
def call(folder,label,args,expected=0):
    if (p/'pause-requested').exists():raise SystemExit(2)
    assert manifest()==frozen
    command=list(map(str,args));start=time.monotonic()
    with (folder/(label+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(command,cwd=runtime,stdout=log,stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        write(folder/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],
            arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code=child.wait()
    row=dict(step=label,exit_code=code,expected_exit_code=expected,elapsed_seconds=time.monotonic()-start)
    records.append(row);write(folder/(label+'-terminal.json'),row);print(json.dumps(row),flush=True);assert code==expected
exe=out/'mc-path-verification.exe'
command=['C:/Users/link_/.cargo/bin/rustc.exe','--edition','2024','-L','dependency='+str(deps),p/'fresh-path-driver.rs','-o',exe]
for name,v in libraries.items():command+=['--extern',name+'='+v['path']]
call(out,'compile-driver',command);image=sha(exe)
write(out/'library-and-image-identity.json',dict(libraries=libraries,driver_sha256=image,
    driver_source_sha256=sha(p/'fresh-path-driver.rs'),source_files=616,verification_only_not_production_cli=True))
cases=[]
for tension in [False,True]:
    folder=out/('rejected-material-prefix' if tension else 'accepted-plastic-prefix');folder.mkdir()
    mesh=read(runtime/'examples/mc-column-mesh.json');ground=read(runtime/'examples/ground-geostatic-mc-column.json')
    case=read(runtime/'examples/stages-geostatic-column.json')
    for m in ground['materials']:
        m['unit_weight_kn_m3']=1e-10;m['mohr_coulomb']=dict(young_kpa=30000.,poisson=.25,
            cohesion_kpa=20. if tension else 40.,friction_deg=30. if tension else 0.,dilation_deg=0.,
            tension_cutoff_kpa=None,flow='non_associated')
    if tension:
        constraints=[]
        for axis in range(3):
            direction=[0.,0.,0.];direction[axis]=1.;gradient=[0.,0.,0.];gradient[axis]=1./6000.
            constraints.append(dict(node_set='all',direction=direction,value_m=0.,gradient=gradient,reference_point_m=[0.,0.,0.]))
        case['boundary']=dict(constraints=constraints)
        case['stages']=[dict(kind='boundary',name='extend',target_factor=1.)]
    else:case['stages']=[dict(kind='pressure',name='load',target_pressure_kpa=100.)]
    target=[dict(m['mohr_coulomb'],cohesion_kpa=0. if tension else 20.) for m in ground['materials']]
    spec=dict(schema='soil-fem-mc-epoch-path/1',target_parameters=target,
        initial_increment=1. if tension else .25,minimum_increment=1./1024.,maximum_attempts=100)
    for n,v in [('mesh.json',mesh),('ground.json',ground),('case.json',case),('spec.json',spec)]:write(folder/n,v)
    call(folder,'origin',[exe,'origin',folder/'mesh.json',folder/'ground.json',folder/'case.json',folder/'origin.json'])
    call(folder,'continuous',[exe,'start-path',folder/'mesh.json',folder/'origin.json',folder/'spec.json','100',folder/'continuous.json'])
    call(folder,'prefix',[exe,'start-path',folder/'mesh.json',folder/'origin.json',folder/'spec.json','1' if tension else '2',folder/'prefix.json'])
    protected={n:sha(folder/n) for n in ['mesh.json','ground.json','case.json','spec.json','origin.json','prefix.json']}
    prefix=read(folder/'prefix.json');state=prefix['payload']['continuation']
    assert state['progress']==(0. if tension else .5)
    if tension:
        assert state['attempts'][0]['failure_kind']=='material' and not state['attempts'][0]['accepted']
        assert state['attempts'][0]['before_state_sha256']==state['attempts'][0]['after_state_sha256']
    call(folder,'resume',[exe,'resume-path',folder/'mesh.json',folder/'origin.json',folder/'prefix.json','100',folder/'resumed.json'])
    assert (folder/'continuous.json').read_bytes()==(folder/'resumed.json').read_bytes()
    final=read(folder/'resumed.json')['payload']['continuation']
    if tension:
        assert final['stop_reason']=='cutback_exhausted'
        limit=1.-10.*math.tan(math.radians(30.))/20.
        assert final['progress']<=limit+1e-10 and limit-final['progress']<2./1024.
    else:assert final['progress']==1. and final['stop_reason'] is None and len(final['accepted'])==4
    for attempt in final['attempts']:
        if not attempt['accepted']:
            assert attempt['before_state_sha256']==attempt['after_state_sha256']
            assert attempt['before_parameters_sha256']==attempt['after_parameters_sha256']
    tamper=read(folder/'prefix.json');tamper['payload']['continuation']['next_increment']*=.75
    write(folder/'tampered.json',tamper)
    call(folder,'reject-tamper',[exe,'resume-path',folder/'mesh.json',folder/'origin.json',folder/'tampered.json','100',folder/'must-not-exist.json'],1)
    assert not (folder/'must-not-exist.json').exists()
    assert all(sha(folder/n)==h for n,h in protected.items())
    assert read(folder/'origin.json')['provenance']['solver_binary_sha256']==image
    cases.append(dict(case=folder.name,passed=True,final_progress=final['progress'],stop_reason=final['stop_reason'],
        accepted_epochs=len(final['accepted']),attempts=len(final['attempts']),fresh_process_checkpoint_byte_exact=True))
old=p/'Fresh process recovery with spaces'
call(out,'reject-old-image',[exe,'advance',old/'mesh.json',old/'origin.json','30',out/'wrong-image-must-not-exist.json'],1)
assert not (out/'wrong-image-must-not-exist.json').exists()
assert all(sha(Path(v['path']))==v['sha256'] for v in libraries.values()) and sha(exe)==image and manifest()==frozen
write(out/'terminal.json',dict(passed=True,cases=cases,verification_executable_sha256=image,
    same_exact_image_for_origin_and_each_new_process=True,no_identity_hash_rewriting=True,
    real_rejected_material_prefix_and_accepted_plastic_prefix_resumed=True,
    checkpoint_byte_exact_to_continuous_for_both_cases=True,tamper_and_old_image_rejected_without_output=True,
    source_files=616,verification_only_not_production_cli=True,full296_171_and_production_cli_pending=True,
    stability_classification='not_evaluated',records=records,p6_complete=False))
print('PASS fresh adaptive path processes: exact accepted and rejected prefixes, no Fs/collapse classification',flush=True)
