from pathlib import Path
import datetime,hashlib,json,os,shutil,sys,time
import meshio
import numpy as np
p=Path(__file__).resolve().parent;root=p.parents[1];runtime=p/'runtime'
sys.path.insert(0,str(runtime/'tools'));import export_mc_epoch as exporter
read=lambda q:json.loads(q.read_text('utf-8'));sha=exporter.sha
write=lambda q,v:q.write_text(json.dumps(v,indent=2,allow_nan=False),'utf-8')
write(p/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,
    script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
base=root/'runs/mc-epoch-cli-development-20261029-v1';out=p/'real-cases-v2';out.mkdir()
cases=[('factor1',base/'cad-cphi-paths/factor1.0/continuous'),
    ('factor1.5',base/'cad-cphi-paths/factor1.5/continuous'),
    ('plastic-pause',base/'cad-cphi-paths/factor1.5/paused'),
    ('plastic-recovered',base/'cad-cphi-paths/factor1.5/recovered'),
    ('state-reset-and-rebirth',base/'cad-controls/gate10/job')]
# The real nonconvergent parameter path must export only its last accepted point.
plastic=base/'cad-controls/plastic-preload'
for candidate in plastic.rglob('job.json'):
    job=candidate.parent
    if (job/'report.json').exists() and not read(job/'report.json')['complete']:
        cases.append(('incomplete-path',job));break
assert len(cases)==6
records=[]
for name,job in cases:
    if (p/'pause-requested').exists():raise SystemExit(2)
    fingerprints={q.relative_to(job).as_posix():sha(q) for q in job.rglob('*') if q.is_file()}
    start=time.monotonic();packet=exporter.export(job,out/name)
    original=read(job/'inputs/origin.json');pointer=read(job/'recovery.json')
    snapshot=read(job/pointer['checkpoint']);cp=snapshot['checkpoint']['payload'];material=cp['material_checkpoint']['payload']
    state=material['committed'];u=np.asarray(state['displacement_m']);base_u=np.asarray(original['last_committed']['displacement_m'])
    continuum=meshio.read(out/name/'current.vtu');gp=meshio.read(out/name/'current-points.vtu')
    nodes=continuum.point_data['original_node_id'];ids=continuum.cell_data['original_element_id'][0]
    points=gp.point_data['integration_point_id'];active=np.flatnonzero(state['active_elements'])
    assert np.array_equal(ids,active)
    assert np.array_equal(continuum.point_data['total_displacement_from_seed_m'],u[nodes])
    assert np.array_equal(continuum.point_data['epoch_displacement_increment_m'],(u-base_u)[nodes])
    for field in ('effective_stress_compression_positive_kpa','plastic_strain_compression_positive','void_ratio','equivalent_plastic_strain','plastic_work_kpa'):
        assert np.array_equal(gp.point_data[field],np.asarray([state['material_states'][i][field] for i in points]))
    assert (out/name/'accepted-snapshot.json').read_bytes()==(job/pointer['checkpoint']).read_bytes()
    assert fingerprints=={q.relative_to(job).as_posix():sha(q) for q in job.rglob('*') if q.is_file()}
    record=dict(name=name,passed=True,complete=packet['complete'],progress=packet['parameter_progress'],
        stop_reason=packet['stop_reason'],active_elements=len(ids),active_gp=len(points),
        maximum_equivalent_plastic_strain=float(gp.point_data['equivalent_plastic_strain'].max()),
        raw_field_and_increment_exact=True,source_all_files_unchanged=True,elapsed_seconds=time.monotonic()-start,
        job=str(job),source_inventory_sha256=hashlib.sha256(json.dumps(fingerprints,sort_keys=True).encode()).hexdigest())
    records.append(record);print(json.dumps(record),flush=True)
assert records[1]['maximum_equivalent_plastic_strain']>0 and records[2]['maximum_equivalent_plastic_strain']>0
assert records[2]['complete'] is False and records[-1]['stop_reason']=='cutback_exhausted'
for filename in ('current.vtu','current-points.vtu','accepted-snapshot.json'):
    assert (out/'factor1.5'/filename).read_bytes()==(out/'plastic-recovered'/filename).read_bytes() if filename!='accepted-snapshot.json' else True
negative=p/'negative-v2';negative.mkdir();rejections=[]
for mode in ('mesh-checksum','solver-checksum','snapshot-checksum','pointer-traversal','forged-state-and-file-checksum'):
    if (p/'pause-requested').exists():raise SystemExit(2)
    job=negative/mode;shutil.copytree(cases[2][1],job)
    pointer=read(job/'recovery.json');target=job/pointer['checkpoint']
    if mode=='mesh-checksum':
        with (job/'inputs/mesh.json').open('ab') as f:f.write(b' ')
    elif mode=='solver-checksum':
        with (job/'solver/soil-fem-3d.exe').open('ab') as f:f.write(b' ')
    elif mode=='snapshot-checksum':
        with target.open('ab') as f:f.write(b' ')
    elif mode=='pointer-traversal':
        pointer['checkpoint']='../outside.json';write(job/'recovery.json',pointer)
    else:
        snapshot=read(target)
        snapshot['checkpoint']['payload']['material_checkpoint']['payload']['committed']['displacement_m'][0][2]+=.01
        write(target,snapshot);pointer['checkpoint_sha256']=sha(target);write(job/'recovery.json',pointer)
    destination=negative/(mode+'-must-not-exist')
    try:exporter.export(job,destination)
    except ValueError as error:
        assert not destination.exists();rejections.append(dict(mode=mode,reason=str(error),no_export_created=True))
    else:raise AssertionError('tampered job accepted: '+mode)
write(p/'real-terminal-v2.json',dict(passed=True,cases=records,rejected=rejections,
    source_files_sha256={str(q.relative_to(runtime)):sha(q) for q in (runtime/'tools/export_mc_epoch.py',runtime/'tests/test_mc_epoch_export.py')},
    accepted_only=True,ordinary_run_not_forged=True,native_reader_pending=True,goal_complete=False,p6_complete=False))
print('PASS actual CAD plastic/paused/incomplete export and immutable/deep-check tamper refusals',flush=True)
