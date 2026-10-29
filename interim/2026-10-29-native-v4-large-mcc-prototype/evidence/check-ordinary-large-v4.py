"""Wait for immutable Native67-v4 then check the actual 535 MB ordinary MCC view."""
from pathlib import Path
import datetime,hashlib,json,os,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];out=p/'ordinary-large-v4';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8'))
def sha(q):
    h=hashlib.sha256()
    with q.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def write(q,v):q.write_text(json.dumps(v,indent=2,allow_nan=False),'utf-8')
write(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,
    script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
source=p/'ordinary-reader-probe.rs';frozen_probe=sha(source)
while not (p/'desktop-and-reader-v4/terminal.json').exists():
    if (p/'pause-requested').exists():raise SystemExit(2)
    time.sleep(5)
assert read(p/'desktop-and-reader-v4/terminal.json')['passed']
identity=read(p/'desktop-and-reader-v4/binary-identities.json')
engine=Path(identity['engine']);serde=Path(identity['serde_json'])
assert sha(engine)==identity['engine_sha256'] and sha(serde)==identity['serde_json_sha256']
assert sha(source)==frozen_probe
case=root/'runs/mcc-large-structures-development-20261028-v1/large-case'
document=case/'source.solidcraft.json';manifest=case/'workflow/view/view-manifest.json'
raw=case/'workflow/analysis/stage-results.json'
before={str(q):sha(q) for q in [document,manifest,raw,source,Path(__file__)]}
assert before[str(raw)]=='750355444ed1c05325e0983dfbd006d90312f5ef049787f315f873043739e9c5'
def call(label,args):
    if (p/'pause-requested').exists():raise SystemExit(2)
    command=list(map(str,args));start=time.monotonic()
    with (out/(label+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(command,cwd=p/'native',stdout=log,stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        write(out/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],
            created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code=child.wait()
    write(out/(label+'-terminal.json'),dict(exit_code=code,elapsed_seconds=time.monotonic()-start))
    assert code==0
binary=out/'ordinary-large-reader.exe'
call('compile',['C:/Users/link_/.cargo/bin/rustc.exe','--edition','2024','-L','dependency='+str(engine.parent),
    '--extern','solidcraft_engine='+str(engine),'--extern','serde_json='+str(serde),source,'-o',binary])
call('read',[binary,document,manifest,out/'reader-receipt.json'])
receipt=read(out/'reader-receipt.json');original=read(raw)
assert len(receipt['frames'])==len(original['frames'])==5
checks=0
for frame,actual in zip(original['frames'],receipt['frames'],strict=True):
    assert actual['configured']['triangles']>0
    assert actual['configured']['frame_name']==frame['name']
    assert 'mc_epoch' not in actual['configured']
    for entry in actual['points']:
        expected=frame['state']['material_states'][entry['id']];value=entry['value']
        for key in ['effective_stress_compression_positive_kpa','plastic_strain_compression_positive','void_ratio','preconsolidation_pressure_kpa']:
            assert value[key]==expected[key],(frame['name'],entry['id'],key)
            checks+=1
    for entry in actual['nodes']:
        assert entry['value']['stored_node_fields']['total_displacement_from_seed_m']==frame['state']['displacement_m'][entry['id']]
        checks+=1
assert before=={q:sha(Path(q)) for q in before}
write(out/'terminal.json',dict(passed=True,frames=5,raw_source_bytes=raw.stat().st_size,
    representative_raw_field_checks=checks,source_and_document_input_identity_verified=True,
    original_source_unchanged=True,reader_sha256=sha(binary),native_identity=identity,
    original_and_new_result_families_remain_distinct=True,
    all_GP_full_field_accuracy_and_GPU_acceptance_separate=True,p6_complete=False))
print('PASS actual ordinary MCC all5 Native67-v4 reader scenes and representative raw fields')
