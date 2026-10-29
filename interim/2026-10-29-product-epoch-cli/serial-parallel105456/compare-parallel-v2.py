"""Sequential full-field hashes; no numerical rerun and no isolated speed claim."""
from pathlib import Path
import datetime,gc,hashlib,json,os,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1]
parallel=root/'runs/mcc-parallel-100k-durable-20261029-v1/n26';serial=p/'n26'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
write=lambda q,v:q.write_text(json.dumps(v,indent=2,allow_nan=False),'utf-8')
def sha(q):
    h=hashlib.sha256()
    with q.open('rb') as f:
        for b in iter(lambda:f.read(4<<20),b''):h.update(b)
    return h.hexdigest()
assert read(serial/'terminal.json')['passed'] and read(parallel/'terminal.json')['passed']
write(p/'comparison-v2-driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,
    script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
inputs={}
for n in ['mesh.json','ground.json','stages.json']:
    assert sha(serial/n)==sha(parallel/n);inputs[n]=sha(serial/n)
groups=['frames','last_committed','reference_displacement_m','integration_points',
    'excavation_audits','transition_audits','geostatic','attempts']
encoder=json.JSONEncoder(ensure_ascii=False,sort_keys=True,separators=(',',':'),allow_nan=False)
def strip_elapsed(value):
    if isinstance(value,dict):
        value.pop('elapsed_seconds',None)
        for v in value.values():strip_elapsed(v)
    elif isinstance(value,list):
        for v in value:strip_elapsed(v)
def hash_value(v):
    h=hashlib.sha256()
    for s in encoder.iterencode(v):h.update(s.encode('utf-8'))
    return h.hexdigest()
records=[]
for label,folder in [('serial',serial),('parallel',parallel)]:
    if (p/'pause-requested').exists():raise SystemExit(2)
    source=folder/'solver-job/result.json';source_hash=sha(source);start=time.monotonic()
    r=read(source);assert r['complete'] and len(r['frames'])==3
    solver=r['provenance']['solver_binary_sha256'];strip_elapsed(r)
    hashes={};defaults=[]
    for group in groups:
        if (p/'pause-requested').exists():raise SystemExit(2)
        # Run::transition_audits declares serde(default, skip_serializing_if=Vec::is_empty).
        # Both old formats legitimately omit that empty vector. No other missing
        # group is accepted by this comparison.
        if group=='transition_audits' and group not in r:
            value=[];defaults.append(group)
        else:value=r.pop(group)
        hashes[group]=hash_value(value);del value;gc.collect()
    del r;gc.collect();assert sha(source)==source_hash
    record=dict(label=label,result_sha256=source_hash,solver_sha256=solver,all8_groups_sha256=hashes,
        elapsed_seconds=time.monotonic()-start,declared_serde_empty_vector_defaults=defaults)
    write(p/(label+'-full-field-hashes.json'),record);records.append(record)
    print('Hashed all8 groups '+label,flush=True)
assert records[0]['all8_groups_sha256']==records[1]['all8_groups_sha256']
write(p/'serial-parallel-comparison.json',dict(passed=True,elements=105456,
    all8_groups_all_nodes_gp_and_attempts_exact_except_elapsed=True,
    different_solver_identities_explicitly_recorded=True,input_sha256=inputs,records=records,
    both_independent_analytical_audits_passed=True,concurrent_other_jobs=True,
    isolated_throughput_benchmark=False,no_speedup_ratio_claimed=True,
    scope='Uniform hydrostatic nonlinear MCC patch; no general CAD, excavation, water or structures.',p6_complete=False))
print('PASS105456 serial/parallel all8 full-field groups exact except elapsed; no isolated speed ratio',flush=True)
