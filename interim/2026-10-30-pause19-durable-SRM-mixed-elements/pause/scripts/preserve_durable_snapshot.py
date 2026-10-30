"""Deep-validate an immutable accepted snapshot before stopping a long Newton trial.

Only a trusted caller supplies a job's already-verified frozen solver. This never
signals a PID or modifies the original job, and never edits checkpoint numbers.
"""
from pathlib import Path
import hashlib,json,subprocess

def read(path,limit):
    with path.open('rb') as f:raw=f.read(limit+1)
    if len(raw)>limit:raise ValueError('snapshot metadata exceeds budget')
    return json.loads(raw)

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for raw in iter(lambda:f.read(1<<20),b''):h.update(raw)
    return h.hexdigest()

def copy(source,target,digest):
    if source.stat().st_size>8<<30 or sha(source)!=digest:raise ValueError('snapshot source changed or exceeds budget')
    target.parent.mkdir(parents=True,exist_ok=True)
    with source.open('rb') as src,target.open('xb') as dst:
        total=0
        for raw in iter(lambda:src.read(1<<20),b''):
            total+=len(raw)
            if total>8<<30:raise ValueError('snapshot grew beyond budget')
            dst.write(raw)
    if sha(source)!=digest or sha(target)!=digest:raise ValueError('snapshot changed while copying')

def preserve_snapshot(job,destination,solver,expected_solver_sha256):
    job,solver,destination=Path(job).resolve(),Path(solver).resolve(),Path(destination).resolve()
    if destination.is_relative_to(job):raise ValueError('snapshot preservation must be outside original job')
    if solver!=job/'solver/soil-fem-3d.exe' or sha(solver)!=expected_solver_sha256:
        raise ValueError('trusted frozen solver differs')
    manifest=read(job/'job.json',256<<10);manifest_sha=sha(job/'job.json')
    pointer=read(job/'recovery.json',8192);pointer_sha=sha(job/'recovery.json')
    binding=manifest.get('binding',{});files=binding.get('files_sha256',{})
    required={'inputs/mesh.json','inputs/ground.json','inputs/case.json','solver/soil-fem-3d.exe'}
    interval=binding.get('checkpoint_interval');generation=pointer.get('generation')
    if (manifest.get('schema')!='soil-fem-stages-job/1' or type(interval) is not int or not 1<=interval<=4096
        or not isinstance(files,dict) or set(files) not in (required,required|{'inputs/checkpoint.json'})
        or files.get('solver/soil-fem-3d.exe')!=expected_solver_sha256
        or pointer.get('schema')!='soil-fem-stages-recovery/1' or pointer.get('job_id')!=manifest.get('job_id')
        or type(generation) is not int or not 1<=generation<=2**64-1
        or pointer.get('checkpoint')!='recovery-checkpoint-'+str(generation%2)+'.json'):
        raise ValueError('periodic accepted snapshot binding differs')
    for digest in [*files.values(),pointer.get('checkpoint_sha256'),pointer.get('run_sha256')]:
        if not isinstance(digest,str) or len(digest)!=64 or any(c not in '0123456789abcdef' for c in digest):
            raise ValueError('invalid snapshot digest')
    if read(job/'recovery.json',8192)!=pointer or sha(job/'recovery.json')!=pointer_sha:
        raise ValueError('snapshot pointer advanced during selection')
    destination.mkdir(parents=True,exist_ok=False)
    source=destination/'immutable-source';source.mkdir()
    for name,digest in files.items():copy(job/name,source/name,digest)
    copy(job/'job.json',source/'job.json',manifest_sha)
    copy(job/pointer['checkpoint'],source/pointer['checkpoint'],pointer['checkpoint_sha256'])
    copy(job/'recovery.json',source/'recovery.json',pointer_sha)
    if read(source/'recovery.json',8192)!=pointer or sha(job/'recovery.json')!=pointer_sha:
        raise ValueError('snapshot pointer advanced while preserving')
    recovered=destination/'prepared-recovery-job'
    # The caller chose this trusted executable by its actual owned process and
    # checksum, never by an executable path supplied in imported job metadata.
    with (destination/'deep-prepare.log').open('xb') as log:
        result=subprocess.run([str(solver),'prepare-stages-recovery',str(source),str(recovered)],
            stdout=log,stderr=subprocess.STDOUT,timeout=60,
            **(dict(creationflags=subprocess.CREATE_NO_WINDOW) if __import__('os').name=='nt' else {}))
    if result.returncode:raise ValueError('Rust checkpoint/material-history validation failed')
    prepared=read(recovered/'job.json',256<<10)
    if (prepared.get('binding',{}).get('parent_job_id')!=manifest['job_id']
        or prepared['binding'].get('identity')!=binding.get('identity')
        or prepared['binding']['files_sha256'].get('inputs/checkpoint.json')!=pointer['checkpoint_sha256']
        or sha(recovered/'inputs/checkpoint.json')!=pointer['checkpoint_sha256']
        or sha(recovered/'solver/soil-fem-3d.exe')!=expected_solver_sha256):
        raise ValueError('prepared immutable recovery job differs from selected accepted snapshot')
    receipt=dict(passed=True,original_job=str(job),accepted_generation=generation,
        original_job_not_modified=True,no_pid_signalled=True,unaccepted_newton_trial_not_preserved=True,
        immutable_source=str(source),prepared_recovery_job=str(recovered),
        pointer_sha256=pointer_sha,checkpoint_sha256=pointer['checkpoint_sha256'],
        solver_sha256=expected_solver_sha256,deep_rust_history_validation_passed=True,p6_complete=False)
    with (destination/'preservation.json').open('x') as f:json.dump(receipt,f,indent=2)
    return receipt
