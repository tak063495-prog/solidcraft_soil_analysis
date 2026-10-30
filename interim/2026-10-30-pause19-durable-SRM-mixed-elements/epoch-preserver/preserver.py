"""Preserve a registered accepted-origin SRM snapshot; never signal a process."""
from pathlib import Path
import hashlib,json,subprocess

def read(path,limit=256*1024):
    with path.open('rb') as f:data=f.read(limit+1)
    if len(data)>limit:raise ValueError('epoch SRM metadata exceeds budget')
    return json.loads(data)

def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def preserve_epoch_search_job(job,destination,solver,expected_solver_sha256,expected_job_id,expected_manifest_sha256):
    job,destination,solver=map(lambda q:Path(q).resolve(),[job,destination,solver])
    if destination.is_relative_to(job):raise ValueError('epoch SRM preservation must be outside source job')
    if solver!=job/'solver/soil-fem-3d.exe' or sha(solver)!=expected_solver_sha256:
        raise ValueError('registered epoch SRM executable differs')
    if sha(job/'job.json')!=expected_manifest_sha256:raise ValueError('registered epoch SRM manifest differs')
    source=read(job/'job.json');binding=source['binding']
    if source['schema']!='soil-fem-epoch-srm-search-job/1' or source['job_id']!=expected_job_id:
        raise ValueError('registered epoch SRM identity differs')
    if binding['identity']['solver_binary_sha256']!=expected_solver_sha256:
        raise ValueError('epoch SRM origin/executable binding differs')
    expected={'inputs/mesh.json','inputs/origin.json','inputs/specification.json','solver/soil-fem-3d.exe'}
    if binding['parent'] is not None:expected.add('inputs/checkpoint.json')
    if set(binding['files_sha256'])!=expected:raise ValueError('epoch SRM immutable inventory differs')
    def verify_original():
        if sha(job/'job.json')!=expected_manifest_sha256:raise ValueError('epoch SRM manifest changed during preservation')
        for name,h in binding['files_sha256'].items():
            if sha(job/name)!=h:raise ValueError('epoch SRM immutable input changed during preservation')
    verify_original()
    destination.mkdir(parents=True,exist_ok=False);frozen=destination/'job'
    with (destination/'prepare.log').open('x',encoding='utf-8') as log:
        process=subprocess.run([str(solver),'prepare-mc-epoch-srm-recovery',str(job),str(frozen)],
            stdout=log,stderr=subprocess.STDOUT,timeout=240,creationflags=subprocess.CREATE_NO_WINDOW)
    if process.returncode:
        tail=(destination/'prepare.log').read_text('utf-8',errors='replace')[-8192:]
        race=any(word in tail.lower() for word in ['changed','advanced'])
        raise ValueError('epoch SRM recovery preparation '+('changed/advanced during copy' if race else 'failed; inspect prepare.log'))
    copied=read(frozen/'job.json');child_binding=copied['binding'];parent=child_binding['parent']
    if copied['schema']!=source['schema'] or child_binding['identity']!=binding['identity'] or parent['job_id']!=expected_job_id:
        raise ValueError('epoch SRM fork changed origin/provenance')
    for name in ['specification_sha256','material_attempts_per_save']:
        if child_binding[name]!=binding[name]:raise ValueError('epoch SRM fork changed specification/save interval')
    if set(child_binding['files_sha256'])!=expected|{'inputs/checkpoint.json'}:
        raise ValueError('epoch SRM fork immutable inventory differs')
    for name,h in binding['files_sha256'].items():
        if name!='inputs/checkpoint.json' and child_binding['files_sha256'].get(name)!=h:
            raise ValueError('epoch SRM fork changed original immutable files')
    child_solver=frozen/'solver/soil-fem-3d.exe'
    if sha(child_solver)!=expected_solver_sha256:raise ValueError('epoch SRM fork executable differs')
    process=subprocess.run([str(child_solver),'verify-mc-epoch-srm-job',str(frozen)],capture_output=True,
        text=True,encoding='utf-8-sig',timeout=240,creationflags=subprocess.CREATE_NO_WINDOW)
    (destination/'verify.stdout').write_text(process.stdout,'utf-8');(destination/'verify.stderr').write_text(process.stderr,'utf-8')
    if process.returncode:raise ValueError('epoch SRM full typed context/checksum verification failed')
    audit=json.loads(process.stdout)
    if audit['schema']!='soil-fem-epoch-srm-search-snapshot-verification/1' or audit['job_id']!=copied['job_id']:
        raise ValueError('epoch SRM fork audit identity differs')
    if audit['generation']!=1 or audit['canonical_snapshot_sha256']!=parent['checkpoint_payload_sha256']:
        raise ValueError('epoch SRM fork changed complete accepted payload')
    verify_original()
    receipt=dict(passed=True,source_job=str(job),source_job_id=expected_job_id,source_manifest_sha256=expected_manifest_sha256,
        selected_generation=parent['generation'],selected_source_snapshot_sha256=parent['snapshot_sha256'],
        selected_canonical_snapshot_sha256=parent['checkpoint_payload_sha256'],frozen_job=str(frozen),
        full_typed_context_and_payload_checksum_audit=audit,solver_sha256=expected_solver_sha256,
        same_unreduced_completed_origin_and_search_and_save_interval_preserved=True,
        same_factor_accepted_history_and_completed_trials_preserved=True,
        saved_material_epoch=bool(audit['saved_material_epoch']),search_complete=bool(audit['complete']),
        initial_prepared_snapshot_not_called_an_equilibrated_state=True,source_immutable_inputs_unchanged=True,
        no_process_signalled=True,unaccepted_trial_not_claimed_saved=True,
        equilibrium_or_safety_not_inferred_from_snapshot_checksum=True,p6_complete=False)
    with (destination/'receipt.json').open('x',encoding='utf-8') as f:json.dump(receipt,f,indent=2,allow_nan=False)
    return receipt
