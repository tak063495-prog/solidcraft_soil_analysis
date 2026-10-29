"""Replay every public comparison row before creating the review archive."""
from pathlib import Path
import datetime,hashlib,json,os,subprocess,sys,time,zipfile
p=Path(__file__).resolve().parent;root=p.parents[1];data=p/'review-package'
read=lambda q:json.loads(q.read_text('utf-8-sig'));sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
with (p/'replay-driver-child.json').open('x') as f:json.dump(dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f)
metadata=read(data/'metadata.json');assert metadata['prepared'] and not metadata['full_gp_histories_in_public_subset']
assert all(sha(data/n)==h for n,h in metadata['files_sha256'].items())
if (p/'pause-requested').exists():raise SystemExit(2)
command=[sys.executable,'-X','utf8',str(data/'review-comparison.py')];started=time.monotonic()
with (p/'replay.log').open('x',encoding='utf-8') as log:
    child=subprocess.Popen(command,cwd=root,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1'),stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
    with (p/'replay-child.json').open('x') as f:json.dump(dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),f)
    code=child.wait()
with (p/'replay-terminal.json').open('x') as f:json.dump(dict(exit_code=code,elapsed_seconds=time.monotonic()-started),f,indent=2)
assert code==0 and read(data/'comparison-recomputed.json')['passed']
assert all(sha(data/n)==h for n,h in metadata['files_sha256'].items())
archive=p/'fixed-wall-comparison-review.zip';inventory={**metadata['files_sha256'],'metadata.json':sha(data/'metadata.json')}
with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED) as z:
    for n,h in inventory.items():assert sha(data/n)==h;z.write(data/n,n)
with zipfile.ZipFile(archive) as z:
    assert set(z.namelist())==set(inventory)
    for n,h in inventory.items():assert hashlib.sha256(z.read(n)).hexdigest()==h
with (p/'terminal.json').open('x') as f:json.dump(dict(passed=True,all5_nodal_and_wall_rows_and_geometry_metadata_exact=True,all_published_files_hash_verified=True,
    original_full_results_bytes=metadata['full_results_bytes'],public_full_gp_histories=False,full_gp_audit_claims_only_recorded_separately=True,
    archive_bytes=archive.stat().st_size,archive_sha256=sha(archive),files=len(inventory),frozen_solver_and_both_independent_audit_sources_included=True,
    a_new_independent_nonlinear_solution_not_claimed=True,spatial_convergence_accepted=False,new_pc_acceptance=False,p6_complete=False),f,indent=2)
print('PASS portable5-frame comparison replay and immutable archive inventory; full GP history scope explicit',flush=True)
