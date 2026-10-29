"""Adopt only the accepted endpoint replay candidate, then actual Root controls.

No CSR candidate is combined here. Full 277/171 QA is bound by byte-identical
non-Markdown sources; actual Root build, MC CLI and MCC fields are checked.
"""
from pathlib import Path
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

p = Path(__file__).resolve().parent
root = p.parents[1]
candidate = root/'runs/mc-stage-gate-development-20261029-v1'
read = lambda q: json.loads(q.read_text('utf-8-sig'))
sha = lambda q: hashlib.sha256(q.read_bytes()).hexdigest()
files = ['src/replay.rs','src/main.rs','src/stages.rs','tests/replay.rs']
env = dict(os.environ,CARGO_BUILD_JOBS='2',RUST_TEST_THREADS='2',PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1')


def write(path,value):
    with path.open('x',encoding='utf-8') as f:
        json.dump(value,f,indent=2,allow_nan=False)


def wait(path):
    while not path.exists():
        if (p/'pause-requested').exists():
            raise SystemExit(2)
        for folder in [candidate/'quality-v2',candidate/'cli-controls',candidate/'cad-controls']:
            for q in folder.rglob('*-terminal.json'):
                r = read(q)
                if 'exit_code' in r:
                    assert r['exit_code']==r.get('expected_exit_code',0),(str(q),r)
        time.sleep(5)
    result = read(path)
    assert result['passed'],str(path)
    return result


def call(label,args,folder=p):
    if (p/'pause-requested').exists():
        raise SystemExit(2)
    command = list(map(str,args))
    started = time.monotonic()
    with (folder/(label+'.log')).open('x',encoding='utf-8') as log:
        child = subprocess.Popen(command,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        write(folder/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code = child.wait()
    row = dict(step=label,exit_code=code,elapsed_seconds=time.monotonic()-started)
    write(folder/(label+'-terminal.json'),row)
    print(json.dumps(row),flush=True)
    assert code==0,row


def stable(v, remove_binary_identity=False):
    if isinstance(v,dict):
        return {k:stable(x,remove_binary_identity) for k,x in v.items() if k!='elapsed_seconds' and not(remove_binary_identity and k=='solver_binary_sha256')}
    if isinstance(v,list):
        return [stable(x,remove_binary_identity) for x in v]
    return v


write(p/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
assert wait(candidate/'terminal.json')['rust_tests']==277
assert wait(candidate/'cad-controls/terminal.json')['independent_3d_elastic_history_and_duffy_integration_and_svd_for_all33_cases']
quality = read(candidate/'quality-v2/terminal.json')
assert quality['passed'] and quality['rust_tests']==277 and quality['python_tests']==171
accepted = read(candidate/'quality-v2/source-snapshot.json')
baseline = read(candidate/'source-snapshot.json')
assert len(accepted)==len(baseline)==611
assert {n for n in accepted if accepted[n]!=baseline[n]}==set(files)
assert all(sha(candidate/'runtime'/n)==h for n,h in accepted.items())
sys.path.insert(0,str(root/'tools'))
from freeze_sources import freeze,manifest
before = p/'before'
before.mkdir()
freeze(before)
initial = manifest()
assert set(initial)==set(accepted)
assert all(initial[n]==baseline[n] for n in initial if not n.endswith('.md'))
old_solver = root/'target/release/soil-fem-3d.exe'
assert sha(old_solver)=='d55ffe42c3b61e0891fec513489c9f5881830733ac9e003a01a0cf69f495905d'
shutil.copy2(old_solver,before/'soil-fem-3d.exe')
for name in files:
    backup = before/'changed-files'/name
    backup.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(root/name,backup)
    shutil.copy2(candidate/'runtime'/name,root/name)
current = manifest()
assert all(current[n]==accepted[n] for n in current if not n.endswith('.md'))
assert {n for n in initial if initial[n]!=current[n]}==set(files)
write(p/'adoption.json',dict(passed=True,adopted_files=files,full277_171_accepted_in_private_runtime=True,
    all_root_non_markdown_sources_byte_exact_to_accepted_candidate=True,root_independent_build_and_actual_controls_pending=True,
    csr_candidate_not_combined=True,private_quality_sha256=sha(candidate/'quality-v2/terminal.json'),cad_controls_sha256=sha(candidate/'cad-controls/terminal.json')))
cargo = 'C:/Users/link_/.cargo/bin/cargo.exe'
call('format-check',[cargo,'fmt','--all','--','--check'])
call('build',[cargo,'build','--locked','--release'])
solver = p/'solver/soil-fem-3d.exe'
solver.parent.mkdir()
shutil.copy2(old_solver,solver)
solver_sha = sha(solver)
assert manifest()==current
frozen = p/'source'
frozen.mkdir()
freeze(frozen)
runtime = p/'runtime'
import zipfile
runtime.mkdir()
with zipfile.ZipFile(frozen/'source-snapshot.zip') as archive:
    archive.extractall(runtime)
cli = p/'cli-controls'
cli.mkdir()
for name in ['mesh.json','ground.json','case.json','spec.json']:
    shutil.copy2(candidate/'cli-controls'/name,cli/name)
paths = [cli/name for name in ['mesh.json','ground.json','case.json','spec.json']]
call('full',[solver,'replay-stages',*paths,cli/'full.json'],cli)
assert stable(read(cli/'full.json'),True)==stable(read(candidate/'cli-controls/full.json'),True)
case = read(cli/'case.json')
for count in [0,1,len(case['stages'])]:
    call('gate'+str(count),[solver,'replay-through-stage',*paths,count,cli/f'gate{count}.json'],cli)
    assert stable(read(cli/f'gate{count}.json'),True)==stable(read(candidate/f'cli-controls/gate{count}.json'),True)
baseline_analysis = root/'runs/mcc-parallel-material-development-20261029-v1/controls/native274/workflow/analysis'
analysis = p/'control/workflow/analysis'
(analysis/'ground').mkdir(parents=True)
for name in ['ground/mesh.json','ground/ground.json','stages.json']:
    shutil.copy2(baseline_analysis/name,analysis/name)
    assert sha(baseline_analysis/name)==sha(analysis/name)
job = analysis/'solver-job'
call('prepare-mcc',[solver,'prepare-durable-stages-job',analysis/'ground/mesh.json',analysis/'ground/ground.json',analysis/'stages.json','4',job])
call('solve-mcc',[job/'solver/soil-fem-3d.exe','run-stages-job',job])
shutil.copy2(job/'result.json',analysis/'stage-results.json')
call('export-mcc',[sys.executable,runtime/'tools/export_stages.py',analysis/'ground/mesh.json',analysis/'stage-results.json',analysis/'results'])
call('verify-mcc',[sys.executable,runtime/'tools/verify_stage_package.py',analysis/'ground/mesh.json',analysis/'ground/ground.json',analysis/'stage-results.json',analysis/'results','--solver',solver])
expected,actual = read(baseline_analysis/'stage-results.json'),read(analysis/'stage-results.json')
groups = ['frames','last_committed','reference_displacement_m','integration_points','excavation_audits','transition_audits','geostatic','attempts']
assert actual['complete'] and len(actual['frames'])==5
for key in groups:
    assert stable(expected[key])==stable(actual[key]),key
assert manifest()==current and sha(solver)==solver_sha
write(p/'terminal.json',dict(passed=True,adopted=True,rust_tests=277,python_tests=171,
    full_quality_reused_only_with_byte_exact_root_non_markdown_sources=True,root_build_passed=True,
    root_actual_mc_full_and0_first_final_cli_outputs_exact_except_binary_hash_and_elapsed=True,
    root_actual_mcc274_all5_stages_and_all8_physical_groups_exact=True,root_mcc_export_and_mechanical_audits_passed=True,
    solver_sha256=solver_sha,adopted_files=files,csr_candidate_not_combined=True,
    fs_search_and_physical_classification_not_implemented=True,github_publication_and_new_windows_ci_pending=True,p6_complete=False))
print('PASS adopted endpoint replay277/171, exact source binding and actual Root MC/MCC controls; publication/new CI pending',flush=True)
