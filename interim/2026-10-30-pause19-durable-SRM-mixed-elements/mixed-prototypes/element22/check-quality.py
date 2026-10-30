"""Build only the small isolated split probe against the closed qualified662 library."""
from pathlib import Path
import datetime,hashlib,json,os,subprocess,time
p=Path(__file__).resolve().parent;root=p.parents[1];base=p/'mixed-element-prototype-v1';out=base/'quality-v1';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
assert not (p/'pause-requested').exists()
parent=p/'epoch-SRM-control-v2';accepted=read(parent/'quality-v1/terminal.json')
assert accepted['passed'] and (accepted['Rust_tests'],accepted['Python_tests'])==(370,219)
sources=read(parent/'quality-v1/source-snapshot.json');assert all(sha(parent/'runtime'/n)==h for n,h in sources.items())
root_sources=read(p/'epoch-SRM-job-v1/Root657-source-snapshot.json');assert all(sha(root/n)==h for n,h in root_sources.items())
current=read(p/'epoch-SRM-job-v1/quality-v2/source-snapshot.json');assert all(sha(p/'epoch-SRM-job-v1/runtime'/n)==h for n,h in current.items())
dependencies=parent/'target/release/deps';libraries=list(dependencies.glob('libsoil_fem_3d-*.rlib'));assert len(libraries)==1
library=libraries[0];library_sha=sha(library)
toolchain=Path('C:/Users/link_/.cargo/bin');records=[]
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=__import__('sys')._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
def call(label,args):
    start=time.monotonic()
    with (out/(label+'.log')).open('x',encoding='utf-8') as log:
        child=subprocess.Popen(list(map(str,args)),cwd=base,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        dump(out/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=str(args[0]),arguments=list(map(str,args[1:])),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code=child.wait()
    row=dict(phase=label,exit_code=code,elapsed_seconds=time.monotonic()-start);records.append(row);dump(out/(label+'-terminal.json'),row)
    print(json.dumps(row),flush=True)
    if code:
        dump(out/'terminal.json',dict(passed=False,failed_phase=label,records=records,Root_and_frozen_candidates_unchanged=True,p6_complete=False));raise SystemExit(code)
call('format',[toolchain/'rustfmt.exe','--edition','2024',base/'mc_split.rs',base/'tests.rs',base/'soil_split.rs'])
frozen={str(q):sha(q) for q in [base/'mc_split.rs',base/'tests.rs',base/'soil_split.rs',base/'typed_dispatch_tests.rs',base/'mixed_element.rs',base/'mixed_element_tests.rs',Path(__file__)]};dump(out/'source-binding.json',frozen)
call('compile',[toolchain/'rustc.exe','--edition=2024','--test','-C','opt-level=2','-L',f'dependency={dependencies}','--extern',f'soil_fem_3d={library}',base/'tests.rs','-o',out/'MC-split-tests.exe'])
call('tests',[out/'MC-split-tests.exe','--test-threads=1','--nocapture'])
assert 'test result: ok. 22 passed; 0 failed; 0 ignored;' in (out/'tests.log').read_text('utf-8')
assert sha(library)==library_sha and all(sha(Path(n))==h for n,h in frozen.items())
assert all(sha(root/n)==h for n,h in root_sources.items())
assert all(sha(p/'epoch-SRM-job-v1/runtime'/n)==h for n,h in current.items())
dump(out/'terminal.json',dict(passed=True,prototype_only=True,tests=22,parent_qualified_Rust370Python219=True,
    parent_manifest=str(parent/'quality-v1/source-snapshot.json'),parent_manifest_sha256=sha(parent/'quality-v1/source-snapshot.json'),
    actual_parent_library=str(library),actual_parent_library_sha256=library_sha,source_binding=frozen,records=records,
    typed_MC_history_and_existing_fixed_partition_nonassociated_map_retained=True,
    independent_elastic_moduli_compliance_directional_derivatives_work_rotation_and_failure_checks=True,
    isolated_element_Jacobian_and_nullspaces_checked_but_global_controller_and_SRM_certificate_NOT_implemented=True,
    all9_return_types_flow_variants_and_explicit_transition_split_covered=True,
    production_or_frozen664_sources_not_changed=True,original2D5percent_and_P6_not_established=True,p6_complete=False))
print('PASS isolated mixed element22 tests; no mixed FEM/5percent/P6 acceptance claimed',flush=True)
