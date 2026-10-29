"""Verify the real Git build identity and the corrected standard disclosure."""
from pathlib import Path
import datetime,hashlib,json,os,subprocess,sys,time,zipfile
p=Path(__file__).resolve().parent;root=p.parents[1];native=p/'native';out=p/'cache-kit-quality-v1';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8'))
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
def dump(q,v):
 with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2)
assert read(p/'tracked-native-preparation-v1/terminal.json')['passed']
prior=read(p/'quality-v5/source72.json');sources={name:sha(native/name) for name in prior['files_sha256']}
assert {name for name,h in prior['files_sha256'].items() if sources[name]!=h}=={'crates/ui-egui/src/soil_results.rs'}
dump(out/'source72.json',dict(base_commit=prior['base_commit'],files_sha256=sources,git_head=prior['base_commit'],working_patch_uncommitted=True,experimental_only=True))
with zipfile.ZipFile(out/'source72.zip','x',compression=zipfile.ZIP_DEFLATED) as z:
 for name in sources:z.write(native/name,name)
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
env=dict(os.environ,CARGO_TARGET_DIR=str(native/'target'),CARGO_BUILD_JOBS='2',RUST_TEST_THREADS='4')
cargo='C:/Users/link_/.cargo/bin/cargo.exe';records=[]
def check():assert all(sha(native/name)==h for name,h in sources.items())
def call(label,args):
 check()
 if (p/'pause-requested').exists():raise SystemExit(2)
 start=time.monotonic();command=[cargo,*args]
 with (out/(label+'.log')).open('x',encoding='utf-8') as log:
  child=subprocess.Popen(command,cwd=native,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
  dump(out/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
  code=child.wait()
 check();row=dict(step=label,exit_code=code,elapsed_seconds=time.monotonic()-start,sources_unchanged=True)
 dump(out/(label+'-terminal.json'),row);records.append(row);print(json.dumps(row),flush=True);assert code==0
call('format',['fmt','--all','--','--check'])
call('io-design',['test','--locked','--profile','ci','-p','solidcraft-io','--lib','design::tests','--','--test-threads=2'])
stamps=[]
for q in (native/'target/ci/build').glob('solidcraft-io-*/output'):
 for line in q.read_text('utf-8').splitlines():
  if line=='cargo:rustc-env=SOLIDCRAFT_SOURCE_COMMIT='+prior['base_commit']:stamps.append(dict(path=str(q),sha256=sha(q),value=prior['base_commit']))
assert stamps
dump(out/'known-build-identity.json',dict(passed=True,git_stamp=stamps,unknown_source_cache_guard_unchanged=True))
call('engine-file',['test','--locked','--profile','ci','-p','solidcraft-engine','--lib','cmd::file::tests','--','--test-threads=2'])
call('save-as',['test','--locked','--profile','ci','-p','solidcraft-engine','--test','save_as_cache','--','--test-threads=2'])
call('cli-model-cache',['test','--locked','--profile','ci','-p','solidcraft-cli','--test','model_cache','--','--test-threads=2'])
call('ui-kit-guard',['test','--locked','--profile','ci','-p','solidcraft-ui-egui','--lib','ui_kit_guard_tests::all_disclosures_use_the_one_kit_chevron_without_allowances','--','--test-threads=2'])
dump(out/'terminal.json',dict(passed=True,records=records,source_files=72,known_build_identity_verified=True,
 original_unknown_cache_guard_kept=True,standard_UI_disclosure_guard_passed=True,
 prior_Native69_failure_preserved=True,full_CI_and_actual_changed_GUI_still_required=True,p6_complete=False))
print('PASS focused model cache and standard disclosure gates; full Native CI remains separate',flush=True)
