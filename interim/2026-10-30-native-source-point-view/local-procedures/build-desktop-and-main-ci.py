"""Build a bound desktop binary, then run the unchanged normal main-policy CI."""
from pathlib import Path
import datetime,hashlib,json,os,shutil,subprocess,sys,time
p=Path(__file__).resolve().parent;root=p.parents[1];native=root/'runs/native-typed-source-integration-20261029-v1/native';out=p/'desktop-and-main-ci-v1';out.mkdir()
read=lambda q:json.loads(q.read_text('utf-8-sig'))
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
def dump(q,v):
 with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2)
assert read(p/'native-quality-v2/terminal.json')['passed'] and read(p/'public-legacy-reader-v1/terminal.json')['passed']
assert read(p/'public8061-reader-v1/compiled-probe-binding.json')['passed']
sources=read(p/'native-quality-v2/source73.json')['files_sha256'];dump(out/'source73.json',read(p/'native-quality-v2/source73.json'))
dump(out/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
env=dict(os.environ,CARGO_TARGET_DIR=str(native/'target'),CARGO_BUILD_JOBS='2',SOLIDCRAFT_TEST_JOBS='1',RUST_TEST_THREADS='4',PYTHONUTF8='1')
def call(label,args):
 if (p/'pause-requested').exists():raise SystemExit(2)
 assert all(sha(native/name)==h for name,h in sources.items());command=['C:/Users/link_/.cargo/bin/cargo.exe',*args];start=time.monotonic()
 with (out/(label+'.log')).open('x',encoding='utf-8') as log:
  child=subprocess.Popen(command,cwd=native,env=env,stdout=log,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
  dump(out/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()));code=child.wait()
 assert all(sha(native/name)==h for name,h in sources.items())
 dump(out/(label+'-terminal.json'),dict(exit_code=code,elapsed_seconds=time.monotonic()-start,sources_unchanged=True));print(json.dumps(dict(step=label,exit_code=code)),flush=True);assert code==0
call('desktop-build',['build','--locked','--profile','ci','-p','solidcraft','--bin','solidcraft','--message-format=json-render-diagnostics'])
artifacts=[]
for line in (out/'desktop-build.log').read_text('utf-8').splitlines():
 try:item=json.loads(line)
 except json.JSONDecodeError:continue
 if item.get('reason')=='compiler-artifact' and item['target']['name']=='solidcraft' and item.get('executable'):artifacts.append(item['executable'])
assert len(set(artifacts))==1;exe=Path(artifacts[-1]);preserved=out/'solidcraft-source73.exe';shutil.copyfile(exe,preserved);assert sha(exe)==sha(preserved)
dump(out/'desktop-preserved.json',dict(passed=True,binary=str(preserved),sha256=sha(preserved),source73_manifest_sha256=sha(out/'source73.json'),
 experimental_CI_profile_not_release=True,actual_GUI_and_full_CI_not_yet_accepted=True,p6_complete=False))
call('main-ci',['xtask','ci'])
log=(out/'main-ci.log').read_text('utf-8')
assert 'ci: all gates passed' in log and 'property_tests=false' in log
assert 'skip kernel test boolean_props: nightly-only property tests as in PR/main CI' in log
dump(out/'terminal.json',dict(passed=True,source_files=73,full_normal_main_policy_CI_passed=True,
 original_failed_Native69_CI_and_nightly_counterexample_preserved=True,nightly_property_success_not_claimed=True,
 relevant_skip_or_ignore_lines=[line for line in log.splitlines() if any(word in line.lower() for word in ['skip','ignored','not found'])],
 official_release_and_changed_GUI_and_distribution_and_all_P6_not_yet_accepted=True,p6_complete=False))
print('PASS source73 full normal main-policy CI; nightly properties, release, GUI and all P6 remain separate',flush=True)
