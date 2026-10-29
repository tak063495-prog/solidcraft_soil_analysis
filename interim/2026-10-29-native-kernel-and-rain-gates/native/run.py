"""Preserve the failed full CI and isolate its exact case without a CI rerun."""
from pathlib import Path
import datetime,hashlib,json,os,shutil,subprocess,sys,time,zipfile
p=Path(__file__).resolve().parent;root=p.parents[1]
prior=root/'runs/native-snapshot-recovery-development-20261029-v1';checkout=prior/'validation-checkout';target=root/'runs/native-large-source-stream-20261029-v1/checkout/target'
read=lambda q:json.loads(q.read_text('utf-8-sig'))
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
def write(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2)
write(p/'driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
terminal=read(prior/'quality-v3/ci-terminal.json');assert terminal['exit_code']==1 and terminal['sources_unchanged']
frozen=read(prior/'source65-v3.json')['files_sha256'];assert len(frozen)==65 and all(sha(checkout/n)==h for n,h in frozen.items())
assert not any(n.startswith('crates/kernel/') for n in frozen)
base='edf57cec293dcea301201f1d5173b33e16351466'
kernel_changes=subprocess.check_output(['git','diff','--name-only',base,'--','crates/kernel/src'],cwd=checkout,encoding='utf-8').strip();assert not kernel_changes
for name in ['ci.log','ci-terminal.json']:shutil.copy2(prior/'quality-v3'/name,p/name)
for name in ['source65-v3.json','source65-v3.zip']:shutil.copy2(prior/name,p/name)
reg='crates/kernel/tests/boolean_props.proptest-regressions';assert 'fe2e0b403aaded8bea37c9f5429ac422c96c7a1d9c012771b25d5a3a093e94ff' in (checkout/reg).read_text('utf-8')
shutil.copy2(checkout/reg,p/'preserved-boolean-proptest-regressions.txt')
names=subprocess.check_output(['git','ls-files','--modified','--others','--exclude-standard'],cwd=checkout,encoding='utf-8').splitlines()
assert set(names)==set(frozen)|{reg},set(names)-set(frozen)
failed_source={name:sha(checkout/name) for name in sorted(names)}
write(p/'failed-source66.json',dict(base=base,files_sha256=failed_source,production_frozen65_unchanged=True,extra_test_generated_seed_file=reg))
with zipfile.ZipFile(p/'failed-source66.zip','x',compression=zipfile.ZIP_DEFLATED) as z:
    for name,h in failed_source.items():assert sha(checkout/name)==h;z.write(checkout/name,name)
write(p/'failure.json',dict(native_full_ci_passed=False,ci_elapsed_seconds=terminal['elapsed_seconds'],engine2038_tests_passed=True,
    failing_binary='kernel test boolean_props',test='identity_on_grid',kernel_production_sources_unchanged_from_official_base=True,
    input_mm=dict(cylinder=dict(base=[25.,0.,0.],axis=[0.,0.,1.],radius=15.,height=20.),sphere=dict(centre=[27.5,0.,0.],radius=12.5)),
    seed='fe2e0b403aaded8bea37c9f5429ac422c96c7a1d9c012771b25d5a3a093e94ff',
    expected_intersection='Upper hemisphere: sphere radial offset2.5+radius12.5 equals cylinder radius15, and base plane bisects sphere; top20 exceeds sphere extent12.5.',
    error='Intersect no result; exact fallback supports planar faces but cylinder/sphere have curved faces; nudged retry also failed.',
    production_code_not_modified=True,full_ci_not_rerun=True,fresh_install_release_gpu_gui_recovery_waiters_aborted_before_start=True,p6_complete=False))
env=dict(os.environ,CARGO_TARGET_DIR=str(target),CARGO_BUILD_JOBS='2',SOLIDCRAFT_TEST_JOBS='1',RUST_TEST_THREADS='4')
def call(label,args,stdout_path=None):
    if (p/'pause-requested').exists():raise SystemExit(2)
    command=list(map(str,args));start=time.monotonic()
    output=stdout_path or p/(label+'.stdout')
    with output.open('x',encoding='utf-8') as out,(p/(label+'.stderr')).open('x',encoding='utf-8') as err:
        child=subprocess.Popen(command,cwd=checkout,env=env,stdout=out,stderr=err,creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        write(p/(label+'-child.json'),dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()));code=child.wait()
    row=dict(step=label,exit_code=code,elapsed_seconds=time.monotonic()-start);write(p/(label+'-terminal.json'),row);print(json.dumps(row),flush=True);assert code==0
# The exact xtask build-only command reports the frozen dependency artifacts;
# no test executable or random property generator is run here.
call('artifact-discovery',['C:/Users/link_/.cargo/bin/cargo.exe','test','--workspace','--profile','ci','--no-run','--message-format=json-render-diagnostics'])
artifacts=[json.loads(s) for s in (p/'artifact-discovery.stdout').read_text('utf-8').splitlines() if s.startswith('{')]
libs={Path(n) for a in artifacts if a.get('reason')=='compiler-artifact' and a.get('target',{}).get('name')=='solidcraft_kernel' for n in a['filenames'] if n.endswith('.rlib')}
assert len(libs)==1,[str(q) for q in libs]
lib=libs.pop();assert lib.is_file()
write(p/'dependency-identity.json',dict(kernel_rlib=str(lib),sha256=sha(lib),base_kernel_sources_unchanged=True,no_test_run_for_discovery=True))
call('compile-diagnostic',['C:/Users/link_/.cargo/bin/rustc.exe','--edition=2024',p/'reproducer.rs','-L','dependency='+str(lib.parent),'--extern','solidcraft_kernel='+str(lib),'-o',p/'reproducer.exe'])
call('run-diagnostic',[p/'reproducer.exe'])
rows=[json.loads(s) for s in (p/'run-diagnostic.stdout').read_text('utf-8').splitlines() if s.strip()];assert len(rows)==9
target_rows=[r for r in rows if r['offset_mm']==2.5]
reproduced=any(r['operation']=='Intersect' and (r['status']!='ok' or abs(r['volume_mm3']-r['expected_volume_mm3'])>.002*r['expected_volume_mm3']) for r in target_rows)
assert all(sha(checkout/n)==h for n,h in failed_source.items())
write(p/'diagnostic-terminal.json',dict(completed=True,kernel_failure_reproduced=reproduced,rows=rows,
    no_full_ci_rerun=True,no_production_code_changes=True,all_frozen65_and_new_seed_bytes_unchanged=True,native_release_not_accepted=True,p6_complete=False))
print('Kernel exact-case diagnostic complete; reproduced:',reproduced,flush=True)
