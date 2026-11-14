import concurrent.futures,hashlib,json,os,pathlib,subprocess,sys,time
root=pathlib.Path(__file__).resolve().parents[1];out=root/'work/b02cause-validation/python-full';assert out.is_dir()
runtime=root/'work/b02cause-powershell'
env=dict(os.environ,PATH=str(runtime/'runtime')+os.pathsep+os.environ['PATH'],XDG_CACHE_HOME=str(runtime/'cache'),XDG_CONFIG_HOME=str(runtime/'config'),XDG_DATA_HOME=str(runtime/'data'))
paths=[p for p in subprocess.check_output(['git','ls-files'],cwd=root,text=True).splitlines() if pathlib.Path(p).suffix in ('.rs','.py','.ps1','.toml') or pathlib.Path(p).name in ('Cargo.lock','rust-toolchain.toml')]
before={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in sorted(paths)};(out/'source-binding-before-core.json').write_text(json.dumps(before,indent=2)+'\n')
groups=['python-core-0','python-core-1'];start=time.monotonic()
def run(g):
 cmd=[sys.executable,'tools/run_ci_python_evidence_shards.py','--group',g,'--output',str(out/g)]
 t=time.monotonic()
 with (out/(g+'-driver.log')).open('x',encoding='utf-8') as f:
  p=subprocess.Popen(cmd,cwd=root,env=env,stdout=f,stderr=subprocess.STDOUT);code=p.wait()
 record=dict(group=g,command=cmd,pid=p.pid,exit_code=code,elapsed_seconds=time.monotonic()-t)
 (out/(g+'-process.json')).write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record),flush=True);return record
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:records=list(pool.map(run,groups))
after={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in before};(out/'source-binding-after.json').write_text(json.dumps(after,indent=2)+'\n')
all_groups=['python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause']
while not all((out/(g+'-process.json')).exists() for g in all_groups):time.sleep(5)
records=[json.loads((out/(g+'-process.json')).read_text()) for g in all_groups]
early=json.loads((out/'source-binding-before.json').read_text())
changes=sorted(p for p in early if early[p]!=after[p])
allowed=['interim/2026-12-06-p6-01-status-reconciliation/audit_status.py','interim/2026-12-10-ci-p601-inventory-fix/inventory.toml','interim/2026-12-24-b02-offset-composite-cause-diagnosis/b02cause_evidence.py','interim/2026-12-24-b02-offset-composite-cause-diagnosis/b02cause_report.py','tests/test_b02_offset_composite_cause_diagnosis.py']
result=dict(groups=records,elapsed_seconds_for_core_phase=time.monotonic()-start,all_groups_attempted=True,source_count=len(before),sources_unchanged=not changes,core_sources_unchanged=before==after,changed_since_early_dedicated_start=changes,only_expected_post_solve_changes=set(changes)<=set(allowed),source_phase_note='Early dedicated evidence groups began before post-solve P6 registration and the post-solve reporter/schema-default correction and packaging metadata correction. No physics, numerical settings, existing test assertions or diagnostic solve code changed; the new small default-load test was extended before either core group ran. The complete core groups use the final scientific sources; all core sources are bound before and after.')
(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
raise SystemExit(int(not result['core_sources_unchanged'] or not result['only_expected_post_solve_changes'] or any(r['exit_code'] for r in records)))
