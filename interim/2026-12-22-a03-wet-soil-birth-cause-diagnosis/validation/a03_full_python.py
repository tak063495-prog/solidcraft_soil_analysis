import concurrent.futures,hashlib,json,pathlib,subprocess,sys,time
root=pathlib.Path.cwd();out=root/'work/a03-validation/python-full';out.mkdir(parents=True,exist_ok=False)
paths=[p for p in subprocess.check_output(['git','ls-files'],text=True).splitlines() if pathlib.Path(p).suffix in ('.rs','.py','.ps1','.toml') or pathlib.Path(p).name in ('Cargo.lock','rust-toolchain.toml')]
before={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in sorted(paths)};(out/'source-binding-before.json').write_text(json.dumps(before,indent=2)+'\n')
groups=['python-mean-evidence','python-core-0','python-core-1','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause'];start=time.monotonic()
def run(g):
 command=[sys.executable,'tools/run_ci_python_evidence_shards.py','--group',g,'--output',str(out/g)]
 t=time.monotonic()
 with (out/(g+'-driver.log')).open('x',encoding='utf-8') as f:
  p=subprocess.Popen(command,stdout=f,stderr=subprocess.STDOUT);code=p.wait()
 record=dict(group=g,command=command,pid=p.pid,exit_code=code,elapsed_seconds=time.monotonic()-t)
 (out/(g+'-process.json')).write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record),flush=True);return record
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:records=list(pool.map(run,groups))
after={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in before};(out/'source-binding-after.json').write_text(json.dumps(after,indent=2)+'\n')
result=dict(groups=records,elapsed_seconds=time.monotonic()-start,all_groups_attempted=True,source_count=len(before),sources_unchanged=before==after)
(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
raise SystemExit(int(not result['sources_unchanged'] or any(r['exit_code'] for r in records)))
