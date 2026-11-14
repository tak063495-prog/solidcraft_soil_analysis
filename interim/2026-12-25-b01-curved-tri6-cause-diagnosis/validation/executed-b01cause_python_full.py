import concurrent.futures,hashlib,json,os,pathlib,subprocess,sys,time
root=pathlib.Path.cwd();out=root/'work/b01cause-validation/python-full';out.mkdir(exist_ok=True);runtime=root/'work/b01cause-powershell'
env=dict(os.environ,PATH=str(runtime/'runtime')+os.pathsep+os.environ['PATH'],XDG_CACHE_HOME=str(runtime/'cache'),XDG_CONFIG_HOME=str(runtime/'config'),XDG_DATA_HOME=str(runtime/'data'))
paths=[p for p in subprocess.check_output(['git','ls-files'],text=True).splitlines() if pathlib.Path(p).suffix in ('.rs','.py','.ps1','.toml') or pathlib.Path(p).name in ('Cargo.lock','rust-toolchain.toml')]
before={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in sorted(paths)};(out/'source-binding-before.json').write_text(json.dumps(before,indent=2)+'\n')
groups=['python-mean-evidence','python-core-0','python-core-1','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause']
def run(g):
 cmd=[sys.executable,'tools/run_ci_python_evidence_shards.py','--group',g,'--output',str(out/g)];t=time.monotonic()
 with (out/(g+'-driver.log')).open('x',encoding='utf-8') as log:p=subprocess.Popen(cmd,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT);code=p.wait()
 row=dict(group=g,command=cmd,pid=p.pid,exit_code=code,elapsed_seconds=time.monotonic()-t);(out/(g+'-process.json')).write_text(json.dumps(row,indent=2)+'\n');print(json.dumps(row),flush=True);return row
t=time.monotonic()
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:rows=list(pool.map(run,groups))
after={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in before};(out/'source-binding-after.json').write_text(json.dumps(after,indent=2)+'\n');report=dict(groups=rows,source_count=len(before),sources_unchanged=before==after,changed_sources=[p for p in before if before[p]!=after[p]],all_groups_attempted=True,elapsed_seconds=time.monotonic()-t)
(out/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
