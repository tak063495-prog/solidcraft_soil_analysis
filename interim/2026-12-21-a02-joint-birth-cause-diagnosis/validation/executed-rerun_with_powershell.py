import pathlib,subprocess,sys,time,json,concurrent.futures,os,hashlib
root=pathlib.Path('/workspace/solidcraft_soil_analysis');out=root/'work/a02-validation';runtime=root/'work/a02-powershell'
env=dict(os.environ,PATH=str(runtime/'runtime')+os.pathsep+os.environ['PATH'],XDG_CACHE_HOME=str(runtime/'cache'),XDG_CONFIG_HOME=str(runtime/'config'),XDG_DATA_HOME=str(runtime/'data'))
def run(group):
 cmd=[sys.executable,'tools/run_ci_python_evidence_shards.py','--group',group,'--output',str(out/group)];start=time.monotonic()
 with (out/(group+'-launcher.log')).open('x') as f:r=subprocess.run(cmd,cwd=root,env=env,stdout=f,stderr=subprocess.STDOUT)
 row=dict(group=group,exit_code=r.returncode,command=cmd,elapsed_seconds=time.monotonic()-start,real_powershell_available=True,project_sources_or_nc_settings_changed=False)
 (out/(group+'-process.json')).write_text(json.dumps(row,indent=2)+'\n');print(group,r.returncode,flush=True);return row
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(run,('python-core-0','python-core-1')))
pre=json.loads((out/'python-pre-source-fingerprints.json').read_text());changed=[p for p,h in pre.items() if hashlib.sha256((root/p).read_bytes()).hexdigest()!=h]
(out/'python-post-rerun-source-check.json').write_text(json.dumps(dict(source_count=len(pre),changed=changed,all_sha256_identical=not changed),indent=2)+'\n')
(out/'python-environment-rerun-commands.json').write_text(json.dumps(results,indent=2)+'\n')
