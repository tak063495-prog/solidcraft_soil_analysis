import pathlib,subprocess,sys,json,time,hashlib,concurrent.futures
root=pathlib.Path('/workspace/solidcraft_soil_analysis');out=root/'work/a02-validation';base='7b254b1cdbe9c288989631944cd698aa861f83ac'
codepaths=sorted(p for name in ('src','tests','tools','examples','interim') for p in (root/name).rglob('*') if p.is_file() and p.suffix in ('.py','.rs') and '__pycache__' not in p.parts)
fingerprints={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in codepaths}
(out/'python-pre-source-fingerprints.json').write_text(json.dumps(fingerprints,indent=2)+'\n')
meanpaths=subprocess.check_output(['git','ls-tree','-r','--name-only',base,'--','interim/2026-11-10-nc-mean-stress-release'],cwd=root,text=True).splitlines()
started=time.monotonic()
def check_ready(paths):return all((root/p).is_file() for p in paths)
def run(group):
 command=[sys.executable,'tools/run_ci_python_evidence_shards.py','--group',group,'--output',str(out/group)]
 start=time.monotonic()
 with (out/(group+'-launcher.log')).open('x') as f:r=subprocess.run(command,cwd=root,stdout=f,stderr=subprocess.STDOUT)
 row=dict(group=group,exit_code=r.returncode,elapsed_seconds=time.monotonic()-start,command=command)
 (out/(group+'-process.json')).write_text(json.dumps(row,indent=2)+'\n');print(group,r.returncode,flush=True);return row
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 while not check_ready(meanpaths):time.sleep(1)
 print('All mean-stress source archives present; starting that dedicated group',flush=True)
 jobs=[pool.submit(run,'python-mean-evidence')]
 while not (root/'work/a02-materialization.json').exists():time.sleep(1)
 print('All original archive bytes restored; starting remaining groups',flush=True)
 for group in ('python-core-0','python-core-1','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause'):jobs.append(pool.submit(run,group))
 results=[job.result() for job in jobs]
changed=[p for p,d in fingerprints.items() if hashlib.sha256((root/p).read_bytes()).hexdigest()!=d]
(out/'python-post-source-check.json').write_text(json.dumps(dict(source_count=len(fingerprints),changed=changed,all_sha256_identical=not changed),indent=2)+'\n')
(out/'python-all-commands.json').write_text(json.dumps(results,indent=2)+'\n')
print('done',time.monotonic()-started,flush=True)
