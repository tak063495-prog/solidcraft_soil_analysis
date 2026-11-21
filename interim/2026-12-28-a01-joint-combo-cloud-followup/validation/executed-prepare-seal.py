from pathlib import Path
import json,hashlib,shutil
root=Path.cwd();here=root/'interim/2026-12-28-a01-joint-combo-cloud-followup';out=root/'work/a01follow-seal-metadata';out.mkdir()
# Recover the actual loaded controller versions byte-for-byte against their
# pre-run hashes; later timeout/Windows handling is not attributed to old runs.
newbatch=(here/'a01follow_batch.py').read_text()
newblock='''            p=subprocess.Popen(args,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
            timed_out=False
            try:code=p.wait(timeout=900)
            except subprocess.TimeoutExpired:
                timed_out=True;p.kill();code=p.wait()
                folder=out/name
                if not (folder/'result.json').exists():
                    states=sum(1 for _ in (folder/'fields.jsonl').open()) if (folder/'fields.jsonl').exists() else 0
                    write(folder/'result.json',dict(complete=False,diagnostic_only=True,pid=p.pid,
                        predictor=pred,repartition=part,accepted_states=states,elapsed_seconds=time.monotonic()-started,
                        failure_reason='900-second diagnostic wall-clock budget exhausted; partials retained',
                        numerical_outcome='incomplete',scope_refusal=False,physical_failure=False))
        result=dict(name=name,task=task,pid=p.pid,command=args,exit_code=code,timed_out=timed_out,elapsed_seconds=time.monotonic()-started)'''
oldblock='''            p=subprocess.Popen(args,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT);code=p.wait(timeout=1000)
        result=dict(name=name,task=task,pid=p.pid,command=args,exit_code=code,elapsed_seconds=time.monotonic()-started)'''
assert newblock in newbatch
extra=newbatch.replace(newblock,oldblock).replace("assert all(r['exit_code']==0 or r['timed_out'] for r in results)","assert all(r['exit_code']==0 for r in results)")
extra=extra.replace("files += [build/'target/release'/(name+('.exe' if os.name=='nt' else '')) for name", "files += [build/'target/release'/name for name")
core=extra.replace("    p.add_argument('--stacked-offdiagonal',action='store_true');a=p.parse_args()", "    a=p.parse_args()")
# Original parser had the parse on its first line.
core=core.replace("p.add_argument('--build',type=Path,required=True)\n    a=p.parse_args()", "p.add_argument('--build',type=Path,required=True);a=p.parse_args()")
core=core.replace("    if a.stacked_offdiagonal:tasks=[('stacked',engine,pred,part,8,False) for engine,pred,part in itertools.product(('native','independent'),repeat=3) if pred!=part]\n",'')
v1=core.replace("    files=[HERE/name for name in ('SPEC.ja.md','a01follow_binding.py','a01follow_build.py','a01follow_service.rs',\n        'a01follow_policies.py','a01follow_reference.py','a01follow_run.py','a01follow_batch.py')]", "    files=sorted(p for p in HERE.iterdir() if p.suffix in ('.py','.rs','.md'))")
build=(here/'a01follow_build.py').read_text()
build=build.replace("            exe=name+('.exe' if os.name=='nt' else '')\n            if target!=out/'target':shutil.copy2(target/'release'/exe,out/'target/release'/exe)","            if target!=out/'target':shutil.copy2(target/'release'/name,out/'target/release'/name)")
build=build.replace("(out/'target/release'/(name+('.exe' if os.name=='nt' else '')))","(out/'target/release'/name)")
run=(here/'a01follow_run.py').read_text()
run=run.replace("        timed_out=False\n        try:p.wait(timeout=900)\n        except subprocess.TimeoutExpired:timed_out=True;p.kill();p.wait()", "        p.wait(timeout=900)")
run=run.replace('pid=p.pid,timed_out=timed_out,','pid=p.pid,')
run=run.replace("    else:result.update(complete=False,accepted_states=0,failure_reason='diagnostic process incomplete; partials retained')","    else:result['complete']=False")
versions={}
for label,folder,batch in [('pre-bridge-fix','a01follow-formal',v1),('formal','a01follow-formal-v3',core),('stacked-extra','a01follow-stacked-extra',extra)]:
 binding=json.load(open(root/'work'/folder/'source-binding-before.json'));manifest={}
 for path,digest in binding.items():
  p=root/path
  if '/target/release/' in path:continue
  raw=p.read_bytes()
  if path.endswith('/a01follow_batch.py'):raw=batch.encode()
  elif path.endswith('/a01follow_build.py'):raw=build.encode()
  elif path.endswith('/a01follow_run.py'):raw=run.encode()
  elif label=='pre-bridge-fix' and path.endswith('/a01follow_reference.py'):raw=raw.replace(b'            m._a01follow_generations=generations.copy()\n',b'')
  elif label=='pre-bridge-fix' and path.endswith('/a01follow_policies.py'):raw=raw.replace(b'born=born,generations=m._a01follow_generations,',b'born=born,')
  elif label=='pre-bridge-fix' and path.endswith('/a01follow_service.rs'):raw=(root/'work/a01follow-final-build/source/src/a01follow_service.rs').read_bytes()
  actual=hashlib.sha256(raw).hexdigest()
  if actual!=digest:raise ValueError((label,path,actual,digest))
  target=out/'source-versions'/label/path;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw);manifest[path]=digest
 versions[label]=manifest
(out/'loaded-source-versions.json').write_text(json.dumps(versions,indent=2)+'\n')
for name in ('a01follow-build','a01follow-final-build','a01follow-build-v3'):
 target=out/'builds'/name;target.mkdir(parents=True)
 for p in (root/'work'/name).iterdir():
  if p.is_file():shutil.copy2(p,target/p.name)
 shutil.copytree(root/'work'/name/'source',target/'source')
for name in ('a01follow-kernel-build','a01follow-kernel-build-v2'):
 target=out/'builds'/name;target.mkdir(parents=True)
 for p in (root/'work'/name).iterdir():
  if p.is_file() and p.suffix in ('.json','.log'):shutil.copy2(p,target/p.name)
current=out/'final-diagnostic-source';current.mkdir()
for p in here.iterdir():
 if p.suffix in ('.py','.rs','.md'):shutil.copy2(p,current/p.name)
shutil.copy2(root/'tests/test_a01_joint_combo_cloud_followup.py',current/'test_a01_joint_combo_cloud_followup.py')
for name in ('archive-acquisition.json','archive-acquisition-zip-gzip.json','acquire_archives.py','a01follow-base-preservation.json','a01follow-preparation-events.json','a01follow-runtime.json','a01follow_timeout_watch.py','p601-required-oids.txt','p601-fetch.log','a01follow-batch.log','a01follow-batch-v3.log','a01follow-stacked-extra.log','base-build.log','a01follow-analytic-test.log','a01follow-analytic-test-final.log','a01follow-independent-smoke.log','a01follow-independent-smoke-2.log','a01follow-separation-preparation.log','a01follow-report.log'):
 p=root/'work'/name
 if p.exists():shutil.copy2(p,out/name)
print('source versions verified', {k:len(v) for k,v in versions.items()})
