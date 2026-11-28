import json,pathlib,subprocess,time
root=pathlib.Path.cwd();out=root/'work/b01cause-validation'
for name,cmd in [('rust-discovery',['cargo','test','--locked','--release','--all-targets','--no-run','--message-format=json']),('rust-main',['cargo','test','--locked','--release']),('rust-examples',['cargo','test','--locked','--release','--examples'])]:
 t=time.monotonic()
 with (out/(name+'.log')).open('x') as log:p=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT);code=p.wait()
 (out/(name+'-process.json')).write_text(json.dumps(dict(command=cmd,pid=p.pid,exit_code=code,elapsed_seconds=time.monotonic()-t),indent=2)+'\n');print(name,code,flush=True)
 if name=='rust-discovery' and code==0:
  artifacts=[];ids=[]
  for l in (out/(name+'.log')).read_text().splitlines():
   try:a=json.loads(l)
   except ValueError:continue
   if a.get('reason')=='compiler-artifact' and a.get('executable') and a['profile']['test']:
    kind=a['target']['kind'][0];target=kind+':'+a['target']['name'];exe=a['executable'];artifacts.append(dict(target=target,executable=exe))
    listing=subprocess.check_output([exe,'--list','--format','terse'],text=True)
    ids += [target+'::'+v[:-6] for v in listing.splitlines() if v.endswith(': test')]
  (out/'rust-discovery').mkdir(exist_ok=True)
  (out/'rust-discovery/inventory.json').write_text(json.dumps(dict(total=len(ids),test_ids=sorted(ids),artifacts=artifacts),indent=2)+'\n')
