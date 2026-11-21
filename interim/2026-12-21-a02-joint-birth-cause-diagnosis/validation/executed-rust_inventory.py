import pathlib,json,time,subprocess,sys,re
root=pathlib.Path('/workspace/solidcraft_soil_analysis');out=root/'work/a02-validation'
while True:
 try:
  rows=json.loads((out/'rust-commands.json').read_text())
  if len(rows)==5:break
 except (OSError,ValueError):pass
 time.sleep(1)
result={}
for label,cmd in [('rust-discovery',['cargo','test','--locked','--release','--','--list']),('rust-examples-discovery',['cargo','test','--locked','--release','--examples','--','--list'])]:
 with (out/(label+'.log')).open('x') as f:r=subprocess.run(cmd,cwd=root,stdout=f,stderr=subprocess.STDOUT)
 lines=(out/(label+'.log')).read_text().splitlines();ids=[line for line in lines if line.endswith(': test')]
 result[label]=dict(command=cmd,exit_code=r.returncode,discovered_count=len(ids),test_names_in_target_order=ids)
(out/'rust-discovery.json').write_text(json.dumps(result,indent=2)+'\n')
print({k:v['discovered_count'] for k,v in result.items()})
