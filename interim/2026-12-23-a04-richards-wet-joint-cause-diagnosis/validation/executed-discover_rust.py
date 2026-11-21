import json,pathlib,subprocess
root=pathlib.Path('/workspace/solidcraft_soil_analysis');out=root/'work/a04cause-validation/rust-discovery';out.mkdir(exist_ok=False)
artifacts=[];ids=[]
for line in (root/'work/a04cause-validation/rust-discovery.log').read_text().splitlines():
 try:m=json.loads(line)
 except ValueError:continue
 if m.get('reason')!='compiler-artifact' or not m.get('executable') or not m['profile']['test']:continue
 target=m['target']['kind'][0]+':'+m['target']['name'];cmd=[m['executable'],'--list','--format','terse'];text=subprocess.check_output(cmd,cwd=root,text=True)
 names=[line.removesuffix(': test') for line in text.splitlines() if line.endswith(': test')]
 artifacts.append(dict(target=target,executable=m['executable'],command=cmd,test_names=names));ids.extend(target+'::'+name for name in names)
assert len(ids)==len(set(ids))
result=dict(total=len(ids),test_ids=sorted(ids),artifacts=artifacts)
(out/'inventory.json').write_text(json.dumps(result,indent=2)+'\n');print('discovered',len(ids),'targets',len(artifacts))
