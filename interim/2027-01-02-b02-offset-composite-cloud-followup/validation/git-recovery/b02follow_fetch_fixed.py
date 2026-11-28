import json,tomllib,subprocess,os
from pathlib import Path
out=Path('work/b02follow-preparation');env=dict(os.environ,GIT_NO_LAZY_FETCH='1')
ledger=json.loads(Path('interim/2026-12-06-p6-01-status-reconciliation/ledger.json').read_text())
sup=tomllib.loads(Path('interim/2026-12-10-ci-p601-inventory-fix/inventory.toml').read_text())
commits=sorted({v['commit'] for v in list(ledger['sources'].values())+sup['sources'] if v.get('origin')!='github_snapshot'})
proof=[]
for n,ids in enumerate([commits[i:i+16] for i in range(0,len(commits),16)]):
 missing=[i for i in ids if subprocess.run(['git','cat-file','-e',i+'^{commit}'],env=env,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL).returncode]
 if not missing:continue
 cmd=['git','-c','fetch.negotiationAlgorithm=noop','-c','gc.auto=0','fetch','--depth=1','--filter=blob:none','--no-tags','--no-write-fetch-head','--stdin','origin']
 p=subprocess.run(cmd,input='\n'.join(missing)+'\n',text=True,capture_output=True,env=env)
 (out/f'fixed-commits-{n}.log').write_text(p.stdout+p.stderr)
 proof.append(dict(command=cmd,ids=missing,exit_code=p.returncode));assert p.returncode==0
m=json.loads(Path('interim/2026-12-07-p6-03-cloud-reconciliation/branch-sources.json').read_text())
ids=sorted({v['git_blob_sha1'] for v in m['files']})
cmd=['git','-c','fetch.negotiationAlgorithm=noop','-c','gc.auto=0','fetch','--filter=blob:none','--no-tags','--no-write-fetch-head','--stdin','origin']
p=subprocess.run(cmd,input='\n'.join(ids)+'\n',text=True,capture_output=True,env=env)
(out/'fixed-p603-blobs.log').write_text(p.stdout+p.stderr);proof.append(dict(command=cmd,ids=ids,exit_code=p.returncode));assert p.returncode==0
(out/'fixed-source-recovery.json').write_text(json.dumps(proof,indent=2)+'\n')
print('fixed commits',len(commits),'p603 blobs',len(ids),flush=True)
