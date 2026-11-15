import os,subprocess,json,hashlib,time
from pathlib import Path
out=Path('interim/2027-01-02-b02-offset-composite-cloud-followup/validation/git-recovery');out.mkdir(exist_ok=True)
m=json.loads(Path('interim/2026-12-07-p6-03-cloud-reconciliation/branch-sources.json').read_text());base=m['base_commit'];env=dict(os.environ,GIT_NO_LAZY_FETCH='1')
cmd=['git','-c','fetch.negotiationAlgorithm=noop','-c','gc.auto=0','fetch','--depth=1','--filter=blob:none','--no-tags','--no-write-fetch-head','origin',base]
start=time.monotonic();p=subprocess.run(cmd,text=True,capture_output=True,env=env)
(out/'p603-base-fetch.log').write_text(p.stdout+p.stderr);assert p.returncode==0
proof=[]
for item in m['files']:
 raw=subprocess.check_output(['git','show',base+':'+item['path']],env=env)
 digest=hashlib.sha256(raw).hexdigest();blob=hashlib.sha1(('blob '+str(len(raw))+'\0').encode()+raw).hexdigest()
 assert digest==item['sha256'] and blob==item['git_blob_sha1']
 proof.append(dict(path=item['path'],sha256=digest,git_blob_sha1=blob,exact=True))
r=dict(command=cmd,base=base,exit_code=p.returncode,elapsed_seconds=time.monotonic()-start,proofs=proof,
 failure='Initial full core-0: P6-03 BASE commit/tree absent in depth=1 checkout; all 29 frozen blobs had already been fetched.',
 restoration='Only fixed BASE commit/tree fetched with blob:none/depth=1; all 29 original SHA and Git blobs matched. No source/reference/tolerance change.')
(out/'p603-base-recovery.json').write_text(json.dumps(r,indent=2)+'\n');print('P603 fixed base recovered and',len(proof),'files verified')
