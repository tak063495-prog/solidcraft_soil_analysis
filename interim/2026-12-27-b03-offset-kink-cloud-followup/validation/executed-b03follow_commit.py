import argparse
import os
import subprocess

p=argparse.ArgumentParser();p.add_argument('message');args=p.parse_args()
env=dict(os.environ,GIT_NO_LAZY_FETCH='1')
def git(*args,**kwargs):return subprocess.check_output(['git',*args],env=env,text=True,**kwargs).strip()
base=git('rev-parse','HEAD');tree=git('write-tree','--missing-ok')
paths=git('diff-tree','--no-commit-id','--name-only','-r',base,tree).splitlines()
allowed=('interim/2026-12-27-b03-offset-kink-cloud-followup/',
         'tests/test_b03_offset_kink_cloud_followup.py',
         'interim/2026-12-10-ci-p601-inventory-fix/inventory.toml',
         'interim/2026-12-06-p6-01-status-reconciliation/audit_status.py')
if not paths or not all(any(name.startswith(prefix) for prefix in allowed) for name in paths):
    raise ValueError(paths)
sha=git('commit-tree',tree,'-p',base,input=args.message+'\n')
subprocess.run(['git','update-ref','refs/heads/codex/plastic-spatial-convergence',sha,base],check=True,env=env)
print(sha)
