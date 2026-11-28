"""Commit the exact staged tree in a verified partial clone, preserving parent."""
import argparse,os,subprocess
p=argparse.ArgumentParser();p.add_argument('message');a=p.parse_args()
env=dict(os.environ,GIT_NO_LAZY_FETCH='1')
def git(*args,**kw):return subprocess.check_output(['git',*args],env=env,text=True,**kw).strip()
parent=git('rev-parse','HEAD');tree=git('write-tree','--missing-ok')
paths=git('diff-tree','--no-commit-id','--name-only','-r',parent,tree).splitlines()
allowed=('interim/2026-12-24-b02-offset-composite-cause-diagnosis/','tests/test_b02_offset_composite_cause_diagnosis.py','interim/2026-12-10-ci-p601-inventory-fix/inventory.toml','interim/2026-12-06-p6-01-status-reconciliation/audit_status.py')
assert paths and all(any(path.startswith(prefix) for prefix in allowed) for path in paths),paths
sha=git('commit-tree',tree,'-p',parent,input=a.message+'\n')
subprocess.run(['git','update-ref','refs/heads/codex/plastic-spatial-convergence',sha,parent],check=True,env=env)
print(sha)
