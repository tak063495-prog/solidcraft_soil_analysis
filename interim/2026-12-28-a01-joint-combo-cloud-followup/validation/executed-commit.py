import argparse,os,subprocess
p=argparse.ArgumentParser();p.add_argument('message');a=p.parse_args();env=dict(os.environ,GIT_NO_LAZY_FETCH='1')
def git(*args,**kw):return subprocess.check_output(['git',*args],env=env,text=True,**kw).strip()
base=git('rev-parse','HEAD');tree=git('write-tree','--missing-ok');paths=git('diff-tree','--no-commit-id','--name-only','-r',base,tree).splitlines();allowed=('interim/2026-12-28-a01-joint-combo-cloud-followup/','tests/test_a01_joint_combo_cloud_followup.py','interim/2026-12-10-ci-p601-inventory-fix/inventory.toml','interim/2026-12-06-p6-01-status-reconciliation/audit_status.py')
assert paths and all(any(p.startswith(prefix) for prefix in allowed) for p in paths),paths
sha=git('commit-tree',tree,'-p',base,input=a.message+'\n');subprocess.run(['git','update-ref','refs/heads/codex/plastic-spatial-convergence',sha,base],env=env,check=True);print(sha)
