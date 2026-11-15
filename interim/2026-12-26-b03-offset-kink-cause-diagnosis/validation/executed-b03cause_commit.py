import argparse,os,subprocess
p=argparse.ArgumentParser();p.add_argument('message');a=p.parse_args();env=dict(os.environ,GIT_NO_LAZY_FETCH='1')
def git(*args,**kw):return subprocess.check_output(['git',*args],env=env,text=True,**kw).strip()
parent=git('rev-parse','HEAD');tree=git('write-tree','--missing-ok');paths=git('diff-tree','--no-commit-id','--name-only','-r',parent,tree).splitlines();allowed=('interim/2026-12-26-b03-offset-kink-cause-diagnosis/','tests/test_b03_offset_kink_cause_diagnosis.py','interim/2026-12-10-ci-p601-inventory-fix/inventory.toml','interim/2026-12-06-p6-01-status-reconciliation/audit_status.py');assert paths and all(any(x.startswith(y) for y in allowed) for x in paths),paths
sha=git('commit-tree',tree,'-p',parent,input=a.message+'\n');subprocess.run(['git','update-ref','refs/heads/codex/plastic-spatial-convergence',sha,parent],check=True,env=env);print(sha)
