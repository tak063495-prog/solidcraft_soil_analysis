"""Single-parent commit from verified index, retaining missing promisor old references."""
import argparse,json,os,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--message',required=True);p.add_argument('--path',action='append',required=True);p.add_argument('--proof',type=Path,required=True);a=p.parse_args()
env=dict(os.environ,GIT_NO_LAZY_FETCH='1')
def git(*args):return subprocess.check_output(['git',*args],env=env,text=True).strip()
old=git('rev-parse','HEAD');ref=git('symbolic-ref','HEAD')
if ref!='refs/heads/codex/plastic-spatial-convergence':raise ValueError('wrong target branch')
subprocess.run(['git','add','--sparse','--',*a.path],env=env,check=True)
paths=git('diff','--cached','--name-only').splitlines()
allowed=('interim/2027-01-02-b02-offset-composite-cloud-followup/','tests/test_b02_offset_composite_cloud_followup.py','interim/2026-12-10-ci-p601-inventory-fix/inventory.toml','interim/2026-12-06-p6-01-status-reconciliation/audit_status.py')
if any(not any(x.startswith(v) if v.endswith('/') else x==v for v in allowed) for x in paths):raise ValueError('unexpected staged path')
tree=git('write-tree','--missing-ok')
sha=subprocess.check_output(['git','commit-tree',tree,'-p',old],input=a.message+'\n',text=True,env=env).strip()
subprocess.run(['git','update-ref',ref,sha,old],check=True,env=env)
subprocess.run(['git','merge-base','--is-ancestor','23c7b0ea48906723620b47eefb50e7481074d875',sha],check=True,env=env)
a.proof.parent.mkdir(parents=True,exist_ok=True)
a.proof.write_text(json.dumps(dict(parent=old,commit=sha,tree=tree,paths=paths,single_parent=True,verified_base_ancestor=True,old_cache_references_retained=True),indent=2)+'\n')
print(sha)
