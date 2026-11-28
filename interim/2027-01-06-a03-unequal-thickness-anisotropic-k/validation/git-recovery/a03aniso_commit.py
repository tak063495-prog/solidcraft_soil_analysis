"""Commit the staged tree without forcing missing historical archive blobs."""
import os,subprocess,sys,json
from pathlib import Path
message=sys.argv[1];env=dict(os.environ,GIT_NO_LAZY_FETCH='1')
def git(*args,input=None):return subprocess.check_output(['git',*args],input=input,env=env,text=True).strip()
head=git('rev-parse','HEAD');branch=git('symbolic-ref','HEAD');assert branch=='refs/heads/codex/plastic-spatial-convergence'
tree=git('write-tree','--missing-ok');commit=git('commit-tree',tree,'-p',head,'-m',message)
git('update-ref',branch,commit,head)
print(json.dumps(dict(commit=commit,parent=head,tree=tree,lazy_fetch_disabled=True,compare_and_swap=True)))
