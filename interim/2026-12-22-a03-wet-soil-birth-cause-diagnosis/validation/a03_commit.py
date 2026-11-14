"""Commit a SHA-verified partial-clone index without refetching old archives."""
import os,pathlib,subprocess,sys
msg=sys.argv[1];env=dict(os.environ,GIT_NO_LAZY_FETCH='1')
def git(*args,input=None):return subprocess.check_output(['git',*args],input=input,env=env,text=True).strip()
branch=git('branch','--show-current');assert branch=='codex/plastic-spatial-convergence'
# This repository has no executable commit hooks to bypass.
hooks=subprocess.run(['git','config','--get','core.hooksPath'],capture_output=True,text=True)
assert not hooks.stdout.strip()
active=[p for p in pathlib.Path('.git/hooks').glob('*') if p.is_file() and not p.name.endswith('.sample') and os.access(p,os.X_OK)]
assert not active,active
old=git('rev-parse','HEAD');tree=git('write-tree','--missing-ok')
commit=git('commit-tree',tree,'-p',old,input=msg+'\n')
git('update-ref','-m','commit: '+msg,'refs/heads/'+branch,commit,old)
print(commit)
