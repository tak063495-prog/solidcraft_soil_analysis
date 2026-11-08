"""Re-run only the existing small case. No large study or production acceptance."""
import argparse,json,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]

def commands(w):
    py=sys.executable;cmds=[]
    def add(label,script,*args):cmds.append((label,[py,str(script),*map(str,args)]))
    for name,n,method in [('independent-input-midpoint',200,'midpoint'),('independent-input-rk4',200,'rk4'),('independent-input-rk4-fine',400,'rk4'),('independent-input-rk4-800',800,'rk4')]:
        add(name,ROOT/'tools/low_preload_nc_reference.py','--output',w/name,'--steps',n,'--method',method,'--budget-seconds',180)
    add('independent-input-rk4-800-history',HERE/'capture_reference_history.py','--output',w/'independent-input-rk4-800-history','--steps',800)
    binary=w/'target/release/examples'/('material_diagnostic.exe' if sys.platform=='win32' else 'material_diagnostic')
    for name,ds,n,resume in [('native-coarse',2e-5,128,None),('native-fine',1e-5,256,None),('native-restart-a',2e-5,64,None),('native-restart-b',2e-5,64,'native-restart-a')]:
        args=['--evidence',w/'prior','--binary',binary,'--output',w/name,'--step-size',ds,'--count',n]
        if resume:args+=['--resume',w/resume/'restart.json']
        add(name,HERE/'native_arc.py',*args)
    add('native-branch-search',HERE/'branch_search.py','--evidence',w/'prior','--binary',binary,'--output',w/'native-branch-search')
    for name,ref,ds,n,resume in [('independent-arc-coarse','independent-input-rk4',2e-5,128,None),('independent-arc-half','independent-input-rk4',1e-5,256,None),('independent-arc-extended200','independent-input-rk4',1e-4,256,None),('independent-arc-history400','independent-input-rk4-fine',2e-5,128,None),('independent-arc-history400-half','independent-input-rk4-fine',1e-5,256,None),('independent-arc-extended400','independent-input-rk4-fine',1e-4,256,None),('independent-arc-history800','independent-input-rk4-800',2e-5,128,None),('independent-restart-a','independent-input-rk4-fine',2e-5,64,None),('independent-restart-b','independent-input-rk4-fine',2e-5,64,'independent-restart-a')]:
        args=['--reference',w/ref,'--output',w/name,'--step-size',ds,'--count',n]
        if resume:args+=['--resume',w/resume/'restart.json']
        add(name,HERE/'independent_arc.py',*args)
    for n in [4,8]:add('independent-branch-search'+str(n),HERE/'independent_branches.py','--reference',w/'independent-input-rk4-fine','--arc',w/'independent-arc-history400','--output',w/('independent-branch-search'+str(n)),'--material-steps',n)
    for name,ds,n,steps in [('independent-map-arc-coarse',2e-5,64,4),('independent-map-arc-fine',1e-5,128,4),('independent-map-arc-material8',2e-5,64,8)]:
        add(name,HERE/'independent_map_arc.py','--reference',w/'independent-input-rk4-fine','--arc',w/'independent-arc-history400','--output',w/name,'--step-size',ds,'--count',n,'--material-steps',steps)
    return cmds

def run(w):
    # Each solver writes diagnostics with passed=false even if its process succeeds.
    records=[]
    def invoke(name,args):
        started=time.monotonic()
        with (w.parent/(w.name+'-'+name+'.log')).open('x') as log:r=subprocess.run(args,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,timeout=900)
        records.append(dict(label=name,command=args,exit_code=r.returncode,elapsed_seconds=time.monotonic()-started));r.check_returncode()
    invoke('prepare',[sys.executable,str(HERE/'prepare_evidence.py'),'--output',str(w)])
    for name,args in commands(w):invoke(name,args)
    invoke('analyze',[sys.executable,str(HERE/'analyze.py'),'--evidence',str(w),'--output',str(w/'analysis')])
    (w/'reproduction-processes.json').write_text(json.dumps(records,indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--list-only',action='store_true');a=p.parse_args();w=a.output.resolve()
    if a.list_only:print(json.dumps(commands(w),indent=2))
    else:run(w)
