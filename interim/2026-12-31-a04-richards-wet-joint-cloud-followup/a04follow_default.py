"""Compare isolated default controls with the ordinary production CLI."""
import argparse
import shutil
from a04follow_common import *

def strip(v):
    if isinstance(v,dict):return {k:strip(x) for k,x in v.items() if k not in ('elapsed_seconds','solver_binary_sha256')}
    if isinstance(v,list):return [strip(x) for x in v]
    return v

def main():
    p=argparse.ArgumentParser();p.add_argument('--soil',type=Path,required=True);p.add_argument('--order',type=Path,required=True);p.add_argument('--exchange',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False);solver=ROOT/('target/release/soil-fem-3d'+('.exe' if os.name=='nt' else ''));result={}
    inputs=[('soil',a.soil/'soil/h2-n4-q3','solve-one-way-construction','native.json'),('dry',a.order/'order/dry/e0f0q0','solve-stages','native/run.json'),('hydro',a.order/'order/hydro/e0f0q0','solve-stages','native/run.json'),('c04',a.exchange/'exchange/native-pNrN','solve-stages','raw/run.json')]
    for name,source,command,old in inputs:
        folder=a.output/name;folder.mkdir()
        for k in ('mesh','ground','case'):
            path=(a.exchange/'exchange'/(k+'.json')) if name=='c04' and k!='case' else source/(k+'.json')
            shutil.copy2(path,folder/(k+'.json'))
        proc=invoke([solver,command,*[folder/(k+'.json') for k in ('mesh','ground','case')],folder/'run.json'],folder,'normal-cli')
        x=read(folder/'run.json');y=read(source/old);same=strip(x)==strip(y)
        write(folder/'comparison.json',dict(whole_run_exact_except_elapsed_and_producer_sha=same,normal_cli_sha256=sha(solver),process=proc))
        if not same:raise ValueError(('isolated default mismatch',name))
        result[name]=dict(exact=True,states=len(x['frames'])+len(x['initialization']['frames'])-1 if name=='soil' else len(x['attempts'])+len(x['geostatic']['attempts']))
    write(a.output/'result.json',dict(execution_complete=True,normal_cli_sha256=sha(solver),results=result,production_changed=False))

if __name__=='__main__':main()
