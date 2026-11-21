"""PC-only fixed-order full-field audit on 3/4 spatial levels; never accepts C-02."""
import argparse
import importlib.util
import math
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import c02diag_driver as d
from c02diag_transfer import dot


def run(source,out,max_level,steps,user_pc):
    if not user_pc:raise ValueError('explicit --user-pc required; this command is documented, not run in cloud')
    out.mkdir(parents=True,exist_ok=False);c=d.legacy()
    legacy_folder=d.ROOT/'interim/2026-11-19-wet-joint-spatial-convergence'
    sys.path.insert(0,str(legacy_folder))
    spec=importlib.util.spec_from_file_location('c02diag_pc_legacy_mapping',legacy_folder/'mapping.py')
    mapping=importlib.util.module_from_spec(spec);spec.loader.exec_module(mapping)
    maps=mapping.Maps([d.read(source/f'mesh-h{h}.json') for h in range(max_level+1)],d.read(source/'stages.json'))
    result=dict(passed=False,original_acceptance_not_replaced=True,all_fields_fixed_order=True,max_level=max_level,steps=steps,spatial={})
    for start in range(max_level-1):
        routes=[(h,steps) for h in range(start,start+3)];statistics={};count=0
        streams=[c.old.states(source/f'h{h}-g{steps}'/'native') for h in range(start,start+3)]
        for states in zip(*streams,strict=True):
            coords=[(r['stage'],r['progress']) for r,f in states]
            if len(set(coords))!=1:raise ValueError('PC accepted times differ')
            values=[maps.apply(h,f) for h,(r,f) in zip(range(start,start+3),states,strict=True)]
            for key in sorted(values[0]):
                a,b,z=(v[key] for v in values);x=a-b;y=b-z
                row=statistics.setdefault(key,dict(d01=[],d12=[],ref1=[],ref2=[],inner=[],maximum01=0.,maximum12=0.))
                for name,p,q in [('d01',x,x),('d12',y,y),('ref1',b,b),('ref2',z,z),('inner',x,y)]:row[name].append(dot(np.asarray(p).ravel(order='C'),np.asarray(q).ravel(order='C')))
                row['maximum01']=max(row['maximum01'],float(np.max(abs(x))));row['maximum12']=max(row['maximum12'],float(np.max(abs(y))))
            count+=1
        if count!=7*steps:raise ValueError('PC full history incomplete')
        fields={}
        for key,row in statistics.items():
            a,b,r1,r2,t=(math.fsum(row[n]) for n in ('d01','d12','ref1','ref2','inner'))
            fields[key]=dict(absolute01=math.sqrt(a),absolute12=math.sqrt(b),relative01=math.sqrt(a)/max(math.sqrt(r1),1e-12),relative12=math.sqrt(b)/max(math.sqrt(r2),1e-12),order=.5*math.log2(a/b) if min(a,b)>0 else None,difference_cosine=t/math.sqrt(a*b) if min(a,b)>0 else None,maximum01=row['maximum01'],maximum12=row['maximum12'])
        result['spatial'][f'h{start}{start+1}{start+2}']=dict(common_states=count,fields=fields)
    d.write(out/'pc-full-field-audit.json',result)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--max-level',type=int,choices=(2,3),default=3);p.add_argument('--steps',type=int,choices=(32,64,128,256,512),required=True);p.add_argument('--user-pc',action='store_true');a=p.parse_args();run(a.source.resolve(),a.output.resolve(),a.max_level,a.steps,a.user_pc)
