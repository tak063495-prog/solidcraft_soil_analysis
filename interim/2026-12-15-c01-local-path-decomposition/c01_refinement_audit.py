"""Rebuild finite-reference budgets from the sealed point trace, no integration."""
import argparse
import json
from math import fsum
from decimal import Decimal,localcontext
from pathlib import Path
import numpy as np
from c01_decomposition import COMPONENTS,FIELDS,STAGES,dot,metric,read
from c01_local_study import write


def norm_ratio(a,b):
    numerator,denominator=dot(a,a),dot(b,b)
    if denominator==0:return None
    with localcontext() as context:
        context.prec=80
        return float(Decimal.from_float(numerator).sqrt()/Decimal.from_float(denominator).sqrt())


def build(report,trace):
    data={}
    with trace.open(encoding='utf-8') as f:
        for line in f:
            r=json.loads(line)
            key=(r['control'],r['kind'],r['steps'])
            coord=(r['stage'],r['progress'])
            series=data.setdefault(key,{})
            if coord in series:raise ValueError('duplicate trace coordinate')
            series[coord]=r['fields']
    output=dict(schema='c01-finite-reference-budget-audit/1',new_global_equilibrium_solved=False,c01_passed=False,pairs={})
    for label,pair in report['pairs'].items():
        control,kind,levels=label.split('/');a,b=map(int,levels.split('__'))
        da,db=data[(control,kind,a)],data[(control,kind,b)]
        coordinates=[(s,i/a) for s in STAGES for i in range(1,a+1)]
        if coordinates!=[(r['stage'],r['progress']) for r in pair['intervals']]:raise ValueError('pair interval order differs')
        rows={}
        for field in FIELDS:
            original_a=np.asarray([da[c][field]['stored'] for c in coordinates]);original_b=np.asarray([db[c][field]['stored'] for c in coordinates]);delta=original_b-original_a
            m=metric(original_a,original_b)
            if m!=pair['fields'][field]['original']:raise ValueError('original metric changed in trace')
            versions={}
            for tag in ('coarser_replay','refined_replay'):
                ra=np.asarray([da[c][field][tag] for c in coordinates]);rb=np.asarray([db[c][field][tag] for c in coordinates]);path=rb-ra
                local=delta-path;switch=np.zeros_like(path[0]);previous=np.zeros_like(switch);sw=[]
                for i,r in enumerate(pair['intervals']):
                    mask=np.zeros(len(switch),dtype=bool);mask[r['switch_timing_mismatch_gp']]=True
                    if switch.ndim==2:mask=mask[:,None]
                    switch+=np.where(mask,path[i]-previous,0.)
                    previous=path[i];sw.append(switch.copy())
                switched=np.stack(sw)
                cs=dict(zip(COMPONENTS,(local,path-switched,switched)))
                total_squared=dot(delta,delta)
                versions[tag]=dict(components=cs,projection={c:dot(v,delta)/total_squared if total_squared else None for c,v in cs.items()})
            change=versions['refined_replay']['components']['local_integration']-versions['coarser_replay']['components']['local_integration']
            local=versions['refined_replay']['components']['local_integration']
            previous_local=np.zeros_like(local[0]);previous_total=np.zeros_like(delta[0]);increments=[]
            for i,(stage,progress) in enumerate(coordinates):
                di=local[i]-previous_local;dt=delta[i]-previous_total
                increments.append(dict(stage=stage,progress=progress,local_defect_increment_squared=dot(di,di),stored_difference_increment_squared=dot(dt,dt),local_defect_increment_signed_sum=fsum(float(x) for x in di.flat)))
                previous_local=local[i];previous_total=delta[i]
            rows[field]=dict(original=m,largest_local_increment_defects=sorted(increments,key=lambda r:r['local_defect_increment_squared'],reverse=True)[:12],coarser_projection=versions['coarser_replay']['projection'],refined_projection=versions['refined_replay']['projection'],reference_change_squared=dot(change,change),reference_change_norm_over_original_difference=norm_ratio(change,delta))
        output['pairs'][label]=rows
    return output


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--report',type=Path,required=True);p.add_argument('--points',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();write(a.output,build(read(a.report),a.points))

if __name__=='__main__':main()
