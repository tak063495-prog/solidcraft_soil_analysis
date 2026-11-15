"""Precision accounting for the same saved rounded matrices, no new physics."""
import argparse
from decimal import Decimal, localcontext
from pathlib import Path
import time
import numpy as np
from b01follow_common import metric, read, rotate, write


def decimal_solve(record,field,digits):
    free,bd=record['free_dofs'],record['boundary_dofs'];row=next(s for s in record['solves'] if s['name']==field)
    with localcontext() as context:
        context.prec=digits;number=lambda v:Decimal.from_float(float(v))
        d=[number(v) for v in row['displacement_rotation']];f=[number(v) for v in row['applied_force_moment']]
        k=record['stiffness'];a=[[number(k[i][j]) for j in free] for i in free]
        b=[f[i]-sum((number(k[i][j])*d[j] for j in bd),Decimal(0)) for i in free]
        for j in range(len(free)):
            pivot=max(range(j,len(free)),key=lambda i:abs(a[i][j]))
            if a[pivot][j]==0:raise ValueError('singular fixed matrix')
            a[j],a[pivot]=a[pivot],a[j];b[j],b[pivot]=b[pivot],b[j]
            for i in range(j+1,len(free)):
                factor=a[i][j]/a[j][j];a[i][j]=Decimal(0)
                if factor:
                    for c in range(j+1,len(free)):a[i][c]-=factor*a[j][c]
                    b[i]-=factor*b[j]
        solved=[Decimal(0)]*len(free)
        for i in range(len(free)-1,-1,-1):
            solved[i]=(b[i]-sum((a[i][j]*solved[j] for j in range(i+1,len(free))),Decimal(0)))/a[i][i]
        for i,v in zip(free,solved):d[i]=v
        return d


def run(formal,measurements,output):
    output.mkdir(parents=True,exist_ok=False);m=read(measurements);jobs=[]
    for label in ('planar','faceted','rigid'):
        candidates=[]
        for row in m['cases']:
            if row['rotation']=='identity':continue
            if label=='planar' and row['shape'].startswith('faceted'):continue
            if label=='faceted' and not row['shape'].startswith('faceted'):continue
            for field in row['covariance']:
                if label=='rigid' and not field.startswith('rigid'):continue
                candidates.append((row['covariance'][field]['displacement_rotation_rad']['maximum_absolute'],row['id'],field))
        maximum,name,field=max(candidates);jobs.append((label,name,field,maximum))
    cache={};rows=[]
    for label,name,field,maximum in jobs:
        shape=name.split('--')[0];base=shape+'--identity';solutions=[]
        for rid in (base,name):
            key=(rid,field);record=read(formal/'assemblies'/(rid+'.json'))
            if key not in cache:
                started=time.monotonic();values=[]
                for digits in (50,80):
                    d=decimal_solve(record,field,digits);values.append(d)
                    write(output/(rid+'--'+field+'--'+str(digits)+'.json'),dict(digits=digits,field=field,
                          exact_decimal_conversion_of_saved_float64=True,displacement_rotation=[str(v) for v in d]))
                a,b=(np.array([float(v) for v in d]) for d in values)
                ordinary=next(s['displacement_rotation'] for s in record['solves'] if s['name']==field)
                cache[key]=dict(d=b,same_float64_roundtrip=bool(np.array_equal(a,b)),
                    digit50_80_difference_max_abs=max(float(abs(a-b)) for a,b in zip(*values)),
                    ordinary_translation_difference_m=metric(b.reshape(-1,6)[:,:3],np.array(ordinary).reshape(-1,6)[:,:3]),
                    ordinary_rotation_difference_rad=metric(b.reshape(-1,6)[:,3:],np.array(ordinary).reshape(-1,6)[:,3:]),
                    elapsed_seconds=time.monotonic()-started)
            q=np.array(record['rotation_matrix']);solutions.append(rotate(cache[key]['d'].reshape(-1,6),q.T))
        rows.append(dict(selection=label,id=name,field=field,ordinary_rotation_covariance_max_abs_rad=maximum,
             decimal_translation_covariance_m=metric(solutions[1][:,:3],solutions[0][:,:3]),
             decimal_rotation_covariance_rad=metric(solutions[1][:,3:],solutions[0][:,3:]),
             limitation='Only RHS formation and solve of fixed rounded K improved. Geometry/basis/kernel rounding remains.'))
    result=dict(passed=False,b01_complete=False,diagnostic_only=True,new_physical_global_solves=0,
                decimal_linear_solves=2*len(cache),cases=rows,
                all_50_80_digit_roundtrips_identical=all(v['same_float64_roundtrip'] for v in cache.values()),
                per_matrix={rid+'--'+f:{k:v for k,v in d.items() if k!='d'} for (rid,f),d in cache.items()})
    write(output/'summary.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--formal',type=Path,required=True)
    p.add_argument('--measurements',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.formal.resolve(),a.measurements.resolve(),a.output.resolve())
