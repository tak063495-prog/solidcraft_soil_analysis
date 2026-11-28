"""Independent basis-change check for oriented section resultants, not a solve."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import argparse
from pathlib import Path
import numpy as np
from b03_cases import write
from plate6_offset_reference import Element
from convergence_b03 import tensor


def run(output):
    points=np.array([[0.,0.,1.],[1.,0.,1.],[.3,1.,1.],[.5,0.,1.],[.65,.5,1.],[.15,.5,1.]])
    sec=dict(young_kpa=3e7,poisson=.2,thickness_m=.1,drilling_factor=1e-4)
    permutation=[0,2,1,5,4,3];l=np.array([.2,.3,.5]);rows=[]
    for e in (0.,.05,-.05):
        a=Element(points,sec,e);b=Element(points[permutation],sec,-e)
        d=np.array([[1e-5*(x*x+2*y),2e-5*(y*y-x),3e-5*x*y,4e-5*x,5e-5*y,6e-5*(x+y)] for x,y,z in points])
        def fields(k,kin,ll):
            local=k.neutral.transform@k.map@kin.ravel();bm,bb,bq,_=k.neutral.operators(ll)
            n=k.neutral.membrane@(bm@local);m=k.neutral.bending@(bb@local);q=k.neutral.shear*(bq@local)
            return dict(membrane=tensor(n,k.axes),bending=tensor(m,k.axes),attachment_moment=tensor(m-k.e*n,k.axes),shear=k.axes.T@np.r_[q,0.])
        f=fields(a,d,l);r=fields(b,d[permutation],l[[0,2,1]])
        converted={k:v if k=='membrane' else -v for k,v in r.items()}
        rows.append(dict(e_m=e,normal_first=a.normal.tolist(),normal_reversed=b.normal.tolist(),
            fields_first={k:v.tolist() for k,v in f.items()},fields_reversed_raw={k:v.tolist() for k,v in r.items()},
            fields_reversed_oriented_to_first={k:v.tolist() for k,v in converted.items()},
            maximum_difference_by_field={k:float(np.max(np.abs(f[k]-converted[k]))) for k in f}))
    write(output,dict(passed=False,b03_complete=False,accepted=None,rows=rows,
        scope='Coordinate/sign transformation diagnostic using independent polynomial operators and synthetic physical kinematics. Does not generate a reference solution or normalize an input. No new engineering acceptance gate.'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args();run(a.output)
