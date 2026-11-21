"""Finite rotated boundary inputs; order audit never changes the native gate."""
import argparse,itertools,os
from pathlib import Path
import numpy as np
from a02follow_binding import study
from a02follow_run import write
from a02follow_cone import Service,orders,ulp

def main(out,build):
    out.mkdir(parents=True,exist_ok=False);_,_,case=study.fixture('transfer_force',4);mat=case['joint_system']['joints'][0]['material']
    services={p:Service(build/'source/target'/p/('a02follow_cone_probe'+('.exe' if os.name=='nt' else '')),out/(p+'-raw.json')) for p in ('release','debug')};rows=[]
    try:
        for normal in itertools.permutations([1.,2.,3.]):
            n=np.array(normal);n/=np.linalg.norm(n)
            for axis in range(3):
                a=np.eye(3)[axis];t=np.cross(n,a);t/=np.linalg.norm(t);base=-100*n+30*t
                for k in (-4,-2,-1,0,1,2,4):
                    q=base.copy();q[axis]=ulp(q[axis],k)
                    request=dict(kind='local',material=mat,normal=list(normal),traction=q.tolist(),target=[1e-5,0.,-1e-5])
                    answers={p:s.call(request) for p,s in services.items()}
                    rows.append(dict(normal=list(normal),axis=axis,ulp_steps=k,request=request,orders=orders(q,normal,mat['friction_coefficient']),answers=answers,builds_exact=answers['release']==answers['debug']))
    finally:
        for s in services.values():s.close()
    write(out/'rounding-results.json',dict(passed=False,diagnostic_only=True,rows=rows,search_is_finite_not_error_bound=True,native_gate_unmodified=True))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--build',type=Path,required=True);a=p.parse_args();main(a.output.resolve(),a.build.resolve())
