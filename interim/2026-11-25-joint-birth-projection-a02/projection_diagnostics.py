"""Uncorrected mathematical projections from diagnostic inputs, not 3D solves."""
import argparse
from pathlib import Path
import numpy as np
from study import read,write,sha
from joint_reference import System


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);a=p.parse_args()
    folder=a.run/'projection';mesh=read(folder/'mesh.json');case=read(folder/'case.json')
    s=System(mesh,case['joint_system']);nd=3*len(mesh['nodes'])
    matrix=np.einsum('qag,q->gqa',s.c,np.sqrt(s.weight)).reshape(nd,-1);rows={}
    for name in ('good','tension','outside-circle','balanced-unreproducible-tip'):
        path=folder/(name+'-request.json');before=sha(path);rhs=np.asarray(read(path)['residual_kn'])
        x=np.linalg.lstsq(matrix,rhs,rcond=1e-12)[0];q=x.reshape(-1,3)/np.sqrt(s.weight)[:,None]
        normal=np.einsum('qa,qa->q',q,s.normal);tau=q-normal[:,None]*s.normal
        excess=np.linalg.norm(tau,axis=1)-s.mu*(-normal)
        residual=float(np.max(np.abs(matrix@x-rhs)))
        assert before==sha(path)
        if name=='outside-circle':assert normal.max()<0 and excess.max()>0 and residual<1e-10
        if name=='tension':assert normal.max()>0 and residual<1e-10
        if name=='balanced-unreproducible-tip':assert residual>1e-10
        rows[name]=dict(input_sha256=before,postsolve_mathematical_diagnostic_not_new_3d_solve=True,
            maximum_nodal_reproduction_error_kn=residual,normal_traction_range_kpa=[float(normal.min()),float(normal.max())],
            maximum_friction_excess_kpa=float(excess.max()),uncorrected_traction_kpa=q.tolist(),
            production_accepted=read(folder/(name+'-result.json'))['accepted'])
    write(a.run/'projection-diagnostics.json',rows)


if __name__=='__main__':main()
