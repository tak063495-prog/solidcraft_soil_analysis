"""Forced fine-grid open intervals for a DIFFERENT equilibrium problem.

No change to production or old reference modules. A forced-open point carries
no force; a forced-closed point uses a bilateral normal spring and Coulomb cap
mu*max(-qn,0). The bilateral branch (including tensile qn if the control opens)
is explicit: this is an intervention, never an accepted original solution.
Event cuts only integrate the prescribed constitutive time history; no new
load-step controller, tolerances, bisection, or global acceptance rule.
"""
from math import fsum
import numpy as np


def open_at(schedule,stage,time):
    if time==1. and 'endpoint_open' in schedule:
        return np.asarray(schedule['endpoint_open'][stage])>0
    return np.array([any(a<time<b for a,b in intervals) for intervals in schedule['intervals'][stage]],dtype=bool)


def project_forced(system,start,old_q,target,counts,schedule,stage,start_progress,end_progress):
    if target.ndim==2:target=target[:,None,:]
    n=system.normal[:,None,:];kn=system.kn[:,None];kt=system.kt[:,None];mu=system.mu[:,None]
    q=np.broadcast_to(old_q[:,None,:],target.shape).copy();x=np.broadcast_to(start[:,None,:],target.shape).copy()
    flow=np.zeros_like(q);acc=np.zeros(target.shape[:2]);pd=np.zeros(target.shape[:2]);nd=pd.copy();work=pd.copy()
    h=end_progress-start_progress
    if h<=0:raise ValueError('forced timeline requires increasing time')
    # Each GP retains exactly its frozen local count plus the prescribed open
    # event roots. Events are common to all Newton iterates and perturbations.
    for gp,count in enumerate(counts):
        cuts={i/int(count) for i in range(int(count)+1)}
        for a,b in schedule['intervals'][stage][gp]:
            for t in (a,b):
                if start_progress<t<end_progress:cuts.add((t-start_progress)/h)
        tau=q[gp]-(q[gp]*n[gp]).sum(axis=1)[:,None]*n[gp]
        def energy(v):
            vn=(v*n[gp]).sum(axis=1);vt=v-vn[:,None]*n[gp]
            return vn*vn/(2*kn[gp])+(vt*vt).sum(axis=1)/(2*kt[gp])
        previous=0.
        for b in sorted(cuts)[1:]:
            end=start[gp]+(target[gp]-start[gp])*b
            opened=any(a<start_progress+h*(previous+b)/2<c for a,c in schedule['intervals'][stage][gp])
            dx=end-x[gp];before=q[gp].copy()
            if opened:
                q[gp]=0.;tau[:]=0.
                # Removed spring energy is intervention work, not plastic D.
                nd[gp]+=energy(before)
            else:
                dt=dx-(dx*n[gp]).sum(axis=1)[:,None]*n[gp]
                trial=tau+kt[gp,:,None]*dt
                qn=kn[gp]*(end*n[gp]).sum(axis=1) # bilateral forced branch
                cap=mu[gp]*np.maximum(-qn,0.)
                length=np.sqrt((trial*trial).sum(axis=1));scale=np.minimum(1.,np.divide(cap,length,out=np.ones_like(cap),where=length>0))
                tau=trial*scale[:,None];dp=(trial-tau)/kt[gp,:,None]
                q[gp]=qn[:,None]*n[gp]+tau
                dd=np.maximum((tau*dp).sum(axis=1),0.);dw=(q[gp]*dx).sum(axis=1)
                flow[gp]+=dp;acc[gp]+=np.sqrt((dp*dp).sum(axis=1));pd[gp]+=dd;work[gp]+=dw
                nd[gp]+=dw-(energy(q[gp])-energy(before))-dd
            x[gp]=end;previous=b
    if target.shape[1]==1:return q,flow,acc,pd[:,0],nd[:,0],work[:,0]
    return q,flow,acc,pd,nd,work
