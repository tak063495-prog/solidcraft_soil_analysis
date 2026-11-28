"""Exact frozen contact arithmetic with optional read-only branch operands."""
import numpy as np

def contact(system,old_jump,old_traction,target,count,trace=None):
    """Independent two-tangent Coulomb disk projection along a straight jump."""
    if target.ndim==2:target=target[:,None,:]
    n=system.normal[:,None,:];kn=system.kn[:,None];kt=system.kt[:,None];mu=system.mu[:,None]
    delta=(target-old_jump[:,None,:])/count;dt=delta-(delta*n).sum(axis=2)[:,:,None]*n
    tau=old_traction[:,None,:]-(old_traction*system.normal).sum(axis=1)[:,None,None]*n
    tau=np.broadcast_to(tau,target.shape).copy();slip=np.zeros_like(tau);cumulative=np.zeros(target.shape[:2])
    gn=(old_jump*system.normal).sum(axis=1)[:,None];dg=(delta*n).sum(axis=2)
    for step in range(count):
        next_g=gn+dg;was=gn<=0;closed=next_g<=0
        fraction=np.ones_like(next_g)
        opening=was & ~closed;closing=~was & closed
        fraction[~was & ~closed]=0
        np.divide(-gn,dg,out=fraction,where=opening & (dg!=0))
        np.divide(next_g,dg,out=fraction,where=closing & (dg!=0))
        tau=np.where((~was)[:,:,None],0,tau)
        trial=tau+kt[:,:,None]*dt*fraction[:,:,None]
        pressure=kn*np.maximum(-next_g,0);limit=mu*pressure
        norm=np.linalg.norm(trial,axis=2);scale=np.minimum(1,np.divide(limit,norm,out=np.ones_like(norm),where=norm>0))
        projected=trial*scale[:,:,None];flow=(trial-projected)/kt[:,:,None]
        if trace is not None:
            trace.append(dict(trial_norm_kpa=norm.copy(),radius_kpa=limit.copy(),flow_m=flow.copy(),gap_m=next_g.copy(),closed=closed.copy()))
        slip+=flow;cumulative+=np.linalg.norm(flow,axis=2)
        tau=np.where(closed[:,:,None],projected,0);gn=next_g
    traction=tau+kn[:,:,None]*np.minimum(gn,0)[:,:,None]*n
    return traction,slip,cumulative
