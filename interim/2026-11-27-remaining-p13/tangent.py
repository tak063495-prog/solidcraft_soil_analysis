"""Integrate a directional lower matrix for the closed projection over each Bernstein leaf."""
import numpy as np


def integrate(red, ball, elastic, parent, targets, kn, kt, mu):
    l,h=ball
    v=(kt/2)*(l[:,:2]+h[:,:2])
    length=np.linalg.norm(v,axis=1)
    radius=kt*np.linalg.norm((h[:,:2]-l[:,:2])/2,axis=1)*(1+1e-12)+1e-90
    e=np.divide(v,length[:,None],out=np.tile([1.,0.],(len(l),1)),where=length[:,None]>1e-80)
    qmin=np.maximum(length-radius,0.)*(1-1e-12)
    slide=qmin>mu*kn*(-l[:,2])*(1+1e-12)
    eerr=np.minimum(2.,2*radius/np.maximum(length,1e-80)+1e-12)
    # For both one-sided tangents the perpendicular spring part remains positive.
    # Spend half the normal stiffness only on error in the signed cross term.
    perp=kt*np.minimum(1.,mu*kn*(-h[:,2])/np.maximum(length+radius,1e-80))*(1-1e-12)
    beta=np.where(slide,1.,.5)
    cross_err=np.where(slide,eerr,.5+eerr)
    beta=np.where(cross_err>=1.,0.,beta)
    cross_err=np.minimum(1.,cross_err)
    loss=perp*np.minimum(1.,eerr)+mu**2*kn*cross_err**2/2
    ex,ey=e.T
    coeff=np.column_stack((perp*ey**2-loss,-perp*ex*ey,perp*ex**2-loss,
        np.full(len(l),kn/2),-mu*kn*beta*ex/2,-mu*kn*beta*ey/2))
    coeff[elastic]=[kt,0.,kt,kn,0.,0.]
    A=np.zeros((75,75))
    local={}
    for sl,pos,_,sgn,mass,_,_ in red.blocks:
        vals=(coeff[sl].T@mass).reshape(6,15,15)
        vals[4:]*=sgn
        for idx,i,j in ((0,0,0),(1,0,1),(1,1,0),(2,1,1),(3,2,2),(4,0,2),(4,2,0),(5,1,2),(5,2,1)):
            A[np.ix_(3*pos+i,3*pos+j)]+=vals[idx]
        for k in np.unique(parent[sl][targets[parent[sl]]]):
            sel=parent[sl]==k
            # Sum of all shape products equals 1. Store local integral entries in
            # interface coordinates; the assembled matrix above uses every shape product.
            local[int(k)]=(coeff[sl][sel].T@mass[sel].sum(axis=1)).tolist()
    margins={str(f):float(np.linalg.eigvalsh(G+(A+A.T)/2)[0]-(.1+1e-8*max(kn,kt)))
             for f,G in red.floors.items()}
    lower=max((f for f in red.floors if margins[str(f)]>0),default=0.)
    return lower,margins,local
