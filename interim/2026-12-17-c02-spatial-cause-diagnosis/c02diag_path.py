"""Committed-jump path accounting. No local replay substitutes for equilibrium.

For a straight accepted jump, integrate only its closed fraction. The vector
plastic history should satisfy P=C-tau/kt for this fixed planar Coulomb law.
We measure any defect rather than assuming this identity holds in saved data.
"""
import numpy as np
from c02diag_transfer import norm


def primitives(jump,traction,plastic,normal,kt):
    jump=np.asarray(jump,dtype=float);traction=np.asarray(traction,dtype=float)
    plastic=np.asarray(plastic,dtype=float);normal=np.asarray(normal,dtype=float)
    if jump.shape!=traction.shape or jump.shape!=plastic.shape or jump.ndim!=3:
        raise ValueError('path shape mismatch')
    before=np.concatenate((np.zeros_like(jump[:1]),jump[:-1]))
    g0=np.sum(before*normal,axis=-1);g1=np.sum(jump*normal,axis=-1)
    dg=g1-g0;delta=jump-before;dt=delta-dg[:,:,None]*normal
    closed0=g0<=0;closed1=g1<=0;fraction=np.ones_like(g0)
    fraction[~closed0&~closed1]=0
    opening=closed0&~closed1;closing=~closed0&closed1
    np.divide(-g0,dg,out=fraction,where=opening&(dg!=0))
    np.divide(g1,dg,out=fraction,where=closing&(dg!=0))
    tangential=traction-np.sum(traction*normal,axis=-1)[:,:,None]*normal
    elastic=tangential/np.asarray(kt)[None,:,None]
    closed=np.cumsum(dt*fraction[:,:,None],axis=0)
    defect=plastic-(closed-elastic)
    return dict(closed_path=closed,elastic_label=elastic,defect=defect,fraction=fraction,open=g1>0)


def path_summary(parts,plastic):
    return dict(defect_absolute_l2=norm(parts['defect']),defect_maximum=float(np.max(abs(parts['defect']))),closed_path_absolute_l2=norm(parts['closed_path']),elastic_label_absolute_l2=norm(parts['elastic_label']),raw_plastic_absolute_l2=norm(plastic),opening_or_closing_segments=int(((parts['fraction']>0)&(parts['fraction']<1)).sum()))


def compare_primitives(coarse,fine,pcoarse,pfine):
    dp=pfine-pcoarse;dc=fine['closed_path']-coarse['closed_path'];de=fine['elastic_label']-coarse['elastic_label']
    mismatch=coarse['open']!=fine['open'];n=norm(dp)
    return dict(plastic_time_difference_l2=n,closed_path_time_difference_l2=norm(dc),elastic_label_time_difference_l2=norm(de),closed_path_over_plastic=norm(dc)/max(n,1e-12),elastic_label_over_plastic=norm(de)/max(n,1e-12),accounting_residual_l2=norm(dp-(dc-de)),open_label_disagreements=int(mismatch.sum()),plastic_difference_l2_at_label_disagreements=norm(dp[mismatch]))
