"""Independent reduction of existing backward-Euler MCC; no production calls.

Eliminate q and pc analytically; scipy HYBR solves two unknowns rather than
production's analytic 4x4 Newton. Same frozen old-state elastic moduli.
"""
import copy
import numpy as np
from scipy.optimize import root


def step(mat, old, de):
    de=np.asarray(de);stress=np.asarray(old['effective_stress_compression_positive_kpa'])
    p0=float(stress[:3].mean());bulk=(1+old['void_ratio'])*p0/mat['kappa']
    shear=3*bulk*(1-2*mat['poisson'])/(2*(1+mat['poisson']))
    trial=stress.copy();ev=float(de[:3].sum())
    trial[:3]+=bulk*ev+2*shear*(de[:3]-ev/3);trial[3:]+=shear*de[3:]
    pt=float(trial[:3].mean());dev=trial.copy();dev[:3]-=pt
    qt=float(np.sqrt(1.5*(dev[:3]@dev[:3]+2*(dev[3:]@dev[3:]))))
    pc0=old['preconsolidation_pressure_kpa'];m2=mat['critical_state_m']**2
    nxt=copy.deepcopy(old);nxt['void_ratio']=float((1+old['void_ratio'])*np.exp(-ev)-1)
    if pt<=0 or nxt['void_ratio']<=0:raise ArithmeticError('MCC tensile mean or invalid volume; no clipping')
    out=trial;pc=pc0
    if ((qt/pc0)**2+m2*(pt/pc0)*(pt/pc0-1))/m2>1e-11:
        scale=max(pc0,pt,qt);ptn=pt/scale;qtn=qt/scale;cn=pc0/scale
        kn=bulk/scale;gn=shear/scale;a=(1+old['void_ratio'])/(mat['lambda']-mat['kappa'])
        def equations(x):
            p,l=x;c=cn*np.exp(a*(ptn-p)/kn);q=qtn/(1+6*gn*l)
            return [p-ptn+kn*l*m2*(2*p-c),q*q+m2*p*(p-c)]
        solved=root(equations,[ptn,0.],tol=1e-11)
        pn,l=solved.x;cn1=cn*np.exp(a*(ptn-pn)/kn)
        if not np.isfinite(solved.x).all() or pn<=0 or l<0 or np.linalg.norm(equations(solved.x))>=1e-11:
            raise ArithmeticError(('independent two-unknown MCC return failed',solved.message,solved.x.tolist()))
        p=pn*scale;q=qt/(1+6*gn*l);pc=cn1*scale
        out=dev*(q/qt if qt>0 else 0.);out[:3]+=p
        dp=np.empty(6);dp[:3]=((trial[:3]-pt)-(out[:3]-p))/(2*shear)+(pt-p)/(3*bulk)
        dp[3:]=(trial[3:]-out[3:])/shear
        nxt['plastic_strain_compression_positive']=(np.asarray(old['plastic_strain_compression_positive'])+dp).tolist()
    nxt['effective_stress_compression_positive_kpa']=out.tolist();nxt['preconsolidation_pressure_kpa']=float(pc)
    return nxt
