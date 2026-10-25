"""Independent compression-positive principal-stress checks (NumPy eigvalsh)."""
import numpy as np


def strength(p):
    c=p['cohesion_kpa'];phi=np.deg2rad(p['friction_deg']);psi=np.deg2rad(p['dilation_deg'])
    if p['flow']=='davis':
        beta=np.cos(phi)*np.cos(psi)/(1-np.sin(phi)*np.sin(psi));c*=beta;phi=np.arctan(beta*np.tan(phi));psi=phi
    return c,phi,psi


def tensor(s):
    s=np.asarray(s);out=np.zeros((*s.shape[:-1],3,3))
    out[...,0,0]=s[...,0];out[...,1,1]=s[...,1];out[...,2,2]=s[...,2]
    out[...,0,1]=out[...,1,0]=s[...,3];out[...,1,2]=out[...,2,1]=s[...,4];out[...,2,0]=out[...,0,2]=s[...,5]
    return out


def normalized_yield(p,s):
    values=np.linalg.eigvalsh(tensor(s));c,phi,_=strength(p);scale=np.maximum(np.max(np.abs(values),axis=-1),max(p['cohesion_kpa'],1.))
    result=((1-np.sin(phi))*values[...,-1]-(1+np.sin(phi))*values[...,0]-2*c*np.cos(phi))/scale
    if p.get('tension_cutoff_kpa') is not None:result=np.maximum(result,(-values[...,0]-p['tension_cutoff_kpa'])/scale)
    return result


def step(p,old,de):
    """Iterative principal-plane complementarity (projected Gauss-Seidel).

    Rust enumerates nonsingular active subsets with LU. This reference instead
    iterates all multipliers and checks complementarity; no Rust data/tangent is
    used as a reference state. It is intentionally limited to small test paths.
    """
    e=p['young_kpa'];nu=p['poisson'];g=e/(2*(1+nu));lam=e*nu/((1+nu)*(1-2*nu));d=np.full((3,3),lam);np.fill_diagonal(d,lam+2*g)
    de=np.asarray(de);ds=np.r_[d@de[:3],g*de[3:]];trial=np.asarray(old['effective_stress_compression_positive_kpa'])+ds
    values,v=np.linalg.eigh(tensor(trial));values=values[::-1];v=v[:,::-1];c,phi,psi=strength(p)
    n=[];m=[];limit=[]
    for i in range(3):
        for j in range(3):
            if i==j:continue
            a=np.zeros(3);b=a.copy();a[i]=1-np.sin(phi);a[j]=-1-np.sin(phi);b[i]=1-np.sin(psi);b[j]=-1-np.sin(psi)
            n.append(a);m.append(b);limit.append(2*c*np.cos(phi))
    if p.get('tension_cutoff_kpa') is not None:
        for i in range(3):a=np.zeros(3);a[i]=-1;n.append(a);m.append(a);limit.append(p['tension_cutoff_kpa'])
    n=np.asarray(n);m=np.asarray(m);limit=np.asarray(limit);dm=d@m.T;den=np.einsum('ij,ji->i',n,dm)
    multipliers=np.zeros(len(n));stress=values.copy();scale=max(np.max(np.abs(values)),p['cohesion_kpa'],1.)
    for iterations in range(20000):
        for i in range(len(n)):
            change=max(0.,multipliers[i]+(n[i]@stress-limit[i])/den[i])-multipliers[i]
            multipliers[i]+=change;stress-=dm[:,i]*change
        f=n@stress-limit
        if np.max(f)<2e-11*scale and np.max(np.abs(multipliers*f))<1e-13*scale:break
    else:raise ArithmeticError('independent projected complementarity iteration did not converge')
    assert np.min(multipliers)>=0 and np.all(stress[:-1]>=stress[1:]-1e-8*scale)
    dp=np.linalg.solve(d,values-stress);plastic=v@np.diag(dp)@v.T;full=v@np.diag(stress)@v.T
    voigt=lambda a:np.array([a[0,0],a[1,1],a[2,2],a[0,1],a[1,2],a[2,0]])
    pv=voigt(plastic);pv[3:]*=2;dev=plastic-np.eye(3)*np.trace(plastic)/3
    out={'effective_stress_compression_positive_kpa':voigt(full).tolist(),
        'plastic_strain_compression_positive':(np.asarray(old['plastic_strain_compression_positive'])+pv).tolist(),
        'void_ratio':float((1+old['void_ratio'])*np.exp(-de[:3].sum())-1),
        'equivalent_plastic_strain':float(old['equivalent_plastic_strain']+np.sqrt(2/3*np.square(dev).sum())),
        'plastic_work_kpa':float(old['plastic_work_kpa']+max(0,stress@dp))}
    assert normalized_yield(p,out['effective_stress_compression_positive_kpa'])<1e-8
    return out


def advance(p,old,de,maximum):
    count=max(1,int(np.ceil(np.max(np.abs(de))/maximum)));state=old
    for _ in range(count):state=step(p,state,np.asarray(de)/count)
    return state
