"""Dry 2D slice equilibrium, implemented here from force/moment balance.

Units kN,m,kPa. See README for derivation and checked sources. No solver imports.
"""
import itertools
import numpy as np
from scipy.optimize import brentq, root


def ground(x):
    return np.minimum(10., (32.-np.asarray(x))/2.)


def arc(a, b, sag):
    p=np.array([a,float(ground(a))]); q=np.array([b,float(ground(b))])
    chord=q-p; length=np.linalg.norm(chord)
    normal=np.array([-chord[1],chord[0]])/length
    center=(p+q)/2+normal*(length*length/(8*sag)-sag/2)
    radius=length*length/(8*sag)+sag/2
    if center[1] <= max(p[1],q[1]):
        raise ValueError('endpoints must lie on lower circle branch')
    def bottom(x):
        return center[1]-np.sqrt(np.maximum(0., radius**2-(np.asarray(x)-center[0])**2))
    return bottom, center, radius


def slices(a,b,n,bottom):
    x=np.linspace(a,b,n+1); y=bottom(x)
    # Validate densely, including endpoints; no slip below fixed base or above terrain.
    probe=np.unique(np.r_[np.linspace(a,b,10*n+1),12.])
    probe=probe[(probe>=a)&(probe<=b)]
    if (np.any(bottom(probe)<-1e-10) or
            np.any(bottom(probe)>ground(probe)+1e-9) or
            abs(y[0]-ground(a))>1e-8 or abs(y[-1]-ground(b))>1e-8):
        raise ValueError('inadmissible slip geometry')
    dx=np.diff(x); dy=np.diff(y); length=np.hypot(dx,dy)
    # Piecewise linear polygonal approximation to arc; 2-point Gauss integration
    # split at terrain corner integrates polygon area and weight first moment exactly.
    weights=[]; xg=[]
    for i in range(n):
        cuts=sorted(set([x[i],x[i+1]]+([12.] if x[i]<12<x[i+1] else [])))
        area=moment=0.
        for left,right in zip(cuts[:-1],cuts[1:]):
            z=(left+right)/2+(right-left)/2*np.array([-1,1])/np.sqrt(3)
            base=y[i]+(z-x[i])*dy[i]/dx[i]
            h=ground(z)-base
            area+=(right-left)/2*np.sum(h)
            moment+=(right-left)/2*np.sum(z*h)
        if area<=0: raise ValueError('empty slice')
        weights.append(20*area); xg.append(moment/area)
    return dict(x=x,y=y,b=dx,l=length,s=-dy/length,c=dx/length,
                W=np.array(weights),xg=np.array(xg),xb=(x[:-1]+x[1:])/2,
                yb=(y[:-1]+y[1:])/2)


def bishop(s):
    """USACE EM1110-2-1902 (2003), C-15/C-16, dry circular case only."""
    w,sa,ca,b=s['W'],s['s'],s['c'],s['b']; mu=np.tan(np.deg2rad(20))
    def residual(f):
        m=ca+sa*mu/f
        if np.min(m)<=0: return np.nan
        return np.sum((10*b+w*mu)/m)-f*np.sum(w*sa)
    # Bracket excludes singular m_alpha; do not turn a failed root into an Fs.
    lower=max(.1,float(np.max(-sa*mu/ca))+1e-6)
    if residual(lower)*residual(10)>=0: raise ValueError('no Bishop root')
    f=brentq(residual,lower,10,xtol=1e-11)
    normal=(w-10*s['l']*sa/f)/(ca+mu*sa/f)
    return dict(fs=float(f),residual=float(abs(residual(f))/sum(w)),
                compression_only=bool(np.min(normal)>=-1e-7),min_normal=float(np.min(normal)))


def propagate(s, f, lam, function='constant'):
    """Slice force balance, X_j=lam*f_j*E_j. Global moment about origin."""
    n=len(s['W']); t=np.linspace(0,1,n+1)
    fun=np.ones(n+1) if function=='constant' else np.sin(np.pi*t)
    mu=np.tan(np.deg2rad(20)); E=[0.]; normals=[]; basal=[]
    for i in range(n):
        sa,ca,w,C=s['s'][i],s['c'][i],s['W'][i],10*s['l'][i]
        aa=sa-mu*ca/f; bb=-C*ca/f
        den=ca+mu*sa/f-lam*fun[i+1]*aa
        if abs(den)<1e-8: raise ValueError('singular slice balance')
        N=(w-C*sa/f-lam*(fun[i]-fun[i+1])*E[-1]+lam*fun[i+1]*bb)/den
        T=(C+mu*N)/f
        E.append(E[-1]+N*sa-T*ca); normals.append(N)
        basal.append([N*sa-T*ca,N*ca+T*sa])
    E=np.array(E); basal=np.array(basal); N=np.array(normals)
    moments=s['xb']*basal[:,1]-s['yb']*basal[:,0]-s['xg']*s['W']
    residual=np.array([E[-1]/sum(s['W']),sum(moments)/(32*sum(s['W']))])
    # Local moment recurrence gives the line of thrust; retain admissibility.
    moment_Ey=0.; heights=[]
    for i in range(n):
        moment_Ey += -moments[i]-s['x'][i]*lam*fun[i]*E[i]+s['x'][i+1]*lam*fun[i+1]*E[i+1]
        heights.append(moment_Ey/E[i+1] if abs(E[i+1])>1e-7 else np.nan)
    return residual,E,N,np.array(heights)


def complete(s, function='constant'):
    def residual(z):
        try: return propagate(s,*z,function)[0]
        except (ValueError, FloatingPointError): return np.array([1e3,1e3])
    sol=root(residual,[1.4,.3],tol=1e-9)
    f,lam=sol.x
    r,E,N,heights=propagate(s,f,lam,function)
    if not sol.success or np.linalg.norm(r,np.inf)>1e-8 or not .1<f<10 or abs(lam)>2:
        raise ValueError('equilibrium root unresolved')
    interior=heights[:-1]; base=s['y'][1:-1]; top=ground(s['x'][1:-1])
    thrust=bool(np.all(np.isfinite(interior)) and np.all(interior>=base-1e-6) and np.all(interior<=top+1e-6))
    compression_only=bool(np.min(E)>=-1e-6 and np.min(N)>=-1e-6 and thrust)
    return dict(fs=float(f),lambda_value=float(lam),residual=float(np.linalg.norm(r,np.inf)),
                compression_only=compression_only,thrust_inside_mass=thrust,min_interslice=float(np.min(E)),min_normal=float(np.min(N)))


def search(grid, n, noncircular=False, bounds=((2.,11.5),(20.,32.),(1.,7.))):
    """Full finite grid, no extrapolated/global-minimum claim.

    Noncircular candidates are arc + amplitude*sin(pi*t)^2; endpoints unchanged.
    """
    methods=['spencer','morgenstern_price'] if noncircular else ['bishop','spencer','morgenstern_price']
    best={k:None for k in methods}; compression_best={k:None for k in methods}; rejected={k:0 for k in methods}; invalid=0; count=0
    for a,b,sag in itertools.product(*(np.linspace(lo,hi,grid) for lo,hi in bounds)):
        for amplitude in ([-.5,0.,.5] if noncircular else [0.]):
            count+=1
            try:
                bot,_,_=arc(a,b,sag)
                def surface(x):
                    return bot(x)+amplitude*np.sin(np.pi*(np.asarray(x)-a)/(b-a))**2
                sl=slices(a,b,n,surface)
            except ValueError:
                invalid+=1; continue
            for method in methods:
                try:
                    v=bishop(sl) if method=='bishop' else complete(sl,'constant' if method=='spencer' else 'sine')
                    v['surface']=[float(a),float(b),float(sag),float(amplitude)]
                    if best[method] is None or v['fs']<best[method]['fs']: best[method]=v
                    if v['compression_only'] and (compression_best[method] is None or v['fs']<compression_best[method]['fs']): compression_best[method]=v
                except ValueError: rejected[method]+=1
    return dict(grid=grid,slices=n,bounds=[list(pair) for pair in bounds],noncircular=noncircular,candidates=count,invalid_geometry=invalid,
                unresolved_roots=rejected,minima=best,compression_only_minima=compression_best,accepted=None)
