"""Continuous RM manufactured load on an exact graph with its analytic director.

Differentiate the continuum virtual-work fluxes in surface parameters. Domain
force/couple and conormal flux on the three natural edges use a separate GL rule.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from standard import connectivity, execute

E, NU, T, AMP = 20000., .25, .1, 1e-4
G = E/(2*(1+NU))
RULE, WEIGHT = leggauss(6)
RULE, WEIGHT = (RULE+1)/2, WEIGHT/2


def outer(a,b):
    return a[..., :, None]*b[..., None, :]


def geometry(p,a,b):
    r,s=p.T
    zero=np.zeros_like(r)
    x=np.array([r,s,a*r*r+b*s*s]).T
    normal=np.stack([-2*a*r,-2*b*s,1+zero],axis=1)
    length=np.sqrt(np.sum(normal*normal,axis=1))
    d=normal/length[:,None]
    dd=np.stack([(v-d*(d@v)[:,None])/length[:,None] for v in (np.array([-2*a,0,0]),np.array([0,-2*b,0]))],axis=-1)
    dx=np.stack([np.array([1+zero,zero,2*a*r]).T,np.array([zero,1+zero,2*b*s]).T],axis=-1)
    return x,dx,d,dd


def exact(p,a,b):
    r,s=p.T
    x,dx,d,dd=geometry(p,a,b)
    u=AMP*np.stack([np.sin(np.pi*r)*np.sin(np.pi*s),r*r*s,r*s*s],axis=1)
    du=AMP*np.stack([np.array([np.pi*np.cos(np.pi*r)*np.sin(np.pi*s),2*r*s,s*s]).T,
                     np.array([np.pi*np.sin(np.pi*r)*np.cos(np.pi*s),r*r,2*r*s]).T],axis=-1)
    basis=np.concatenate([dx,d[:,:,None]],axis=-1)
    grad=du@np.linalg.inv(basis)[:,:2,:]
    curl=np.array([grad[:,2,1]-grad[:,1,2],grad[:,0,2]-grad[:,2,0],grad[:,1,0]-grad[:,0,1]]).T
    trial=AMP*np.array([r*r+s,r*s,s*s]).T
    theta=trial+d*np.sum(d*(curl/2-trial),axis=1)[:,None]
    return u,du,theta


def derivative(fn,p,k):
    # Propagate the selected parameter tangent in the imaginary component.
    # All continuum operations, including normalization, are analytic here.
    shift=np.eye(2)[k]*1e-24j
    return fn(p+shift).imag/1e-24


def continuum(p,z,a,b):
    _,dx,d,dd=geometry(p,a,b)
    _,du,_=exact(p,a,b)
    r,s=p.T
    trial=AMP*np.column_stack((r*r+s,r*s,s*s))
    zero=np.zeros_like(r)
    dtrial=(AMP*np.column_stack((2*r,s,zero)),AMP*np.column_stack((zero+1,r,2*s)))
    # The drilling term of theta is parallel to d and cancels identically.
    delta=np.cross(trial,d)
    deriv=np.stack([np.cross(dtrial[k],d)+np.cross(trial,dd[:,:,k]) for k in range(2)],axis=-1)
    basis=np.concatenate([dx+z*dd,d[:,:,None]],axis=2)
    inv=np.linalg.inv(basis)
    gradient=np.concatenate([du+z*deriv,delta[:,:,None]],axis=2)@inv
    eps=(gradient+gradient.swapaxes(1,2))/2
    proj=np.eye(3)-outer(d,d)
    tangent=proj@eps@proj
    shear=2*(proj@eps@d[:,:,None])[:,:,0]
    sigma=2*G*tangent+proj*(E*NU/(1-NU**2)*np.trace(tangent,axis1=1,axis2=2))[:,None,None]
    sigma+=5*G/6*(outer(shear,d)+outer(d,shear))
    piola=np.linalg.det(basis)[:,None,None]*sigma@inv.swapaxes(1,2)
    return eps,piola


def fluxes(p,a,b):
    _,_,d,dd=geometry(p,a,b)
    flux=np.zeros((len(p),6,2),dtype=p.dtype)
    direct=np.zeros((len(p),3),dtype=p.dtype)
    for unit,w in zip(RULE,WEIGHT):
        z=T*(unit-.5)
        _,stress=continuum(p,z,a,b)
        flux[:,:3]+=T*w*stress[:,:,:2]
        direct+=T*w*np.cross(d,stress[:,:,2])
        for k in range(2):
            flux[:,3:,k]+=T*w*z*np.cross(d,stress[:,:,k])
            direct+=T*w*z*np.cross(dd[:,:,k],stress[:,:,k])
    return flux,direct


def shape(r,s):
    l=1-r-s
    return np.array([l*(2*l-1),r*(2*r-1),s*(2*s-1),4*l*r,4*r*s,4*s*l])


def mesh_case(n,a,b,distorted):
    params=np.array([[i/(2*n),j/(2*n)] for j in range(2*n+1) for i in range(2*n+1)])
    for j in range(1,2*n):
        for i in range(1,2*n):
            if distorted:
                # Perturb macro vertices, then put midsides at parametric midpoints.
                if i%2==0 and j%2==0:
                    params[j*(2*n+1)+i,0]+=.12*(-1)**((i+j)//2)/n
    faces=connectivity(n)
    for f in faces:
        v=f['nodes']
        for k in range(3): params[v[3+k]]=(params[v[k]]+params[v[(k+1)%3]])/2
    x,_,d,_=geometry(params,a,b)
    nodes=len(x)
    load=np.zeros((nodes,6))
    ids=np.array([f['nodes'] for f in faces])
    tri=params[ids[:,:3]]
    jac=np.stack([tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]],axis=2)
    det=np.linalg.det(jac)
    # Integrate continuous body and distributed couple against nodal test fields.
    ref=np.array([[r,(1-r)*v] for r in RULE for v in RULE])
    weights=np.array([wr*wv*(1-r) for r,wr in zip(RULE,WEIGHT) for wv in WEIGHT])
    samples=tri[:,None,0,:]+np.einsum('eij,qj->eqi',jac,ref)
    p=samples.reshape(-1,2)
    body=-sum(derivative(lambda t:fluxes(t,a,b)[0][:,:,k],p,k) for k in range(2))
    body[:,3:]+=fluxes(p,a,b)[1]
    local=np.einsum('eqc,iq,q,e->eic',body.reshape(len(tri),len(ref),6),shape(*ref.T),weights,det)
    np.add.at(load,ids.ravel(),local.reshape(-1,6))
    for v in ids:
        for k in range(3):
            h,j=v[k],v[(k+1)%3]
            start,end=params[h],params[j]
            if not (np.allclose([start[0],end[0]],1.) or start[1]==end[1] in (0.,1.)):
                continue
            edge=end-start
            p=start+RULE[:,None]*edge
            traction=fluxes(p,a,b)[0]@np.array([edge[1],-edge[0]])
            nedge=np.array([(1-RULE)*(1-2*RULE),RULE*(2*RULE-1),4*RULE*(1-RULE)])
            load[[h,j,v[k+3]]]+=np.einsum('iq,qc,q->ic',nedge,traction,WEIGHT)
    # Physical nodal rotations tilt their director. Their axial component
    # is the independent drilling unknown and has no applied MMS couple.
    load[:,3:]-=d*np.sum(load[:,3:]*d,axis=1)[:,None]
    u,_,theta=exact(params,a,b)
    q=np.hstack((u,theta))
    fixed=np.repeat(params[:,0]==0.,6)
    prescribed=np.where(fixed,q.ravel(),0.)
    case=dict(nodes=x.tolist(),faces=faces,directors=[dict(group='shell',node=i,direction=z.tolist()) for i,z in enumerate(d)],
        section=dict(id='mms',young_kpa=E,poisson=NU,thickness_m=T,unit_weight_kn_m3=0.,drilling_factor=1e-4),
        steps=[dict(event='install',holder='edge',contacts={}),
               dict(event='load',prescribed_configuration=prescribed.tolist(),free_dofs=np.flatnonzero(~fixed).tolist(),
                    nodal_load=load.ravel().tolist(),plate_traction_kpa=[0.,0.,0.])])
    return case,params


def errors(frame,case,params,a,b):
    total=np.zeros(5)
    q=np.array(frame['configuration']).reshape(-1,6)
    h=0.
    for element in frame['active']['elements']:
        ids=element['face']['nodes']
        xx=np.array(case['nodes'])[ids]
        h=max(h,max(np.linalg.norm(x-y) for x in xx for y in xx))
        fields=element['fields']
        local=np.array([f['rs'] for f in fields[:6]])
        N=shape(*local.T).T
        p=N@params[ids]
        u,_,theta=exact(p,a,b)
        err=N@q[ids]-np.hstack([u,theta])
        _,dx,_,_=geometry(p,a,b)
        surface=np.linalg.norm(np.cross(dx[:,:,0],dx[:,:,1]),axis=1)
        delta=params[ids[1:3]]-params[ids[0]]
        area=surface*abs(np.linalg.det(delta))*np.array([.109951743655322]*3+[.223381589678011]*3)/2
        total[:2]+=[np.sum(area*np.sum(err[:,sl]**2,axis=1)) for sl in (slice(3),slice(3,6))]
        mid=None
        for layer in [1,0,2]:
            ff=fields[6*layer:6*layer+6]
            z=ff[0]['z_m']
            eps,_=continuum(p,z,a,b)
            axes=np.array([f['material_axes'] for f in ff])
            eps=axes@eps@axes.swapaxes(1,2)
            elastic_error=np.array([f['strain'] for f in ff])-np.array([eps[:,0,0],eps[:,1,1],2*eps[:,0,1]]).T
            shear_error=np.array([f['shear_strain'] for f in ff])-2*eps[:,:2,2]
            if layer==1: mid=elastic_error
            weight=np.array([f['weight'] for f in ff])
            def energy(v):
                return E/(1-NU**2)*(v[:,0]**2+v[:,1]**2+2*NU*v[:,0]*v[:,1]+(1-NU)*v[:,2]**2/2)
            total[2:]+=[np.sum(weight*energy(mid)),np.sum(weight*energy(elastic_error-mid)),
                         np.sum(weight*5*G/6*np.sum(shear_error**2,axis=1))]
    return h,dict(zip(['u_L2','rotation_L2','membrane_energy_norm','bending_energy_norm','shear_energy_norm'],np.sqrt(total).tolist()))


def continuum_work(a,b):
    # Independent tensor GL integration on the exact square and thickness.
    x,w=leggauss(12)
    x,w=(x+1)/2,w/2
    p=np.array([[r,s] for r in x for s in x])
    area=np.outer(w,w).ravel()
    _,dx,d,dd=geometry(p,a,b)
    energy=0.
    for t,wt in zip(x,w):
        z=T*(t-.5)
        eps,piola=continuum(p,z,a,b)
        basis=np.concatenate([dx+z*dd,d[:,:,None]],axis=2)
        # P A^T = J sigma for the reference chart.
        energy+=T*wt*np.sum(area*np.einsum('nij,nij->n',eps,piola@basis.swapaxes(1,2)))/2
    edge=np.column_stack((np.zeros_like(x),x))
    force=-np.einsum('q,qi->i',w,fluxes(edge,a,b)[0][:,:3,0])
    return float(energy),force


def extrapolate(series):
    """Fit error = limit + C h**p to three monotonically decreasing errors."""
    if len(series)!=3 or any('errors' not in row for row in series):
        return None
    h0,h1,h2=[row['h'] for row in series]
    result={}
    for name in series[0]['errors']:
        a,b,c=[row['errors'][name] for row in series]
        fit=None
        if a>b>c>0 and h0>h1>h2:
            ratio=(a-b)/(b-c)
            def equation(p):
                return np.expm1(p*np.log(h0/h1))/(-np.expm1(p*np.log(h2/h1)))-ratio
            if equation(1e-6)*equation(20.)<0:
                p=float(brentq(equation,1e-6,20.))
                fit=dict(exponent=p,extrapolated_error=float(c-(b-c)/np.expm1(p*np.log(h1/h2))))
        result[name]=fit
    return result


def run(cli,out):
    out.mkdir(parents=True,exist_ok=False)
    rows=[]
    for a,b in [(0.,0.),(.18,0.),(.18,-.12)]:
        work,edge_force=continuum_work(a,b)
        for distorted in [False,True]:
            prev=None
            for n in [2,4,8]:
                case,p=mesh_case(n,a,b,distorted)
                frame,error=execute(cli,out/f'mms-{a}-{b}-{distorted}-{n}.json',case)
                row=dict(a=a,b=b,distorted=distorted,n=n,diagnostic=error,passed=False)
                if frame:
                    h,norm=errors(frame,case,p,a,b)
                    row.update(h=h,errors=norm,interval_order=None if prev is None else
                        {k:float(np.log(prev['errors'][k]/val)/np.log(prev['h']/h)) for k,val in norm.items()},
                        applied_work=frame['cumulative_applied_work_kn_m'],
                        reaction_work=frame['cumulative_reaction_work_kn_m'],
                        free_residual=float(np.linalg.norm(np.array(frame['reactions'])[case['steps'][1]['free_dofs']])))
                    reactions=np.reshape(frame['reactions'],(-1,6))[p[:,0]==0,:3].sum(axis=0)
                    row.update(continuum_work=work,work_difference=row['applied_work']+row['reaction_work']-work,
                               exact_support_force=edge_force.tolist(),support_force=reactions.tolist(),
                               support_force_error=float(np.linalg.norm(reactions-edge_force)))
                    prev=row
                rows.append(row)
                if n==8: row['three_level_fit_diagnostic']=extrapolate(rows[-3:])
                (out/'measurements.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')
                print(a,b,distorted,n,row.get('errors'),error,flush=True)
    return rows


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('cli',type=Path)
    parser.add_argument('folder',type=Path)
    args=parser.parse_args()
    run(args.cli.resolve(),args.folder)
