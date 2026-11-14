"""Input-only independent reference variants; callable only by explicit diagnosis.

No mortar, no production configuration field, no saved history transfer.
"""
import copy
import numpy as np


def alternate_system(original, rule, interpolation):
    if rule not in ('jacobi4','legendre4') or interpolation not in ('tri6','corner3'):
        raise ValueError('unknown contact observation variant')
    result=copy.copy(original)
    points=[]
    mesh=original.mesh
    xyz=np.asarray(mesh['nodes'])
    from fem_reference import jacobi
    if rule=='jacobi4':
        us,uw=jacobi(4,1);vs,vw=jacobi(4,0)
        rule_rs=[(np.array([u,(1-u)*v]),wu*wv) for u,wu in zip(us,uw) for v,wv in zip(vs,vw)]
    else:
        us,uw=np.polynomial.legendre.leggauss(4);us=(us+1)/2;uw=uw/2
        rule_rs=[(np.array([u,(1-u)*v]),wu*wv*(1-u)) for u,wu in zip(us,uw) for v,wv in zip(us,uw)]
    seen=set()
    for old in original.points:
        identity=(old['joint_id'],old['face_index'])
        if identity in seen:
            continue
        seen.add(identity)
        pa=xyz[old['nodes_a']]
        if len(pa)!=6:
            raise ValueError('bounded affine Tri6 diagnostic only')
        for rs,w in rule_rs:
            l=np.r_[1-rs.sum(),rs];dl=np.array([[-1.,-1.],[1.,0.],[0.,1.]])
            shape=np.r_[l*(2*l-1),[4*l[a]*l[b] for a,b in ((0,1),(1,2),(2,0))]]
            gradient=np.vstack(((4*l-1)[:,None]*dl,[4*(l[a]*dl[b]+l[b]*dl[a]) for a,b in ((0,1),(1,2),(2,0))]))
            jac=(pa-pa[0]).T@gradient
            area=np.cross(jac[:,0],jac[:,1])*w
            weight=np.linalg.norm(area)
            kinematic=shape if interpolation=='tri6' else np.r_[l,[0.,0.,0.]]
            points.append({**old,'position_m':(pa[0]+shape@(pa-pa[0])).tolist(),'normal':(area/weight).tolist(),'weight_m2':float(weight),'shape':kinematic.tolist()})
    result.points=points
    result.normal=np.asarray([p['normal'] for p in points]);result.weight=np.asarray([p['weight_m2'] for p in points])
    for name,key in [('kn','normal_stiffness_kpa_m'),('kt','tangential_stiffness_kpa_m'),('mu','friction_coefficient')]:
        setattr(result,name,np.asarray([p['material'][key] for p in points]))
    result.c=np.zeros((len(points),3,xyz.size))
    for i,p in enumerate(points):
        for a,b,w in zip(p['nodes_a'],p['nodes_b'],p['shape'],strict=True):
            for k in range(3):
                result.c[i,k,3*a+k]-=w;result.c[i,k,3*b+k]+=w
    return result


def water_factory(default_water,rule,interpolation):
    class C02DiagnosticWater(default_water):
        def __init__(self,mesh,ground,case,layout=None):
            if layout is not None:
                raise ValueError('do not import production layout/history into independent variant')
            super().__init__(mesh,ground,case)
            self.contact=alternate_system(self.contact,rule,interpolation)
            self.faces=list({(p['joint_id'],p['face_index']):p for p in self.contact.points}.values())
    return C02DiagnosticWater
