"""Independent polynomial Tri6 Mindlin reference.

Quadratic Lagrange polynomials are solved from physical nodal coordinates.
Shear edge moments and interior means are integrated with GL5; the full energy
uses a separate 5x5 Duffy rule (Rust uses point tying and Dunavant integration).
"""
import numpy as np
from plate_reference import polynomial

SAMPLES=((.816847572980459,.091576213509771,.091576213509771),
 (.091576213509771,.816847572980459,.091576213509771),
 (.091576213509771,.091576213509771,.816847572980459),
 (.108103018168070,.445948490915965,.445948490915965),
 (.445948490915965,.108103018168070,.445948490915965),
 (.445948490915965,.445948490915965,.108103018168070))

class Element:
    def __init__(self,points,section):
        points=np.asarray(points,float)
        if points.shape!=(6,3) or not np.isfinite(points).all():raise ValueError('finite Tri6 geometry required')
        scale=max(np.linalg.norm(points[i]-points[j]) for i,j in ((0,1),(1,2),(2,0)))
        if scale<1e-12 or any(np.linalg.norm(points[3+e]-(points[i]+points[j])/2)>64*np.finfo(float).eps*scale for e,(i,j) in enumerate(((0,1),(1,2),(2,0)))):
            raise ValueError('straight planar Tri6 required; offsets unsupported')
        ex=points[1]-points[0];ex/=np.linalg.norm(ex)
        ez=np.cross(ex,points[2]-points[0]);ez/=np.linalg.norm(ez)
        self.axes=np.array([ex,np.cross(ez,ex),ez]);self.transform=np.kron(np.eye(12),self.axes)
        self.xy=((points-points[0])@self.axes.T)[:,:2]
        self.area=abs(np.linalg.det(self.xy[1:3]))/2
        self.lagrange=np.linalg.inv(np.array([polynomial(*p) for p in self.xy]))
        self.j=self.xy[1:3].T
        e=section['young_kpa'];nu=section['poisson'];h=section['thickness_m'];g=e/(2*(1+nu))
        c=e/(1-nu*nu)*np.array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]])
        self.membrane=c*h;self.bending=c*h**3/12;self.shear=5/6*g*h;self.drilling=section['drilling_factor']*g*h
        xi,wi=np.polynomial.legendre.leggauss(5);xi=(xi+1)/2;wi=wi/2
        self.quadrature=[(np.array([(1-r)*(1-s),r,(1-r)*s]),wr*ws*(1-r)*2*self.area) for r,wr in zip(xi,wi) for s,ws in zip(xi,wi)]
        # Same Nedelec subspace, in physical coordinates, with integrated edge
        # moments instead of the Rust covariant point-evaluation functionals.
        def basis(x,y):return np.array([[1,x,y,0,0,0,-x*y,-y*y],[0,0,0,1,x,y,x*x,x*y]])
        self.basis=basis
        rows=[];targets=[]
        for i,j in ((0,1),(1,2),(2,0)):
            edge=self.xy[j]-self.xy[i]
            for power in (0,1):
                row=np.zeros(8);target=np.zeros(36)
                for z,w in zip(xi,wi):
                    x,y=self.xy[i]+z*edge
                    row+=w*z**power*(edge@basis(x,y))
                    target+=w*z**power*(edge@self.raw_shear(x,y))
                rows.append(row);targets.append(target)
        for a in range(2):
            rows.append(sum(w*basis(*(l@self.xy[:3]))[a] for l,w in self.quadrature))
            targets.append(sum(w*self.raw_shear(*(l@self.xy[:3]))[a] for l,w in self.quadrature))
        self.assumed=np.linalg.solve(np.array(rows),np.array(targets))
        local=np.zeros((36,36))
        for l,w in self.quadrature:
            bm,bb,bs,bd=self.operators(l)
            local+=w*(bm.T@self.membrane@bm+bb.T@self.bending@bb+self.shear*bs.T@bs+self.drilling*np.outer(bd,bd))
        self.stiffness=self.transform.T@local@self.transform

    def shapes(self,x,y):
        n=polynomial(x,y)@self.lagrange
        dx=np.array([0,1,0,2*x,y,0])@self.lagrange
        dy=np.array([0,0,1,0,x,2*y])@self.lagrange
        return n,dx,dy

    def raw_shear(self,x,y):
        n,dx,dy=self.shapes(x,y);bs=np.zeros((2,36))
        bs[0,2::6]=dx;bs[1,2::6]=dy;bs[0,4::6]=n;bs[1,3::6]=-n
        return bs

    def operators(self,l):
        x,y=np.asarray(l)@self.xy[:3];n,dx,dy=self.shapes(x,y)
        bm=np.zeros((3,36));bb=np.zeros((3,36));bd=np.zeros(36)
        bm[0,0::6]=dx;bm[1,1::6]=dy;bm[2,0::6]=dy;bm[2,1::6]=dx
        bb[0,4::6]=dx;bb[1,3::6]=-dy;bb[2,4::6]=dy;bb[2,3::6]=-dx
        bd[0::6]=dy/2;bd[1::6]=-dx/2;bd[5::6]=n
        return bm,bb,self.basis(x,y)@self.assumed,bd

    def load(self,traction):
        result=np.zeros(36)
        for l,w in self.quadrature:
            n,*_=self.shapes(*(l@self.xy[:3]))
            for i in range(6):result[6*i:6*i+3]+=w*n[i]*np.asarray(traction)
        return result

    def recover(self,global_relative):
        d=self.transform@global_relative;samples=[];energy=np.zeros(4)
        for l in SAMPLES:
            bm,bb,bs,bd=self.operators(l);strain=bm@d;curve=bb@d;gamma=bs@d
            samples.append(dict(barycentric=l,membrane_strain=strain,curvature_per_m=curve,
                membrane_force_kn_m=self.membrane@strain,bending_moment_kn=self.bending@curve,
                transverse_shear_strain=gamma,transverse_shear_force_kn_m=gamma*self.shear,drilling_mismatch_rad=float(bd@d)))
        for l,w in self.quadrature:
            bm,bb,bs,bd=self.operators(l);m=bm@d;b=bb@d;s=bs@d
            energy+=w/2*np.array([m@self.membrane@m,b@self.bending@b,self.shear*(s@s),self.drilling*(bd@d)**2])
        return dict(area_m2=self.area,local_axes=self.axes,samples=samples,**dict(zip(
            ('membrane_energy_kn_m','bending_energy_kn_m','transverse_shear_energy_kn_m','drilling_penalty_energy_kn_m'),energy)))
