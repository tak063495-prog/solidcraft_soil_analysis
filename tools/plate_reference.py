"""Independent flat CST/DKT: solve 12 polynomial rotation constraints.

No Rust shape functions or solved fields are used. Quadratic beta is found
from vertex rotations, edge-normal interpolation and integrated Kirchhoff
edge constraints; a separate 5x5 Gauss/Duffy rule integrates its energy.
"""
import numpy as np


def polynomial(x, y):
    return np.array([1., x, y, x*x, x*y, y*y])


class Element:
    def __init__(self, points, section):
        points=np.asarray(points,float)
        ex=points[1]-points[0];ex/=np.linalg.norm(ex)
        ez=np.cross(ex,points[2]-points[0]);ez/=np.linalg.norm(ez)
        self.axes=np.array([ex,np.cross(ez,ex),ez])
        self.transform=np.kron(np.eye(6),self.axes)
        self.xy=((points-points[0])@self.axes.T)[:,:2]
        self.area=abs(np.linalg.det(np.column_stack((self.xy[1]-self.xy[0],self.xy[2]-self.xy[0]))))/2
        self.affine=np.linalg.inv(np.column_stack((np.ones(3),self.xy)))
        self.gradient=self.affine[1:].T
        equations=[];targets=[]
        for i,(x,y) in enumerate(self.xy):
            for component in range(2):
                row=np.zeros(12);row[6*component:6*component+6]=polynomial(x,y)
                target=np.zeros(18);target[6*i+(4 if component==0 else 3)]=1 if component==0 else -1
                equations.append(row);targets.append(target)
        xi,wi=np.polynomial.legendre.leggauss(5)
        for i,j in ((0,1),(1,2),(2,0)):
            edge=self.xy[j]-self.xy[i];length=np.linalg.norm(edge);tangent=edge/length;normal=np.array([-tangent[1],tangent[0]])
            mid=(self.xy[i]+self.xy[j])/2
            equations.append(np.kron(normal,polynomial(*mid)))
            target=np.zeros(18)
            for n in (i,j):target[6*n+4]=normal[0]/2;target[6*n+3]=-normal[1]/2
            targets.append(target)
            average=sum((w/2*polynomial(*(self.xy[i]+(z+1)/2*edge)) for z,w in zip(xi,wi)),np.zeros(6))
            equations.append(np.kron(tangent,average))
            target=np.zeros(18);target[6*i+2]=1/length;target[6*j+2]=-1/length;targets.append(target)
        self.coefficients=np.linalg.solve(np.array(equations),np.array(targets))
        e=section['young_kpa'];nu=section['poisson'];h=section['thickness_m']
        c=e/(1-nu*nu)*np.array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]])
        self.membrane=c*h;self.bending=c*h**3/12;self.drilling=section['drilling_factor']*e/(2*(1+nu))*h
        self.bm=np.zeros((3,18))
        for i,(dx,dy) in enumerate(self.gradient):
            self.bm[0,6*i]=dx;self.bm[1,6*i+1]=dy;self.bm[2,6*i:6*i+2]=[dy,dx]
        self.quadrature=[]
        local=np.zeros((18,18))
        for r,wr in zip((xi+1)/2,wi/2):
            for s,ws in zip((xi+1)/2,wi/2):
                bary=np.array([(1-r)*(1-s),r,(1-r)*s]);weight=wr*ws*(1-r)*self.area*2
                bb,bd=self.operators(bary);self.quadrature.append((bary,weight,bb,bd))
                local+=weight*(self.bm.T@self.membrane@self.bm+bb.T@self.bending@bb+self.drilling*np.outer(bd,bd))
        self.stiffness=self.transform.T@local@self.transform

    def operators(self,bary):
        x,y=np.asarray(bary)@self.xy
        dx=np.array([0.,1.,0.,2*x,y,0.]);dy=np.array([0.,0.,1.,0.,x,2*y])
        derivatives=np.zeros((3,12));derivatives[0,:6]=dx;derivatives[1,6:]=dy
        derivatives[2,:6]=dy;derivatives[2,6:]=dx
        bb=derivatives@self.coefficients;bd=np.zeros(18)
        for i,(gx,gy) in enumerate(self.gradient):
            bd[6*i:6*i+2]=[gy/2,-gx/2];bd[6*i+5]=bary[i]
        return bb,bd

    def load(self,traction):
        out=np.zeros(18)
        for bary,weight,_,_ in self.quadrature:
            for i,shape in enumerate(bary):out[6*i:6*i+3]+=weight*shape*np.asarray(traction)
        return out

    def recover(self,global_relative):
        d=self.transform@global_relative;strain=self.bm@d;force=self.membrane@strain
        samples=[]
        for bary in ((2/3,1/6,1/6),(1/6,2/3,1/6),(1/6,1/6,2/3)):
            bb,bd=self.operators(bary);curvature=bb@d
            samples.append(dict(barycentric=bary,membrane_strain=strain,curvature_per_m=curvature,
                membrane_force_kn_m=force,bending_moment_kn=self.bending@curvature,drilling_mismatch_rad=float(bd@d)))
        bending_energy=sum(float((bb@d)@self.bending@(bb@d))*w/2 for _,w,bb,_ in self.quadrature)
        drilling_energy=sum(float(bd@d)**2*self.drilling*w/2 for _,w,_,bd in self.quadrature)
        return dict(area_m2=self.area,local_axes=self.axes,samples=samples,
            membrane_energy_kn_m=float(strain@force)*self.area/2,bending_energy_kn_m=bending_energy,
            drilling_penalty_energy_kn_m=drilling_energy)
