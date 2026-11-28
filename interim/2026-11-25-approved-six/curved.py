"""Initially circular planar Kirchhoff elements and directional tangent derivatives."""
import numpy as np
from arch_elements import Arch
from continuation import Evaluation
from references import A, EI, R, span, theta0


class CurvedArch(Arch):
    """Quadratic rotation with a condensed bubble and constant section force vector."""

    def __init__(self, n):
        super().__init__(n)
        x, w = np.polynomial.legendre.leggauss(8)
        s = (1 + x) / 2
        length = span / n
        self.weights = w * length / 2
        self.shape = np.stack([1 - s, s, 4 * s * (1 - s)], axis=1)
        gradients = np.stack([-np.ones_like(s), np.ones_like(s), 4 - 8*s], axis=1) / length
        self.bending = gradients.T @ (self.weights[:, None] * gradients)
        self.angle = theta0 + np.pi/2 + (np.arange(n)[:, None] + s) * length
        self.axial_compliance = EI / (1e6 * A * R**2)

    def __call__(self, parent, z):
        u = np.zeros_like(self.load)
        u[self.free] = z[:-1]
        e = u[self.ids]
        chord = self.d0 + e @ self.L.T
        # Stationarity in the bubble rotation and the two end-force components
        # enforces bending equilibrium and integrated chord compatibility.
        rotations = np.column_stack((e[:, 2], e[:, 5], np.zeros(len(e))))
        force = np.zeros_like(chord)
        c = self.axial_compliance
        w = self.weights
        b = self.shape
        for _ in range(35):
            angle = self.angle + rotations @ b.T
            tangent = np.stack((np.cos(angle), np.sin(angle)), axis=2)
            normal = np.stack((-tangent[:,:,1], tangent[:,:,0]), axis=2)
            axial = np.einsum("eqi,ei->eq", tangent, force)
            shear = np.einsum("eqi,ei->eq", normal, force)
            stretch = 1 + axial*c
            moment = rotations @ self.bending - (w * stretch * shear) @ b
            closure = chord - np.einsum("eq,eqi->ei", w * stretch, tangent)
            hrr = self.bending + np.einsum("eq,qi,qj->eij", w*(stretch*axial-c*shear**2), b, b)
            hrp = -np.einsum("q,qi,eqj->eij", w,b,stretch[:,:,None]*normal + (c*shear)[:,:,None]*tangent)
            hpp = -c*np.einsum("q,eqi,eqj->eij", w,tangent,tangent)
            local = np.concatenate((np.concatenate((hrr[:,2:,2:],hrp[:,2:]),axis=2),
                                    np.concatenate((hrp[:,2:].transpose(0,2,1),hpp),axis=2)),axis=1)
            error = np.column_stack((moment[:,2],closure))
            if np.max(np.abs(error)) < 5e-14:
                break
            step = np.linalg.solve(local, error[:,:,None])[:,:,0]
            rotations[:,2] -= step[:,0]
            force -= step[:,1:]
        else:
            raise ArithmeticError("curved beam element did not equilibrate")

        # Hessian of the stationary element potential condensed to the six
        # nodal coordinates. Axial extension is eliminated by its linear law.
        coupling = np.zeros((len(e),6,3))
        coupling[:,:,1:] = self.L.T
        coupling[:,[2,5],0] = hrr[:,:2,2]
        coupling[:,[2,5],1:] = hrp[:,:2,:]
        ke = -coupling @ np.linalg.solve(local, coupling.transpose(0,2,1))
        for i, ir in enumerate((2,5)):
            for j, jr in enumerate((2,5)):
                ke[:,ir,jr] += hrr[:,i,j]
        fe = force @ self.L
        fe[:,[2,5]] += moment[:,:2]
        internal = np.zeros_like(u)
        stiff = np.zeros((len(u),len(u)))
        np.add.at(internal,self.ids,fe)
        np.add.at(stiff,(self.ids[:,:,None],self.ids[:,None,:]),ke)
        q = -self.load[self.free]
        return Evaluation(internal[self.free]+q*z[-1],stiff[np.ix_(self.free,self.free)],q,{},
                          max(1.,abs(z[-1])),details=dict(crown=-u[3*self.apex+1]*R))


class AnalyticChordArch(Arch):
    """Differentiate the original chord-element force twice along a displacement."""

    def tangent_action_derivative(self, z, v):
        u = np.zeros_like(self.load)
        direction = np.zeros_like(u)
        u[self.free], direction[self.free] = z[:-1], v
        e, de = u[self.ids], direction[self.ids]
        d = self.d0 + e @ self.L.T
        length = np.linalg.norm(d, axis=1)
        t = d / length[:,None]
        n = np.column_stack((-t[:,1],t[:,0]))
        g, h = t @ self.L, n @ self.L
        dl = np.einsum("ei,ei->e",g,de)
        da = np.einsum("ei,ei->e",h,de) / length
        ddl = length*da**2
        dda = -2*da*dl/length
        rot = np.arctan2(self.d0[:,0]*d[:,1]-self.d0[:,1]*d[:,0], (self.d0*d).sum(axis=1))
        basic = np.column_stack((length-self.h,e[:,2]-rot,e[:,5]-rot))
        dbasic = np.column_stack((dl,de[:,2]-da,de[:,5]-da))
        ddbasic = np.column_stack((ddl,-dda,-dda))
        f, df, ddf = (np.einsum("eij,ej->ei",self.D,s) for s in (basic,dbasic,ddbasic))
        # Internal force = axial*g - (m1+m2)*h/l + nodal moments.
        dg, dh = da[:,None]*h, -da[:,None]*g
        ddg = dda[:,None]*h-da[:,None]**2*g
        ddh = -dda[:,None]*g-da[:,None]**2*h
        m = f[:,1:].sum(axis=1)/length
        dm = (df[:,1:].sum(axis=1)-m*dl)/length
        ddm = (ddf[:,1:].sum(axis=1)-2*dm*dl-m*ddl)/length
        fe = (ddf[:,:1]*g + 2*df[:,:1]*dg + f[:,:1]*ddg
              - ddm[:,None]*h - 2*dm[:,None]*dh - m[:,None]*ddh)
        fe[:,[2,5]] += ddf[:,1:]
        result = np.zeros_like(u)
        np.add.at(result,self.ids,fe)
        return result[self.free]
