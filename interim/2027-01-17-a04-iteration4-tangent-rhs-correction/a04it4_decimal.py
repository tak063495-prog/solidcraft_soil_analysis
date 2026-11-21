"""Affine Tet10 elastic reference and pivoted Decimal solves in fixed DOF order."""
from decimal import Decimal
import numpy as np

D = Decimal
EDGES = ((0,1),(1,2),(2,0),(0,3),(1,3),(2,3))


def decimal_array(values):
    return np.vectorize(lambda x: D.from_float(float(x)), otypes=[object])(values)


def zeros(shape):
    return np.full(shape, D(0), dtype=object)


def cross(a,b):
    return np.array([a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]],dtype=object)


def strain_matrix(g):
    b = zeros((6,30))
    for i,(x,y,z) in enumerate(g):
        j = 3*i
        b[0,j], b[1,j+1], b[2,j+2] = x,y,z
        b[3,j], b[3,j+1] = y,x
        b[4,j+1], b[4,j+2] = z,y
        b[5,j], b[5,j+2] = z,x
    return b


def elastic(mesh, ground):
    xyz = decimal_array(mesh['nodes'])
    k, seed = zeros((xyz.size,xyz.size)), zeros(xyz.size)
    materials = {m['id']:m['mohr_coulomb'] for m in ground['materials']}
    bodies = {r['id']:r['body'] for r in mesh['regions']}
    states = ground['initialization']['geostatic']['seed']['states']
    for nodes,region,mid in zip(mesh['tetrahedra10'],mesh['region_ids'],mesh['material_ids'],strict=True):
        p = xyz[nodes]
        for i,(a,b) in enumerate(EDGES,4):
            assert np.array_equal(2*p[i],p[a]+p[b])
        a,b,c = p[1:4]-p[0]
        det = a @ cross(b,c)
        assert det > 0
        grads = np.array([cross(b,c),cross(c,a),cross(a,b)])/det
        dl = np.vstack((-grads.sum(axis=0),grads))
        young,nu = decimal_array([materials[mid]['young_kpa'],materials[mid]['poisson']])
        mu = young/(2*(1+nu))
        lam = young*nu/((1+nu)*(1-2*nu))
        tensor = zeros((6,6))
        tensor[:3,:3] = lam
        for j in range(3):
            tensor[j,j] += 2*mu
            tensor[j+3,j+3] = mu
        vertices = []
        for vertex in range(4):
            l = [int(j==vertex) for j in range(4)]
            grad = [(4*l[i]-1)*dl[i] for i in range(4)]
            grad.extend(4*(l[i]*dl[j]+l[j]*dl[i]) for i,j in EDGES)
            vertices.append(strain_matrix(grad))
        mean = sum(vertices)/4
        # E[L_i L_j]=(1+delta_ij)/20. B is affine in barycentric coordinates.
        local = 16*(mean.T @ tensor @ mean)
        for v in vertices:
            local += v.T @ tensor @ v
        ids = [3*n+i for n in nodes for i in range(3)]
        k[np.ix_(ids,ids)] += local*(det/120)
        stress = decimal_array(states[bodies[region]]['effective_stress_compression_positive_kpa'])
        seed[ids] -= (mean.T @ stress)*(det/6)
    return k,seed


def solve(a,b):
    """Gaussian elimination with deterministic partial pivoting, including multiple RHS."""
    m = a.copy()
    rhs = b.reshape(len(a),-1).copy()
    for j in range(len(a)):
        p = j+np.argmax(abs(m[j:,j]))
        if m[p,j] == 0:
            raise ArithmeticError('singular reference matrix')
        if p != j:
            m[[j,p]], rhs[[j,p]] = m[[p,j]],rhs[[p,j]]
        factors = m[j+1:,j]/m[j,j]
        m[j+1:,j+1:] -= factors[:,None]*m[j,j+1:]
        rhs[j+1:] -= factors[:,None]*rhs[j]
    for j in reversed(range(len(a))):
        rhs[j] = (rhs[j]-m[j,j+1:]@rhs[j+1:])/m[j,j]
    return rhs.reshape(b.shape)
