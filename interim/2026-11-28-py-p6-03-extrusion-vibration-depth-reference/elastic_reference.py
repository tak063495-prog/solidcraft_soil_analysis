"""Independent small-strain Q4/Hex8 elasticity; SI, tension-positive stresses.

Only NumPy/SciPy, no production imports. Full 2^d quadrature, engineering shear.
The Lamé, foundation and all constitutive/assembly paths here are local to this file.
"""
from itertools import product
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


def elastic_matrix(E, nu):
    if not (E > 0 and -1 < nu < .5):
        raise ValueError("invalid isotropic material")
    mu = E / (2 * (1 + nu))
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    D[np.arange(3), np.arange(3)] += 2 * mu
    D[3:, 3:] = np.eye(3) * mu
    return D


def signs(dim):
    q = np.array([[-1,-1], [1,-1], [1,1], [-1,1]], dtype=float)
    return q if dim == 2 else np.vstack([np.c_[q, -np.ones(4)], np.c_[q, np.ones(4)]])


def element_points(x):
    dim = x.shape[1]
    sg = signs(dim)
    for xi in product([-1 / np.sqrt(3), 1 / np.sqrt(3)], repeat=dim):
        factors = (1 + sg * xi) / 2
        N = np.prod(factors, axis=1)
        dN = np.column_stack([sg[:, a] / 2 * np.prod(np.delete(factors, a, axis=1), axis=1)
                              for a in range(dim)])
        J = x.T @ dN
        volume = np.linalg.det(J)
        if volume <= 0:
            raise ValueError("nonpositive element Jacobian")
        grad = dN @ np.linalg.inv(J)
        B = np.zeros((6, len(x) * dim))
        for i, g in enumerate(grad):
            B[0, dim*i] = g[0]
            B[1, dim*i+1] = g[1]
            B[3, dim*i] = g[1]
            B[3, dim*i+1] = g[0]
            if dim == 3:
                B[2, dim*i+2] = g[2]
                B[4, dim*i+1] = g[2]
                B[4, dim*i+2] = g[1]
                B[5, dim*i] = g[2]
                B[5, dim*i+2] = g[0]
        yield N, B, float(volume)


def solve(nodes, cells, force, prescribed, E=30e6, nu=.25, ezz=0.):
    dim = nodes.shape[1]
    D = elastic_matrix(E, nu)
    rows, cols, values = [], [], []
    initial = np.zeros(6)
    if dim == 2:
        initial[2] = ezz
    rhs = np.asarray(force, float).reshape(-1).copy()
    records = []
    for cell in cells:
        dof = (cell[:, None] * dim + np.arange(dim)).ravel()
        pts = list(element_points(nodes[cell]))
        ke = sum(B.T @ D @ B * w for _, B, w in pts)
        rhs[dof] -= sum(B.T @ D @ initial * w for _, B, w in pts)
        rows.extend(np.repeat(dof, len(dof)))
        cols.extend(np.tile(dof, len(dof)))
        values.extend(ke.ravel())
        records.append((cell, dof, pts))
    K = coo_matrix((values, (rows, cols)), shape=(len(rhs), len(rhs))).tocsr()
    fixed = np.array(sorted(prescribed), dtype=int)
    free = np.setdiff1d(np.arange(len(rhs)), fixed)
    u = np.zeros_like(rhs)
    u[fixed] = [prescribed[i] for i in fixed]
    u[free] = spsolve(K[free][:, free], rhs[free] - K[free][:, fixed] @ u[fixed])
    if not np.isfinite(u).all():
        raise ArithmeticError("nonfinite elasticity solution")
    residual = K @ u - rhs
    gp = []
    for eid, (cell, dof, pts) in enumerate(records):
        for gid, (N, B, w) in enumerate(pts):
            strain = B @ u[dof] + initial
            gp.append({"element": eid, "gp": gid, "x": N @ nodes[cell],
                       "u": N @ u.reshape(-1, dim)[cell], "stress": D @ strain,
                       "strain": strain, "weight": w})
    applied = np.asarray(force).reshape(-1, dim).sum(axis=0)
    reactions = np.zeros_like(u)
    reactions[fixed] = residual[fixed]
    boundary_half_work = .5 * u @ (np.asarray(force).ravel()+reactions)
    if dim == 2:
        boundary_half_work += .5 * ezz * sum(g["stress"][2]*g["weight"] for g in gp)
    return {"u": u.reshape(-1, dim), "gp": gp,
            "free_residual_relative": float(np.linalg.norm(residual[free]) / max(np.linalg.norm(rhs), 1.)),
            "force_balance_N": (reactions.reshape(-1, dim).sum(axis=0) + applied).tolist(),
            "strain_energy_J": float(sum(.5 * g["strain"] @ g["stress"] * g["weight"] for g in gp)),
            "external_work_half_J": float(.5 * u @ np.asarray(force).ravel()),
            "boundary_work_half_J": float(boundary_half_work),
            "nodes": len(nodes), "elements": len(cells)}


def lame(r, E=30e6, nu=.25, a=1., b=2., pi=1e5, po=0., ezz=0.):
    """Radial equilibrium + 3D Hooke law, constant imposed axial strain.

    A=(pi*a²-po*b²)/(b²-a²), B=(pi-po)*a²*b²/(b²-a²).
    No end-cap force condition is implied by imposed ezz.
    """
    A = (pi*a*a - po*b*b) / (b*b-a*a)
    B = (pi-po)*a*a*b*b / (b*b-a*a)
    ur = ((1+nu)*(1-2*nu)*A/E - nu*ezz)*r + (1+nu)*B/(E*r)
    return ur, A-B/r**2, A+B/r**2, 2*nu*A+E*ezz


def ring(nr=2, nt=4, dim=2, ezz=0.):
    """Quarter annulus a=1,b=2, symmetry sides; extruded length=1 in 3D.

    Straight polygonal Q4/Hex8 geometry, 2x2(/x2) Gauss. Pressure acts normal
    to each polygonal face. All 3D uz=ezz*z, equivalent to imposed GPS.
    """
    layers = [0.] if dim == 2 else [0., 1.]
    nodes = np.array([[r*np.cos(t), r*np.sin(t), z][:dim]
                      for z in layers for r in np.linspace(1, 2, nr+1)
                      for t in np.linspace(0, np.pi/2, nt+1)])
    nn = (nr+1)*(nt+1)
    quads = np.array([[i*(nt+1)+j, (i+1)*(nt+1)+j,
                       (i+1)*(nt+1)+j+1, i*(nt+1)+j+1]
                      for i in range(nr) for j in range(nt)])
    cells = quads if dim == 2 else np.c_[quads, quads+nn]
    f = np.zeros_like(nodes)
    for j in range(nt):
        edge = nodes[j+1, :2] - nodes[j, :2]
        load = 1e5 * np.array([edge[1], -edge[0]])
        if dim == 2:
            f[[j, j+1]] += load/2
        else:
            f[[j, j+1, j+nn, j+1+nn], :2] += load/4
    fixed = {}
    for i, x in enumerate(nodes):
        if abs(x[0]) < 1e-12:
            fixed[dim*i] = 0.
        if abs(x[1]) < 1e-12:
            fixed[dim*i+1] = 0.
        if dim == 3:
            fixed[dim*i+2] = ezz*x[2]
    out = solve(nodes, cells, f, fixed, ezz=ezz)
    eu, es, normu, norms, ez, normz = 0., 0., 0., 0., 0., 0.
    relation = []
    for g in out["gp"]:
        x = g["x"][:2]
        r = np.linalg.norm(x)
        ur, sr, st, sz = lame(r, ezz=ezz)
        c, s = x/r
        exact = np.array([sr*c*c+st*s*s, sr*s*s+st*c*c, sz, (sr-st)*c*s, 0., 0.])
        ue = ur*x/r
        w = g["weight"]
        eu += np.sum((g["u"][:2]-ue)**2)*w
        es += np.sum((g["stress"]-exact)**2*np.array([1,1,1,2,2,2]))*w
        normu += np.sum(ue**2)*w
        norms += np.sum(exact**2*np.array([1,1,1,2,2,2]))*w
        ez += (g["stress"][2]-sz)**2*w
        normz += sz**2*w
        sig = g["stress"]
        relation.append(abs(sig[2] - .25*(sig[0]+sig[1]) - 30e6*ezz))
    return {**{k:v for k,v in out.items() if k not in ("u", "gp")},
            "dim": dim, "nr": nr, "nt": nt, "ezz": ezz,
            "displacement_relative_L2": float(np.sqrt(eu/normu)),
            "stress_relative_L2": float(np.sqrt(es/norms)),
            "sigma_zz_relative_L2": float(np.sqrt(ez/normz)),
            "sigma_zz_relation_max_Pa": float(max(relation)),
            "sigma_zz_exact_Pa": float(lame(1., ezz=ezz)[3])}


def foundation(B=1., L=2., level=1, dim=3):
    """Finite elastic box, x in [-2,2], y in [-2,0], z in [-L/2,L/2].

    Top pressure 1e5 Pa over |x|<=B/2, all z. Bottom fixed, x-side ux=0,
    z-ends traction free. Plane strain 2D has identical cross-section/support.
    Cross-section boundaries are finite; not a halfspace solution.
    """
    nx, ny = 8*level, 4*level
    nz = int(round(2*L*level)) if dim == 3 else 0
    zs = np.linspace(-L/2, L/2, nz+1) if dim == 3 else [0.]
    nodes = np.array([[x,y,z][:dim] for z in zs for y in np.linspace(-2,0,ny+1)
                      for x in np.linspace(-2,2,nx+1)])
    nn = (nx+1)*(ny+1)
    quads = np.array([[j*(nx+1)+i, j*(nx+1)+i+1,
                       (j+1)*(nx+1)+i+1, (j+1)*(nx+1)+i]
                      for j in range(ny) for i in range(nx)])
    cells = quads if dim == 2 else np.vstack([np.c_[quads+k*nn, quads+(k+1)*nn] for k in range(nz)])
    f = np.zeros_like(nodes)
    dx = 4/nx
    for i in range(nx):
        xm = -2+(i+.5)*dx
        if abs(xm) < B/2:
            ids = [ny*(nx+1)+i, ny*(nx+1)+i+1]
            if dim == 2:
                f[ids,1] -= 1e5*dx/2
            else:
                for k in range(nz):
                    face = np.array(ids+list(np.array(ids)+nn))+k*nn
                    f[face,1] -= 1e5*dx*(L/nz)/4
    fixed = {}
    for i,x in enumerate(nodes):
        if abs(x[1]+2) < 1e-12:
            for a in range(dim):
                fixed[dim*i+a] = 0.
        if abs(abs(x[0])-2) < 1e-12:
            fixed[dim*i] = 0.
    out = solve(nodes, cells, f, fixed)
    # Evaluate the FE basis at fixed physical points (0,-1,z), one-sided
    # +x,+y trace at mesh edges. No stress smoothing or cross-mesh interpolation.
    profile = []
    stations = np.linspace(0,L/2,9) if dim == 3 else [0.]
    for z in stations:
        target = np.array([0.,-1.,z][:dim])
        matching = []
        for cell in cells:
            xyz = nodes[cell]
            low, high = xyz.min(axis=0), xyz.max(axis=0)
            if np.all(target >= low-1e-12) and np.all(target <= high+1e-12):
                matching.append((tuple(low[::-1]), cell))
        cell = max(matching, key=lambda p:p[0])[1]
        xyz = nodes[cell]
        xi = 2*(target-xyz.min(axis=0))/(xyz.max(axis=0)-xyz.min(axis=0))-1
        sg = signs(dim)
        factors = (1+sg*xi)/2
        grad_nat = np.column_stack([sg[:,a]/2*np.prod(np.delete(factors,a,axis=1),axis=1)
                                    for a in range(dim)])
        grad = grad_nat @ np.linalg.inv(xyz.T @ grad_nat)
        displacement = out["u"][cell]
        du = displacement.T @ grad
        strain = np.array([du[0,0], du[1,1], du[2,2] if dim == 3 else 0.,
                           du[0,1]+du[1,0], du[1,2]+du[2,1] if dim == 3 else 0.,
                           du[0,2]+du[2,0] if dim == 3 else 0.])
        stress = elastic_matrix(30e6,.25) @ strain
        profile.append({"x_m":float(target[0]), "y_m":float(target[1]), "z_m":float(z),
                        "sigma_zz_Pa":float(stress[2]), "sigma_yy_Pa":float(stress[1]),
                        "plane_strain_relation_defect_Pa":float(stress[2]-.25*(stress[0]+stress[1]))})
    return {**{k:v for k,v in out.items() if k not in ("u", "gp")},
            "B_m": B, "L_m": L, "level": level, "dim": dim, "profile":profile}
