"""Input-only hybrid RT0/P0 Richards reference on affine Tet4 cells.

Conservative head/storage Newton follows the independent-observation approach
of PR #19. Its exponential law cannot represent the application's VG law; this
extension uses direct VG/Mualem formulae and a full cell/trace nonlinear system.
No production code, condensed tangent, output, face ordering or iteration is an
input. Quadrature is assembled with SciPy Jacobi roots from the mesh order.
There is no interface adapter or strain/storage feedback here.
"""
import itertools
import numpy as np
from scipy.special import roots_jacobi


def retention(psi, law):
    psi = np.asarray(psi)
    n = law['n']; m = 1-1/n
    x = law['alpha_per_m']*np.maximum(-psi, 0.)
    se = (1+x**n)**(-m)
    theta = law['residual_water_content']+(law['saturated_water_content']-law['residual_water_content'])*se
    kr = np.sqrt(se)*(1-(1-se**(1/m))**m)**2
    storage = (law['saturated_water_content']-law['residual_water_content'])*np.expm1(-m*np.log1p(x**n))+law['saturated_specific_storage_per_m']*np.maximum(psi, 0.)
    if np.any(kr <= 0) or not np.isfinite(storage).all():
        raise ArithmeticError('invalid retention; no conductivity floor')
    return storage, kr, theta


def local_solve(mass, rhs):
    """Four-flux SPD solves in extended precision, without LAPACK downcast."""
    a = np.asarray(mass, dtype=np.longdouble).copy()
    b = np.asarray(rhs, dtype=np.longdouble).copy()
    for k in range(4):
        if np.any(a[:, k, k] <= 0): raise ArithmeticError('RT0 nonpositive pivot')
        for i in range(k+1, 4):
            factor = a[:, i, k]/a[:, k, k]
            a[:, i, k+1:] -= factor[:, None]*a[:, k, k+1:]
            b[:, i] -= factor*b[:, k]
    x = b.copy()
    for k in range(3, -1, -1):
        x[:, k] = (b[:, k]-(a[:, k, k+1:]*x[:, k+1:]).sum(axis=1))/a[:, k, k]
    return x


def rule(order):
    axes = []
    for alpha in (2, 1, 0):
        r, w = roots_jacobi(order, alpha, 0)
        axes.append(((r+1)/2, w/(2**(alpha+1))))
    p = []; weights = []
    for i, j, k in itertools.product(range(order), repeat=3):
        u, v, w = axes[0][0][i], axes[1][0][j], axes[2][0][k]
        p.append([1-u-(1-u)*v-(1-u)*(1-v)*w, u, (1-u)*v, (1-u)*(1-v)*w])
        weights.append(axes[0][1][i]*axes[1][1][j]*axes[2][1][k])
    return np.array(p), np.array(weights)


class Model:
    def __init__(self, mesh, case):
        if mesh.get('tetrahedra10') or case['schema'] != 'soil-fem-richards/1':
            raise ValueError('affine Tet4 Richards only')
        self.mesh = mesh; self.case = case
        xyz = np.asarray(mesh['nodes']); tet = np.asarray(mesh['tetrahedra'])
        vertices = xyz[tet]; self.nc = len(tet)
        det = np.linalg.det(np.transpose(vertices[:, 1:]-vertices[:, :1], (0, 2, 1)))
        if not np.all(det > 0):
            raise ValueError('positive tetrahedron orientation required')
        self.volume = det/6
        bary, w = rule(mesh.get('quadrature_order', 3))
        self.position = np.einsum('qa,eai->eqi', bary, vertices)
        self.weight = det[:, None]*w
        self.basis = (self.position[:, :, None, :]-vertices[:, None, :, :])/(3*self.volume[:, None, None, None])
        faces = {}; self.cell_keys = []
        for e, cell in enumerate(tet):
            keys = []
            for i in range(4):
                key = tuple(sorted(np.delete(cell, i))); keys.append(key)
                faces.setdefault(key, []).append((e, i))
            self.cell_keys.append(keys)
        self.keys = sorted(faces); self.owners = [faces[k] for k in self.keys]
        index = {k: i for i, k in enumerate(self.keys)}
        self.cf = np.array([[index[k] for k in row] for row in self.cell_keys])
        self.face_centres = np.array([xyz[list(k)].mean(axis=0) for k in self.keys])
        self.area = np.array([np.linalg.norm(np.cross(xyz[k[1]]-xyz[k[0]], xyz[k[2]]-xyz[k[0]]))/2 for k in self.keys])
        mats = {m['id']: m for m in case['materials']}
        self.materials = [mats[i] for i in mesh['material_ids']]
        self.inverse_k = np.array([np.linalg.inv(m['saturated_conductivity_m_s']) for m in self.materials])
        self.source = np.array([m.get('volume_source_per_s', 0.) for m in self.materials])*self.volume
        self.datum = case.get('head_datum_m', 0.)
        self.seed = case['initial_head']['head_m']+np.einsum('ei,i->e', vertices.mean(axis=1)-case['initial_head'].get('reference_point_m', [0., 0., 0.]), case['initial_head'].get('gradient', [0., 0., 0.]))-self.datum

    def local(self, h):
        values = [retention(h[e]+self.datum-self.position[e, :, 2], m['retention']) for e, m in enumerate(self.materials)]
        storage = np.array([self.weight[e]@v[0] for e, v in enumerate(values)])
        kr = np.array([v[1] for v in values])
        mass = np.einsum('eqia,eab,eqjb,eq->eij', self.basis, self.inverse_k, self.basis, self.weight/kr, optimize=True)
        return storage, mass

    def boundaries(self, period):
        head = {}; flux = {}
        sets = {name: {tuple(sorted(f[:3])) for f in faces} for name, faces in self.mesh['face_sets'].items()}
        for bc in period['boundaries']:
            c = bc['condition']
            for f, key in enumerate(self.keys):
                if key not in sets[bc['face_set']]:
                    continue
                if len(self.owners[f]) != 1 or f in head or f in flux:
                    raise ValueError('boundary overlap or interior boundary')
                if c['kind'] == 'head':
                    a = c['value']
                    head[f] = a['head_m']-self.datum+np.dot(a.get('gradient', [0., 0., 0.]), self.face_centres[f]-a.get('reference_point_m', [0., 0., 0.]))
                elif c['kind'] == 'outward_flux':
                    flux[f] = c['flux_m_s']*self.area[f]
                else:
                    raise ValueError('head/outward-flux boundaries only')
        for f, owners in enumerate(self.owners):
            if len(owners) == 1 and f not in head:
                flux.setdefault(f, 0.)
        return head, flux

    def step(self, h, trace, old_storage, period, dt):
        heads, prescribed_flux = self.boundaries(period)
        free = np.array([f for f in range(len(self.keys)) if f not in heads], dtype=int)
        trace = trace.copy()
        for f, value in heads.items(): trace[f] = value
        x = np.concatenate((h, trace[free]))

        def evaluate(x):
            hh = x[:self.nc]; ll = trace.copy(); ll[free] = x[self.nc:]
            storage, mass = self.local(hh)
            q = local_solve(mass, hh[:, None]-ll[self.cf])
            rate = (storage-old_storage)/dt
            cell = q.sum(axis=1)+rate-self.source
            scale = np.abs(q).sum(axis=1)+np.abs(rate)+np.abs(self.source)
            face = np.array([sum(q[e, j] for e, j in self.owners[f])-prescribed_flux.get(f, 0.) for f in free])
            fscale = np.array([sum(abs(q[e, j]) for e, j in self.owners[f])+abs(prescribed_flux.get(f, 0.)) for f in free])
            r = np.concatenate((cell, face)); scales = np.concatenate((scale, fscale))
            return r, scales, hh, ll, storage, q

        max_linear = 0.
        for iteration in range(60):
            r, scale, hh, ll, storage, q = evaluate(x)
            merit = float(np.max(np.abs(r)/(1e-16+1e-8*scale)))
            if merit < .01:
                break
            jac = np.empty((len(x), len(x)))
            for col in range(len(x)):
                delta = 1e-6*(1+abs(x[col])); xp = x.copy(); xm = x.copy()
                xp[col] += delta; xm[col] -= delta
                jac[:, col] = (evaluate(xp)[0]-evaluate(xm)[0])/(2*delta)
            dx = np.linalg.solve(jac, np.asarray(-r, dtype=float))
            lin = float(np.linalg.norm(jac@dx+r)/np.linalg.norm(r))
            if lin > 1e-10:
                raise ArithmeticError(('true hydraulic linear residual', lin))
            max_linear = max(max_linear, lin)
            for cut in range(20):
                candidate = x+dx*2.**(-cut)
                if np.linalg.norm(evaluate(candidate)[0]) < np.linalg.norm(r):
                    x = candidate; break
            else:
                raise ArithmeticError(('hydraulic line search', merit))
        else:
            raise ArithmeticError('hydraulic Newton budget')
        out = sum(q[e, j] for f, owners in enumerate(self.owners) if len(owners) == 1 for e, j in owners)
        return dict(cell_head_m=hh+self.datum, face_head_m=ll+self.datum, outflow_m3_s=q, storage_offset_m3=storage,
                    boundary_outflow_m3_s=float(out), source_m3_s=float(self.source.sum()),
                    storage_rate_m3_s=float((storage-old_storage).sum()/dt),
                    maximum_scaled_residual=merit, maximum_absolute_cell_residual_m3_s=float(np.max(np.abs(r[:self.nc]))),
                    maximum_true_linear_residual=max_linear, iterations=iteration)

    def solve(self):
        h = self.seed.astype(np.longdouble); storage = self.local(h)[0]
        a = self.case['initial_head']
        trace = (a['head_m']-self.datum+(self.face_centres-a.get('reference_point_m', [0., 0., 0.]))@np.array(a.get('gradient', [0., 0., 0.]))).astype(np.longdouble)
        t = 0.; cumulative = 0.; initial = storage.sum(); rows = []
        max_dt = self.case['settings']['maximum_time_step_s']
        for period in self.case['analysis']['periods']:
            end = period['end_time_s']
            while t < end:
                dt = min(max_dt, end-t)
                f = self.step(h, trace, storage, period, dt)
                t += dt
                cumulative += dt*(f['boundary_outflow_m3_s']-f['source_m3_s'])
                f.update(time_s=t, dt_s=dt, period=period['name'], cumulative_net_outflow_m3=cumulative,
                         global_mass_defect_m3=float(f['storage_offset_m3'].sum()-initial+cumulative))
                rows.append(f)
                h = f['cell_head_m']-self.datum; trace = f['face_head_m']-self.datum; storage = f['storage_offset_m3']
        return dict(complete=True, reference_uses_native_results=False, face_keys=self.keys, rows=rows,
                    initial_storage_offset_m3=float(initial), input_cells=self.nc,
                    primary_head_unknown='cell P0 and hybrid face trace; no nodal head unknown')
