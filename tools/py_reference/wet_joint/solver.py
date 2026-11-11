"""Independent small-strain 3D elastic bulk + unilateral Coulomb reference.

Inputs alone construct operators, forces, and a uniform load schedule. Native
fields, tangents, converged states and accepted increments are never inputs.
The bulk MC elastic-domain condition is checked at every committed point;
plastic bulk, transient flow, structures and curved cells are rejected.
"""
import numpy as np
from scipy.linalg import null_space, solve
from .geometry import (basis, tet_geometry, hydro_volume, triangle_rule, topology,
                       face_geometry, pressure_face)


def elastic(young, nu):
    shear = young/(2*(1+nu)); lame = young*nu/((1+nu)*(1-2*nu))
    d = np.zeros((6, 6)); d[:3, :3] = lame
    d[np.arange(3), np.arange(3)] += 2*shear
    d[np.arange(3, 6), np.arange(3, 6)] = shear
    return d


def linear_correction(matrix, rhs):
    """LU with the true original-system relative residual, including tiny RHS."""
    correction = solve(matrix, rhs, assume_a='gen', check_finite=True)
    relative = np.linalg.norm(matrix@correction-rhs)/max(np.linalg.norm(rhs), 1e-30)
    if not np.isfinite(relative) or relative > 1e-10:
        raise RuntimeError(('true linear residual', float(relative)))
    return correction


def coulomb(start_jump, start_traction, target, normal, kn, kt, mu, partitions):
    """Integrate a straight jump path by projection onto the Coulomb disk.

    Includes closed fractions of opening/reclosing segments. Returns net
    plastic vector and accumulated *length* separately. Trials never mutate
    committed history. A candidate axis permits batched tangent evaluations.
    """
    target = np.asarray(target)
    if target.ndim == 2:
        target = target[:, None, :]
    n = normal[:, None, :]; kn = kn[:, None]; kt = kt[:, None]; mu = mu[:, None]
    old = start_jump[:, None, :]
    delta = (target-old)/partitions
    dn = (delta*n).sum(axis=-1)
    tangential = delta-dn[..., None]*n
    traction = start_traction[:, None, :]
    tau = np.broadcast_to(traction-(traction*n).sum(axis=-1)[..., None]*n, target.shape).copy()
    plastic = np.zeros_like(target); accumulated = np.zeros(target.shape[:-1])
    for k in range(partitions):
        g0 = ((old+k*delta)*n).sum(axis=-1)
        g1 = ((old+(k+1)*delta)*n).sum(axis=-1)
        fraction = np.where((g0 <= 0) & (g1 <= 0), 1., 0.)
        opening = (g0 <= 0) & (g1 > 0)
        closing = (g0 > 0) & (g1 <= 0)
        np.divide(-g0, dn, out=fraction, where=opening)
        np.divide(g1, dn, out=fraction, where=closing)
        trial = np.where((g0 > 0)[..., None], 0., tau)+kt[..., None]*tangential*fraction[..., None]
        radius = mu*kn*np.maximum(-g1, 0.)
        length = np.linalg.norm(trial, axis=-1)
        ratio = np.minimum(1., np.divide(radius, length, out=np.ones_like(length), where=length > 0))
        projected = trial*ratio[..., None]
        flow = (trial-projected)/kt[..., None]
        plastic += flow; accumulated += np.linalg.norm(flow, axis=-1)
        tau = np.where((g1 <= 0)[..., None], projected, 0.)
    gap = (target*n).sum(axis=-1)
    return tau+kn[..., None]*np.minimum(gap, 0.)[..., None]*n, plastic, accumulated


class Model:
    def __init__(self, mesh, ground, case, *, pc=False):
        self.mesh = mesh; self.ground = ground; self.case = case
        self.node_body = [set() for _ in mesh['nodes']]
        for ei, cell in enumerate(mesh['tetrahedra10']):
            for node in cell:
                self.node_body[node].add(mesh['region_ids'][ei])
        self.node_body = [','.join(map(str, sorted(ids))) for ids in self.node_body]
        self.xyz = np.asarray(mesh['nodes'], dtype=float); nd = self.xyz.size
        if nd > (15000 if pc else 900) or len(mesh['tetrahedra']) > (768 if pc else 96):
            raise ValueError('specimen budget exceeded; larger studies require --user-pc')
        cfg = ground['initialization']; geo = cfg['geostatic']
        if (cfg['method'] != 'equilibrate' or geo['initial_boundary_factor'] != 0
                or geo['seed']['kind'] != 'body_states' or case['settings']['kinematics'] != 'explicit'
                or any(s['kind'] != 'boundary' for s in case['stages'])
                or any(case.get(k) for k in ('beam_system', 'anchor_system', 'initially_inactive_bodies'))
                or case['joint_system'].get('structure_bodies')):
            raise ValueError('unsupported physical scope')
        self.level = float(cfg['water_level_m']); gamma = float(cfg['water_unit_weight_kn_m3'])
        self.k = np.zeros((nd, nd)); self.seed = np.zeros(nd); self.base = np.zeros(nd)
        materials = {m['id']: m for m in ground['materials']}
        regions = {r['id']: r['body'] for r in mesh['regions']}
        self.volumes = []
        for ei, cell in enumerate(mesh['tetrahedra10']):
            material = materials[mesh['material_ids'][ei]]
            if 'mcc' in material or 'invariant_double_hardening' in material:
                raise ValueError('plastic bulk is outside reference scope')
            mat = material['mohr_coulomb']; d = elastic(mat['young_kpa'], mat['poisson'])
            seed = geo['seed']['states'][regions[mesh['region_ids'][ei]]]
            if (seed['equivalent_plastic_strain'] != 0 or seed['plastic_work_kpa'] != 0
                    or np.any(seed['plastic_strain_compression_positive'])):
                raise ValueError('nonelastic bulk seed')
            s0 = np.asarray(seed['effective_stress_compression_positive_kpa'])
            ids = (3*np.array(cell)[:, None]+np.arange(3)).ravel()
            shape, b, positions, volume = tet_geometry(self.xyz[cell])
            self.k[np.ix_(ids, ids)] += np.einsum('qai,ab,qbj->ij', b, d, b)*volume/4
            self.seed[ids] -= np.einsum('qai,a->i', b, s0)*volume/4
            self.seed[ids] += hydro_volume(self.xyz[cell], self.level, gamma)
            local = np.zeros((10, 3)); local[:, 2] = -material['unit_weight_kn_m3']*shape.sum(axis=0)*volume/4
            self.base[ids] += local.ravel()
            # Degree-two bulk sampling includes extrema at the vertices because
            # stress is affine and the elastic admissible set is convex.
            lam = np.vstack((np.eye(4), np.full((1, 4), .25)))
            _, check_b, check_x, _ = tet_geometry(self.xyz[cell], lam)
            self.volumes.append((ids, check_b, d, s0, mat, check_x))
        topo = topology(mesh)
        # Only explicitly prescribed geostatic base pressure is supported here.
        loads = geo.get('base_loads', {})
        if loads.get('nodal_loads') or geo.get('hydrostatic_face_sets') or cfg['surface_surcharge_kpa'] != 0:
            raise ValueError('unsupported exterior load policy')
        for load in loads.get('surface_loads', []):
            if load['kind'] != 'pressure' or np.any(load.get('gradient_kpa_per_m', [0, 0, 0])):
                raise ValueError('only constant exterior pressure supported')
            for face in mesh['face_sets'][load['face_set']]:
                (ei, nodes), = topo[tuple(sorted(face))]
                ids = (3*np.array(nodes)[:, None]+np.arange(3)).ravel()
                self.base[ids] += pressure_face(self.xyz, nodes, load['pressure_kpa']).ravel()
        if case['boundary'].get('nodal_loads') or case['boundary'].get('surface_loads'):
            raise ValueError('boundary loading must be displacement only')
        self.c = []; self.normal = []; self.weight = []; self.position = []; self.parent_face = []
        self.kn = []; self.kt = []; self.mu = []; pressure = []; self.water_load = np.zeros(nd)
        lam, weights = triangle_rule(4); shape, _ = basis(lam)
        paired_b = set()
        for definition in case['joint_system']['joints']:
            if definition.get('water', {}).get('kind') != 'continuous':
                raise ValueError('fixed CONTINUOUS hydrostatic scope only')
            for fi, face in enumerate(mesh['face_sets'][definition['side_a']]):
                (_, na), = topo[tuple(sorted(face))]
                pa, area_vec = face_geometry(self.xyz, na)
                matches = []
                for other in mesh['face_sets'][definition['side_b']]:
                    (_, nb), = topo[tuple(sorted(other))]
                    distance = np.linalg.norm(pa[:, None]-self.xyz[nb][None, :], axis=-1)
                    permutation = distance.argmin(axis=1)
                    if len(set(permutation)) == 6 and distance[np.arange(6), permutation].max() < 1e-10:
                        matches.append((tuple(sorted(other)), np.array(nb)[permutation]))
                if len(matches) != 1 or matches[0][0] in paired_b:
                    raise ValueError('nonunique paired face')
                key, nb = matches[0]; paired_b.add(key)
                normal = area_vec/np.linalg.norm(area_vec)
                if np.any(np.array(na) == nb):
                    raise ValueError('mechanical interface nodes must be separate')
                # This specimen has a horizontal, completely wet interface.
                if np.ptp(pa[:, 2]) > 1e-12 or pa[0, 2] >= self.level:
                    raise ValueError('only fully wet horizontal interface supported')
                face_p = gamma*(self.level-pa[0, 2])
                local_water = pressure_face(self.xyz, na, face_p)
                for node, other, value in zip(na, nb, local_water):
                    self.water_load[3*node:3*node+3] += value
                    self.water_load[3*other:3*other+3] -= value
                for n, w, point in zip(shape, weights, lam @ pa[:3]):
                    c = np.zeros((3, nd))
                    for a, bnode, value in zip(na, nb, n):
                        c[:, 3*a:3*a+3] -= value*np.eye(3)
                        c[:, 3*bnode:3*bnode+3] += value*np.eye(3)
                    self.c.append(c); self.normal.append(normal); self.position.append(point)
                    self.weight.append(w*np.linalg.norm(area_vec)); pressure.append(face_p)
                    self.parent_face.append(fi)
                    for name, key2 in (('kn', 'normal_stiffness_kpa_m'), ('kt', 'tangential_stiffness_kpa_m'), ('mu', 'friction_coefficient')):
                        getattr(self, name).append(definition['material'][key2])
        for name in ('c', 'normal', 'weight', 'position', 'parent_face', 'kn', 'kt', 'mu'):
            setattr(self, name, np.asarray(getattr(self, name)))
        self.pressure = np.asarray(pressure); self.base += self.water_load
        self.pore = np.concatenate([gamma*np.maximum(self.level-v[-1][:, 2], 0) for v in self.volumes])
        # Build an orthonormal nullspace independently at each constrained node.
        equations = {}
        for cond in case['boundary']['constraints']:
            direction = np.asarray(cond['direction'], dtype=float)
            for node in mesh['node_sets'][cond['node_set']]:
                a, b = equations.setdefault(node, ([], []))
                a.append(direction)
                b.append(cond['value_m']+np.dot(cond.get('gradient', [0, 0, 0]),
                                              self.xyz[node]-cond.get('reference_point_m', [0, 0, 0])))
        self.fixed = np.zeros(nd); columns = []
        for node in range(len(self.xyz)):
            if node in equations:
                a, b = map(np.asarray, equations[node]); value = np.linalg.lstsq(a, b, rcond=None)[0]
                if np.linalg.norm(a@value-b) > 1e-12:
                    raise ValueError('inconsistent displacement constraints')
                self.fixed[3*node:3*node+3] = value; space = null_space(a)
            else:
                space = np.eye(3)
            for col in space.T:
                v = np.zeros(nd); v[3*node:3*node+3] = col; columns.append(v)
        self.t = np.array(columns).T
        self.kr = self.t.T @ self.k @ self.t
        self.cr = self.c @ self.t

    def check_elastic(self, u):
        stress = []
        for ids, b, d, s0, mat, _ in self.volumes:
            s = s0-np.einsum('ab,qbi,i->qa', d, b, u[ids])
            tensor = np.zeros((len(s), 3, 3))
            tensor[:, range(3), range(3)] = s[:, :3]
            tensor[:, 0, 1] = tensor[:, 1, 0] = s[:, 3]
            tensor[:, 1, 2] = tensor[:, 2, 1] = s[:, 4]
            tensor[:, 0, 2] = tensor[:, 2, 0] = s[:, 5]
            principal = np.linalg.eigvalsh(tensor)
            phi = np.deg2rad(mat['friction_deg'])
            limit = 2*mat['cohesion_kpa']*np.cos(phi)
            yf = ((1-np.sin(phi))*principal[:, -1]-(1+np.sin(phi))*principal[:, 0]-limit)
            if np.max(yf) > 1e-8*max(limit, 1):
                raise ValueError('bulk yielded: elastic reference cannot continue')
            cutoff = mat.get('tension_cutoff_kpa')
            if cutoff is not None and np.min(principal) < -cutoff-1e-8:
                raise ValueError('bulk tensile yield outside scope')
            stress.extend(s)
        return np.asarray(stress)


def run(model, steps=16, local=8, tolerance=1e-10, sink=None):
    if not isinstance(steps, int) or steps < 2 or not isinstance(local, int) or local < 1:
        raise ValueError('positive partitions and at least two global steps required')
    if not 0 < tolerance <= 1e-10:
        raise ValueError('reference force tolerance may only be tightened')
    m = model; nq = len(m.weight)
    u = np.zeros(m.xyz.size); old_jump = np.zeros((nq, 3)); old_q = old_jump.copy()
    plastic = old_jump.copy(); accumulated = np.zeros(nq); rows = []; fields = []
    stages = [('geostatic', 0.)]+[(s['name'], s['target_factor']) for s in m.case['stages']]
    origin = 0.
    for phase, (stage, end) in enumerate(stages):
        for index in range(1, steps+1):
            progress = index/steps; factor = origin+(end-origin)*progress
            load = m.base+(1-progress)*(m.seed-m.base) if phase == 0 else m.base
            fixed = factor*m.fixed; free = m.t.T@(u-fixed)
            scale = max(float(np.linalg.norm(load)), 1.)

            def evaluate(x, jacobian=False):
                candidate = fixed+m.t@x; jump = np.einsum('qai,i->qa', m.c, candidate)
                targets = jump[:, None, :]
                if jacobian:
                    targets = np.repeat(targets, 7, axis=1)
                    for a in range(3):
                        targets[:, 1+2*a, a] += 2e-9; targets[:, 2+2*a, a] -= 2e-9
                q, dp, da = coulomb(old_jump, old_q, targets, m.normal, m.kn, m.kt, m.mu, local)
                residual = m.t.T@(m.k@candidate+m.seed-load)+np.einsum('qai,qa,q->i', m.cr, q[:, 0], m.weight)
                if jacobian:
                    tangent = np.stack([(q[:, 1+2*a]-q[:, 2+2*a])/(4e-9) for a in range(3)], axis=-1)
                    jac = m.kr+np.einsum('qai,qab,qbj,q->ij', m.cr, tangent, m.cr, m.weight, optimize=True)
                    return residual, jac
                return residual, candidate, jump, q[:, 0], dp[:, 0], da[:, 0]

            for iteration in range(60):
                residual, candidate, jump, q, dp, da = evaluate(free)
                relative = float(np.linalg.norm(residual)/scale)
                if relative < tolerance:
                    stress = m.check_elastic(candidate)
                    # Commit only after equilibrium AND the physical-scope check.
                    u = candidate; old_jump = jump; old_q = q; plastic += dp; accumulated += da
                    row = dict(stage=stage, progress=progress, factor=factor, iterations=iteration,
                               relative_force_residual=relative, local_partitions=local)
                    value = dict(u=u.reshape(-1, 3).copy(), jump=jump.copy(), traction=q.copy(),
                                 plastic_slip=plastic.copy(), accumulated=accumulated.copy(),
                                 gap=(jump*m.normal).sum(axis=1), stress=stress,
                                 physical_pressure=m.pressure.copy(), pore=m.pore.copy(),
                                 water_load=m.water_load.reshape(-1, 3).copy(),
                                 slip_increment=da.copy())
                    rows.append(row); fields.append(value)
                    if sink is not None:
                        sink(row, value)
                    break
                residual, jac = evaluate(free, True)
                correction = linear_correction(jac, -residual)
                for cut in range(18):
                    trial = free+correction*2.**(-cut)
                    if np.linalg.norm(evaluate(trial)[0]) < np.linalg.norm(residual):
                        free = trial; break
                else:
                    raise RuntimeError(('line search', stage, progress, iteration, relative))
            else:
                raise RuntimeError(('Newton limit', stage, progress, relative))
        origin = end
    return dict(rows=rows, fields={key: np.stack([f[key] for f in fields]) for key in fields[0]},
                metadata=dict(steps=steps, local=local, force_tolerance=tolerance,
                              reference_uses_native_results=False, complete=True,
                              numerical_acceptance=False, nodes=len(m.xyz), elements=len(m.mesh['tetrahedra']),
                              contact_position=m.position.tolist(), contact_weight=m.weight.tolist(),
                              normal=m.normal.tolist(), node_position=m.xyz.tolist(), node_body=m.node_body, water_level=m.level))
