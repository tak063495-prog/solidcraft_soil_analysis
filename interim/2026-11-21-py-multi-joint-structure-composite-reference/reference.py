"""Independent dry small-strain simplex/contact reference (NumPy/SciPy only).

No production modules, executable, stored production solution or water model.
The complete model and raw uncorrected fields accompany every result.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, permutations, product
import json
from pathlib import Path
import time
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import solve

SCHEMA = 'dry-multi-joint-reference/1'
UNITS = {'length': 'm', 'force': 'kN', 'stress': 'kPa', 'work': 'kN m',
         'moment': 'kN m', 'angle': 'rad'}
SIGNS = {'jump': 'u_plus-u_minus', 'normal': 'minus_to_plus',
         'traction': 'internal_on_plus; compression_negative_normal',
         'reaction': 'internal-applied', 'axial_force': 'tension_positive',
         'moment': 'right_hand_global_xyz'}
CASES = ('parallel2d', 'cross2d', 't2d', 'shared2d', 'parallel3d',
         'cross3d', 't3d', 'shared3d', 'beam2d', 'anchor3d')
# Operational Newton controls only, NOT new numerical acceptance criteria.
NEWTON_FORCE_KN = 2e-11
MAX_ITERATIONS = 60


def dump(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2,
                                   allow_nan=False)+'\n', encoding='utf-8')


def elasticity(dim, young=3000., poisson=.25):
    mu = young/(2*(1+poisson))
    lam = young*poisson/((1+poisson)*(1-2*poisson))
    c = np.zeros((3 if dim == 2 else 6,)*2)
    c[:dim, :dim] = lam
    for i in range(dim):
        c[i, i] += 2*mu
    for i in range(dim, len(c)):
        c[i, i] = mu
    return c


def simplex(x):
    """Affine shape gradients, engineering-strain B, exact volume/area."""
    dim = x.shape[1]
    mat = np.column_stack((np.ones(dim+1), x))
    grad = np.linalg.inv(mat)[1:, :].T
    volume = abs(np.linalg.det(x[1:]-x[0]))/(2 if dim == 2 else 6)
    if volume <= 0:
        raise ValueError('degenerate simplex')
    b = np.zeros((3 if dim == 2 else 6, dim*(dim+1)))
    for a, g in enumerate(grad):
        for i in range(dim):
            b[i, dim*a+i] = g[i]
        for k, (i, j) in enumerate(((0, 1),) if dim == 2 else ((0, 1), (1, 2), (0, 2))):
            b[dim+k, dim*a+i] = g[j]
            b[dim+k, dim*a+j] = g[i]
    return b, volume


def face_rule(x, order):
    """GL on segment, positive Duffy-product quadrature on Tri3."""
    q, w = leggauss(order)
    q, w = (q+1)/2, w/2
    if len(x) == 2:
        length = np.linalg.norm(x[1]-x[0])
        return [(np.array([1-a, a]), length*wa) for a, wa in zip(q, w)]
    jac = np.linalg.norm(np.cross(x[1]-x[0], x[2]-x[0]))
    return [(np.array([(1-a)*(1-b), a, (1-a)*b]), jac*(1-a)*wa*wb)
            for a, wa in zip(q, w) for b, wb in zip(q, w)]


@dataclass
class Model:
    name: str
    x: np.ndarray
    node_ids: list
    cells: list
    contacts: list
    members: list
    fixed: np.ndarray
    bulk_k: np.ndarray
    member_k: np.ndarray
    dim: int
    metadata: dict


def make_model(name='cross2d', gp=2, mesh=2, variant='identity'):
    if name not in CASES or gp not in range(1, 9) or mesh not in range(2, 13):
        raise ValueError('unsupported case/GP/mesh; no arbitrary CAD or water input')
    if name.startswith(('beam', 'anchor')) and mesh % 2:
        raise ValueError('axial attachment coordinates require an even mesh; no snapping')
    dim = 3 if name.endswith('3d') else 2
    topology = 'parallel' if name.startswith(('parallel', 'beam', 'anchor')) else name[:-2]
    boxes = ([(0, 1, 0, 1), (1, 2, 0, 1), (2, 3, 0, 1)] if topology == 'parallel'
             else [(0, 1, 0, 1), (1, 2, 0, 1), (0, 1, 1, 2), (1, 2, 1, 2)])
    if topology == 't':
        boxes = [(0, 2, 0, 1), (0, 1, 1, 2), (1, 2, 1, 2)]
    shared = topology == 'shared'
    coords, node_ids, cells, faces = [], [], [], defaultdict(list)
    common = {}
    for body, box in enumerate(boxes):
        lo = np.array([box[0], box[2]]+([0] if dim == 3 else []), float)
        hi = np.array([box[1], box[3]]+([1] if dim == 3 else []), float)
        count = np.rint((hi-lo)*mesh).astype(int)
        local = {}
        for ij in product(*(range(n+1) for n in count)):
            point = lo+np.array(ij)/mesh
            key = tuple(point)
            merge = shared and point[0] == 1 and point[1] == 1
            if merge and key in common:
                local[ij] = common[key]
            else:
                local[ij] = len(coords)
                coords.append(point)
                node_ids.append(f'B{body}:'+','.join(map(str, ij)))
                if merge:
                    common[key] = local[ij]
        body_cells = []
        for ij in product(*(range(n) for n in count)):
            for perm in permutations(range(dim)):
                index = np.array(ij)
                indices = [local[tuple(index)]]
                for axis in perm:
                    index = index.copy(); index[axis] += 1
                    indices.append(local[tuple(index)])
                body_cells.append(indices)
        counts = Counter(tuple(sorted(f)) for c in body_cells for f in combinations(c, dim))
        for f, count_f in counts.items():
            if count_f != 1:
                continue
            key = tuple(sorted(tuple(coords[i]) for i in f))
            faces[key].append((body, f))
        cells.extend(body_cells)
    x = np.array(coords)
    raw_contacts = []
    for key, owners in sorted(faces.items()):
        if len(owners) != 2:
            continue
        (bm, fm), (bp, fp) = owners
        # Canonical geometric order, independent of local/global numbering.
        fm = sorted(fm, key=lambda i: tuple(x[i])); fp = sorted(fp, key=lambda i: tuple(x[i]))
        cm = np.array([(boxes[bm][0]+boxes[bm][1])/2, (boxes[bm][2]+boxes[bm][3])/2])
        cp = np.array([(boxes[bp][0]+boxes[bp][1])/2, (boxes[bp][2]+boxes[bp][3])/2])
        face_x = x[fm]
        axis = next(a for a in range(2) if np.ptp(face_x[:, a]) == 0)
        normal = np.zeros(dim); normal[axis] = np.sign(cp[axis]-cm[axis])
        joint_id = f'J{"xy"[axis]}={face_x[0, axis]:g}'
        facet_id = '|'.join(','.join(f'{v:g}' for v in p) for p in face_x)
        for gi, (shape, weight) in enumerate(face_rule(face_x, gp)):
            raw_contacts.append(dict(joint_id=joint_id, gp_id=facet_id+f'/q{gi}',
                                     minus=fm, plus=fp, shape=shape, weight=weight,
                                     normal=normal, position=shape@face_x))
    if variant == 'nodes':
        order = np.random.default_rng(20261121).permutation(len(x))
        inverse = np.argsort(order)
        x = x[order]; node_ids = [node_ids[i] for i in order]
        cells = [[int(inverse[i]) for i in c] for c in cells]
        for p in raw_contacts:
            p['minus'] = [int(inverse[i]) for i in p['minus']]
            p['plus'] = [int(inverse[i]) for i in p['plus']]
    nd = dim*len(x)
    k = np.zeros((nd, nd)); cmat = elasticity(dim)
    element_data = []
    for ci, c in enumerate(cells):
        dofs = np.array([dim*n+a for n in c for a in range(dim)])
        b, volume = simplex(x[c])
        k[np.ix_(dofs, dofs)] += volume*b.T@cmat@b
        element_data.append(dict(id=f'E{ci}', nodes=c, dofs=dofs, b=b, volume=volume))
    for p in raw_contacts:
        ids = np.array(sorted({dim*n+a for n in p['minus']+p['plus'] for a in range(dim)}))
        jump = np.zeros((dim, len(ids)))
        lookup = {int(v): i for i, v in enumerate(ids)}
        for side, sign in [('minus', -1), ('plus', 1)]:
            for n, sh in zip(p[side], p['shape']):
                for a in range(dim):
                    jump[a, lookup[dim*n+a]] += sign*sh
        p.update(dofs=ids, b=jump)
    members = []; mk = np.zeros_like(k)
    if name.startswith(('beam', 'anchor')):
        # Two linear axial members bridge the joints at shared interface nodes.
        # Nearby attachment uses the existing body's elastic triangles/tetrahedra,
        # never nearest-node snapping or a hidden rigid link.
        for mi, yy in enumerate((.5, 0.)):
            a = next(i for i, p in enumerate(x) if tuple(p) == tuple([1., yy]+([.5] if dim == 3 else [])) and node_ids[i].startswith('B0:'))
            bnode = next(i for i, p in enumerate(x) if tuple(p) == tuple([2., yy]+([.5] if dim == 3 else [])) and node_ids[i].startswith('B2:'))
            vec = x[bnode]-x[a]; length = np.linalg.norm(vec); direction = vec/length
            dofs = np.array([dim*n+i for n in (a, bnode) for i in range(dim)])
            b = np.r_[-direction, direction]
            ea = 600. if name.startswith('beam') else 300.
            mk[np.ix_(dofs, dofs)] += ea/length*np.outer(b, b)
            members.append(dict(id=f'M{mi}', kind='beam_axial' if name.startswith('beam') else 'anchor_axial',
                                nodes=[a, bnode], dofs=dofs, b=b, ea=ea, length=length))
    lo, hi = x.min(axis=0), x.max(axis=0)
    boundary = np.any((x == lo)|(x == hi), axis=1)
    fixed = np.flatnonzero(np.repeat(boundary, dim))
    if variant == 'contacts':
        raw_contacts.reverse()
    elif variant == 'joints':
        raw_contacts.sort(key=lambda p: p['joint_id'], reverse=True)
    elif variant not in ('identity', 'nodes'):
        raise ValueError('unknown order variant')
    metadata = dict(case=name, dimension=dim, bulk='Tri3_plane_strain' if dim == 2 else 'Tet4',
                    E_kpa=3000., nu=.25, kn_kpa_per_m=50000., kt_kpa_per_m=30000.,
                    mu=.25, mesh_per_m=mesh, gp_order=gp, thickness_m=1., generation=0,
                    water='absent', lifecycle='fixed; no BIRTH/DEATH/projection',
                    structure='linear axial only; no plate/shell/offset/rotation DOFs')
    return Model(name, x, node_ids, element_data, raw_contacts, members, fixed, k, mk, dim, metadata)


def contact_return(jump, old_slip, normal, kn, kt, mu):
    """Unilateral penalty + Coulomb radial return; no smoothing or force clipping.

    The min/branch is the stated constitutive law, never an output correction.
    Trial history is read-only throughout Newton; commit only after equilibrium.
    """
    gap = float(normal@jump)
    eye_t = np.eye(len(jump))-np.outer(normal, normal)
    jt = eye_t@jump
    if gap >= 0:
        return np.zeros_like(jump), np.zeros((len(jump),)*2), old_slip.copy(), 'open', 0.
    pn = -kn*gap
    trial = kt*(jt-old_slip); norm = np.linalg.norm(trial)
    limit = mu*pn
    if norm <= limit:
        tt = trial; new_slip = old_slip.copy(); state = 'stick'
        tangent = kn*np.outer(normal, normal)+kt*eye_t
    else:
        direction = trial/norm; tt = limit*direction
        new_slip = jt-tt/kt; state = 'slip'
        tangent = (kn*np.outer(normal, normal)-mu*kn*np.outer(direction, normal)
                   +limit*kt/norm*(eye_t-np.outer(direction, direction)))
    traction = -pn*normal+tt
    return traction, tangent, new_slip, state, float(.5*(kn*gap*gap+tt@tt/kt))


def assemble(model, u, history, tangent=True):
    k = model.bulk_k+model.member_k
    force = k@u
    jac = k.copy() if tangent else None
    values = []
    for p in model.contacts:
        key = (p['joint_id'], p['gp_id'])
        old = history[key]
        # Pair differences first preserve an exactly identical displacement on
        # coincident sides without subtracting two large interpolated sums.
        # This is algebraic assembly, not a zeroing/rounding correction.
        nodal = u.reshape(-1, model.dim)
        jump = p['shape']@(nodal[p['plus']]-nodal[p['minus']])
        q, c, slip, state, energy = contact_return(jump, old, p['normal'], 50000., 30000., .25)
        force[p['dofs']] += p['weight']*p['b'].T@q
        if tangent:
            jac[np.ix_(p['dofs'], p['dofs'])] += p['weight']*p['b'].T@c@p['b']
        values.append((jump, q, slip, state, energy))
    return force, jac, values


def boundary_displacement(model, phase, t):
    x = model.x
    u = -.004*x
    if phase == 0:
        return (u*t).ravel()
    # preload -> shear -> reversal; every phase is continuous.
    factor = t if phase == 1 else 1-1.5*t
    u[:, 0] += .032*factor*x[:, 1]
    u[:, 1] += .022*factor*x[:, 0]
    if model.dim == 3:
        u[:, 2] += .013*factor*x[:, 0]
    return u.ravel()


def moments(x, f):
    if x.shape[1] == 2:
        return [0., 0., float(np.sum(x[:, 0]*f[:, 1]-x[:, 1]*f[:, 0]))]
    return np.cross(x, f).sum(axis=0).tolist()


def audit_frame(model, u, values, applied, history, accumulated, step, phase, progress,
                external_work, dissipation, previous_u, previous_force, previous_values):
    """Reassemble stress/traction/axial end forces, separate from matrix residual."""
    dim = model.dim; nd = len(u)
    bulk = np.zeros(nd); cf = np.zeros(nd); sf = np.zeros(nd)
    energy_bulk = 0.; element_rows = []
    cmat = elasticity(dim)
    for e in model.cells:
        strain = e['b']@u[e['dofs']]; stress = cmat@strain
        bulk[e['dofs']] += e['volume']*e['b'].T@stress
        energy_bulk += .5*e['volume']*strain@stress
        element_rows.append(dict(id=e['id'], strain=strain.tolist(), stress_kpa=stress.tolist()))
    member_rows = []; energy_member = 0.
    for m in model.members:
        extension = m['b']@u[m['dofs']]
        force = m['ea']/m['length']*extension
        ends = m['b']*force
        sf[m['dofs']] += ends
        energy = .5*force*extension; energy_member += energy
        member_rows.append(dict(id=m['id'], kind=m['kind'], node_ids=[model.node_ids[i] for i in m['nodes']],
                                axial_force_kn=float(force), end_internal_kn=ends.reshape(2, dim).tolist(),
                                end_moment_kn_m=[[0., 0., 0.]]*2, elastic_energy_kn_m=float(energy)))
    contact_rows = []; energy_contact = 0.; local_work_error = 0.
    increment_diss = 0.
    for i, (p, value) in enumerate(zip(model.contacts, values)):
        jump, q, slip, state, energy = value
        key = (p['joint_id'], p['gp_id'])
        # Nodal scatter is explicit (not reused from assemble or its B matrix).
        for side, sign in [('minus', -1), ('plus', 1)]:
            for n, sh in zip(p[side], p['shape']):
                cf[dim*n:dim*(n+1)] += sign*p['weight']*sh*q
        ds = slip-history[key]
        delta = np.linalg.norm(ds)
        pn = -q@p['normal']
        dd = p['weight']*.25*pn*delta
        increment_diss += dd
        old_jump, old_q, _, _, old_energy = previous_values[i]
        work = p['weight']*.5*(old_q+q)@(jump-old_jump)
        error = work-p['weight']*(energy-old_energy)-dd
        local_work_error += error
        energy_contact += p['weight']*energy
        contact_rows.append(dict(joint_id=p['joint_id'], gp_id=p['gp_id'], generation=0,
                                 position_m=p['position'].tolist(), normal=p['normal'].tolist(),
                                 weight_m2=float(p['weight']), jump_m=jump.tolist(),
                                 traction_kpa=q.tolist(), force_kn=(p['weight']*q).tolist(),
                                 plastic_slip_m=slip.tolist(), accumulated_slip_m=float(accumulated[key]+delta),
                                 gap_m=float(jump@p['normal']), state=state, open=state == 'open',
                                 elastic_energy_kn_m=float(p['weight']*energy),
                                 friction_increment_kn_m=float(dd), work_increment_kn_m=float(work),
                                 work_defect_increment_kn_m=float(error)))
    total = bulk+cf+sf
    reaction = np.zeros(nd); reaction[model.fixed] = (total-applied)[model.fixed]
    residual = total-applied-reaction
    energy = energy_bulk+energy_contact+energy_member
    external = external_work+.5*(previous_force+total)@(u-previous_u)
    dissip = dissipation+increment_diss
    node_rows = [dict(id=nid, position_m=model.x[i].tolist(), displacement_m=u[dim*i:dim*(i+1)].tolist(),
                      bulk_internal_kn=bulk[dim*i:dim*(i+1)].tolist(), contact_internal_kn=cf[dim*i:dim*(i+1)].tolist(),
                      structure_internal_kn=sf[dim*i:dim*(i+1)].tolist(), applied_kn=applied[dim*i:dim*(i+1)].tolist(),
                      reaction_kn=reaction[dim*i:dim*(i+1)].tolist(), residual_kn=residual[dim*i:dim*(i+1)].tolist())
                 for i, nid in enumerate(model.node_ids)]
    independent_force, _, _ = assemble(model, u, history, False)
    diagnostics = dict(max_node_residual_kn=float(np.max(abs(residual))),
                       matrix_vs_stress_scatter_kn=float(np.max(abs(total-independent_force))),
                       resultant_kn=(applied+reaction).reshape(-1, dim).sum(axis=0).tolist(),
                       external_moment_kn_m=moments(model.x, (applied+reaction).reshape(-1, dim)),
                       contact_resultant_kn=cf.reshape(-1, dim).sum(axis=0).tolist(),
                       contact_moment_kn_m=moments(model.x, cf.reshape(-1, dim)),
                       member_resultant_kn=sf.reshape(-1, dim).sum(axis=0).tolist(),
                       member_moment_kn_m=moments(model.x, sf.reshape(-1, dim)),
                       work_defect_kn_m=float(external-energy-dissip),
                       local_work_defect_increment_kn_m=float(local_work_error),
                       energy_matrix_difference_kn_m=float(energy_bulk+energy_member-.5*u@(model.bulk_k+model.member_k)@u))
    return dict(step=step, phase=phase, progress=progress, nodes=node_rows, contacts=contact_rows,
                members=member_rows, elements=element_rows, diagnostics=diagnostics,
                external_work_kn_m=float(external), bulk_energy_kn_m=float(energy_bulk),
                contact_energy_kn_m=float(energy_contact), member_energy_kn_m=float(energy_member),
                friction_dissipation_kn_m=float(dissip)), total


def run(name='cross2d', steps=16, gp=2, mesh=2, variant='identity'):
    if steps < 2 or steps > 4096:
        raise ValueError('steps outside explicit bounded range')
    model = make_model(name, gp, mesh, variant)
    dim = model.dim; nd = dim*len(model.x)
    free = np.setdiff1d(np.arange(nd), model.fixed)
    history = {(p['joint_id'], p['gp_id']): np.zeros(dim) for p in model.contacts}
    accumulated = {key: 0. for key in history}
    u = np.zeros(nd); previous_u = u.copy(); previous_force = u.copy()
    previous_values = assemble(model, u, history, False)[2]
    work = diss = 0.; frames = []; max_iterations = 0
    started = time.monotonic()
    for phase in range(3):
        for increment in range(1, steps+1):
            t = increment/steps
            target = boundary_displacement(model, phase, t)
            u[model.fixed] = target[model.fixed]
            for iteration in range(MAX_ITERATIONS):
                force, jac, values = assemble(model, u, history)
                residual = force[free]
                if np.max(abs(residual), initial=0.) <= NEWTON_FORCE_KN:
                    break
                du = solve(jac[np.ix_(free, free)], -residual, assume_a='gen', check_finite=True)
                base = np.linalg.norm(residual)
                for backtrack in range(18):
                    candidate = u.copy(); candidate[free] += 2.**(-backtrack)*du
                    trial = assemble(model, candidate, history, False)[0][free]
                    if np.linalg.norm(trial) < base:
                        u = candidate
                        break
                else:
                    raise RuntimeError(f'line search failed: {name}/{phase}/{increment}; residual={base}')
            else:
                raise RuntimeError('Newton iteration budget exhausted; no partial success')
            max_iterations = max(max_iterations, iteration+1)
            frame, total = audit_frame(model, u, values, np.zeros(nd), history, accumulated,
                                       len(frames)+1, phase, t, work, diss, previous_u,
                                       previous_force, previous_values)
            frames.append(frame)
            work, diss = frame['external_work_kn_m'], frame['friction_dissipation_kn_m']
            for p, value, row in zip(model.contacts, values, frame['contacts']):
                key = (p['joint_id'], p['gp_id'])
                history[key] = value[2].copy(); accumulated[key] = row['accumulated_slip_m']
            previous_u, previous_force, previous_values = u.copy(), total.copy(), values
    return dict(schema=SCHEMA, producer='independent_numpy_scipy', complete=True, accepted=None,
                acceptance_reason='A-01 expansion criteria undefined; diagnostic reference only',
                units=UNITS, signs=SIGNS, transformations=[], raw=True, metadata=model.metadata,
                node_ids=model.node_ids, joint_ids=sorted({p['joint_id'] for p in model.contacts}),
                member_ids=[m['id'] for m in model.members],
                evaluation_order=[p['joint_id']+'/'+p['gp_id'] for p in model.contacts],
                bulk_connectivity=[[model.node_ids[i] for i in e['nodes']] for e in model.cells],
                frames=frames, operational=dict(newton_force_kn=NEWTON_FORCE_KN,
                                               max_iterations_observed=max_iterations),
                elapsed_seconds=time.monotonic()-started)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case', choices=CASES, default='cross2d')
    p.add_argument('--steps', type=int, default=16)
    p.add_argument('--gp', type=int, default=2)
    p.add_argument('--mesh', type=int, default=2)
    p.add_argument('--variant', choices=('identity', 'contacts', 'joints', 'nodes'), default='identity')
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    result = run(args.case, args.steps, args.gp, args.mesh, args.variant)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    dump(args.output, result)
    print(json.dumps(dict(case=args.case, frames=len(result['frames']),
                          elapsed_seconds=result['elapsed_seconds'], accepted=None)))


if __name__ == '__main__':
    main()
