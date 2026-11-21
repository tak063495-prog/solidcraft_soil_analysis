"""Instrumented existing input-only wet JOINT reference (verification only).

Restricted to fixed hydrostatic water and boundary-factor cycles. No native
result or accepted increment sequence is an argument to the reference solve.
The existing independent Coulomb projector is shared with the dry verifier.
"""
import copy
from types import SimpleNamespace
import numpy as np
import fem_reference as fem
from joint_water_reference import Water
from mc_reference import normalized_yield
from verify_high_order import physical
from verify_joint_coupling import contact


def contact_work(system, old_jump, old_q, target, count):
    """Endpoint work on closed subsegments, including exact plane crossings.

    Numerical work depends on the integration partition and is not a
    continuum dissipation. Keep it separate from plastic dissipation.
    """
    jump = old_jump.copy(); q = old_q.copy(); n = system.normal
    delta = (target-old_jump)/count
    plastic_work = np.zeros(len(jump)); numerical = np.zeros(len(jump)); work = np.zeros(len(jump))
    def energy(q):
        qn = np.einsum('qa,qa->q', q, n); tau = q-qn[:, None]*n
        return qn**2/(2*system.kn)+np.einsum('qa,qa->q', tau, tau)/(2*system.kt)
    for step in range(1, count+1):
        end = old_jump+step*delta if step < count else target
        g0 = np.einsum('qa,qa->q', jump, n); g1 = np.einsum('qa,qa->q', end, n)
        fraction = np.ones(len(jump)); fraction[(g0 > 0) & (g1 > 0)] = 0
        opening = (g0 <= 0) & (g1 > 0); closing = (g0 > 0) & (g1 <= 0)
        np.divide(-g0, g1-g0, out=fraction, where=opening)
        np.divide(g1, g1-g0, out=fraction, where=closing)
        updated, flow, _ = contact(system, jump, q, end, 1); updated = updated[:, 0]
        dw = np.einsum('qa,qa->q', updated, (end-jump)*fraction[:, None])
        pd = np.maximum(np.einsum('qa,qa->q', updated, flow[:, 0]), 0)
        nd = dw-(energy(updated)-energy(q))-pd
        assert np.all(nd >= -1e-10*(1+np.abs(dw)+energy(q))), 'independent work identity'
        plastic_work += pd; numerical += np.maximum(nd, 0); work += dw
        jump = end; q = updated
    return plastic_work, numerical, work


def reference(mesh, ground, case, steps=32, contact_steps=20, order=5):
    return _reference(mesh, ground, case, steps, contact_steps, order, False)


def discrete_volume_reference(mesh, ground, case, steps=32, contact_steps=20):
    """Input-only reference for the specified Gauss-Jacobi volume discretization.

    This isolates path integration from partial-wetting quadrature error. It
    is not a continuum-volume reference or a new acceptance criterion.
    """
    return _reference(mesh, ground, case, steps, contact_steps, mesh.get('quadrature_order', 3), True)


def cut_volume_reference(mesh, ground, case, steps=32, contact_steps=20, order=5, history_sink=None):
    """Affine-cell exact hydrostatic force; retain the existing GL bulk rule.

    Only the pore-volume term changes. This is a verification reference, not
    production adoption or a certificate of contact-history convergence.
    """
    return _reference(mesh, ground, case, steps, contact_steps, order, False, True, history_sink)


def _reference(mesh, ground, case, steps, contact_steps, order, discrete_volume, cut_volume=False, history_sink=None, *, maximum_m=None, maximum_dofs=300, progress_grids=None):
    assert isinstance(steps, int) and steps >= 2
    assert isinstance(contact_steps, int) and contact_steps >= 1
    assert order >= 3
    assert maximum_m is None or (np.isfinite(maximum_m) and maximum_m > 0)
    assert case['settings']['kinematics'] == 'explicit'
    assert all(s['kind'] == 'boundary' for s in case['stages']), 'boundary cycle only'
    assert not any(case.get(k) for k in ('beam_system', 'anchor_system', 'initially_inactive_bodies'))
    cfg = ground['initialization']['geostatic']
    assert cfg['seed']['kind'] == 'body_states' and cfg['initial_boundary_factor'] == 0
    assert ground['initialization']['method'] == 'equilibrate'
    assert all(j.get('water', {}).get('kind') == 'continuous'
               for j in case['joint_system']['joints']), 'fixed CONTINUOUS hydrostatic scope'
    assert not case['joint_system'].get('structure_bodies')
    xyz = np.asarray(mesh['nodes']); nd = xyz.size
    assert nd <= maximum_dofs <= 15000, 'bounded verification specimen only'
    active = np.ones(len(mesh['tetrahedra']), dtype=bool)
    water = Water(mesh, ground, case); system = water.contact
    materials = {m['id']: m['mohr_coulomb'] for m in ground['materials']}
    regions = {r['id']: r['body'] for r in mesh['regions']}
    level = ground['initialization']['water_level_m']; gamma = water.gamma
    assert level is not None
    k = np.zeros((nd, nd)); seed_force = np.zeros(nd); operators = []
    for ei, mid in enumerate(mesh['material_ids']):
        mat = materials[mid]; d = fem.elastic_tensor(mat['young_kpa'], mat['poisson'])
        seed = cfg['seed']['states'][regions[mesh['region_ids'][ei]]]
        assert seed['equivalent_plastic_strain'] == seed['plastic_work_kpa'] == 0
        assert not np.asarray(seed['plastic_strain_compression_positive']).any()
        s0 = np.asarray(seed['effective_stress_compression_positive_kpa'])
        q = fem.volume(mesh, ei, independent=not discrete_volume, independent_order=order)
        ids = (3*q['nodes'][:, None]+np.arange(3)).ravel()
        k[np.ix_(ids, ids)] += np.einsum('qai,ab,qbj,q->ij', q['b'], d, q['b'], q['weight'], optimize=True)
        pore = gamma*np.maximum(level-q['position'][:, 2], 0)
        total = np.broadcast_to(s0, (len(pore), 6)).copy()
        if not cut_volume:
            total[:, :3] += pore[:, None]
        seed_force[ids] -= np.einsum('qai,qa,q->i', q['b'], total, q['weight'])
        if cut_volume:
            from hydrostatic_cut_volume_reference import force
            # Replace just the pressure term. Material/contact rules stay intact.
            seed_force[ids] += force(mesh, ei, level, gamma)['force'].ravel()
        operators.append((fem.volume(mesh, ei), ids, d, s0, seed, mat))
    basis, error = fem.staged_constraint_space(mesh, case, dict(
        active_elements=active, boundary_factor=0., displacement_m=np.zeros_like(xyz)))
    assert error == 0
    nf = sum(b.shape[1] for b in basis.values()); t = np.zeros((nd, nf)); start = 0
    for node, b in basis.items():
        t[3*node:3*node+3, start:start+b.shape[1]] = b; start += b.shape[1]
    kr = t.T@k@t; cr = np.einsum('qag,gi->qai', system.c, t)
    interface = water.hydrostatic(active, level)
    base, _ = physical(mesh, ground, active, ground['initialization']['surface_surcharge_kpa'],
                       case.get('pressure_on_original_surface_only', False), joint_system=system)
    base = (base+interface['load']).ravel()
    boundary = fem.conditions_load(mesh, case['boundary'], active).ravel()
    holding = seed_force-base
    u = np.zeros(nd); jump = np.zeros((len(system.points), 3)); traction = np.zeros_like(jump)
    plastic = np.zeros_like(jump); accumulated = np.zeros(len(jump))
    plastic_work = np.zeros(len(jump)); numerical = np.zeros(len(jump)); endpoint_work = np.zeros(len(jump))
    maximum_residual = 0.; path = []; frames = []

    def offset(factor):
        equations = {}
        for cond in case['boundary']['constraints']:
            direction = np.asarray(cond['direction'])
            gradient = np.asarray(cond.get('gradient', [0., 0., 0.]))
            for node in mesh['node_sets'][cond['node_set']]:
                a, b = equations.setdefault(node, ([], [])); a.append(direction)
                b.append(factor*(cond['value_m']+gradient@(xyz[node]-cond.get('reference_point_m', [0., 0., 0.]))))
        fixed = np.zeros_like(xyz)
        for node, (a, b) in equations.items():
            fixed[node] = np.linalg.lstsq(a, b, rcond=None)[0]
            np.testing.assert_allclose(np.asarray(a)@fixed[node], b, atol=1e-12, rtol=0)
        return fixed.ravel()

    def partition(target):
        if maximum_m is None:
            return np.full(len(jump), contact_steps, dtype=int)
        delta = target-jump[:, None, :] if target.ndim == 3 else target-jump
        length = np.linalg.norm(delta, axis=-1)
        if length.ndim == 2:
            length = length.max(axis=1)
        counts = np.maximum(1, np.ceil(length/maximum_m).astype(int))
        assert np.max(counts) <= 50000, 'existing JOINT partition budget exceeded'
        return counts

    def subset(mask):
        return SimpleNamespace(**{k: getattr(system, k)[mask] for k in ('normal', 'kn', 'kt', 'mu')})

    frozen_counts = np.ones(len(jump), dtype=int)

    def project(target):
        if maximum_m is None:
            return contact(system, jump, traction, target, contact_steps)
        counts = frozen_counts
        shape = target.shape if target.ndim == 3 else (len(jump), 1, 3)
        q = np.zeros(shape); flow = np.zeros(shape); cum = np.zeros(shape[:2])
        for count in np.unique(counts):
            mask = counts == count
            q[mask], flow[mask], cum[mask] = contact(subset(mask), jump[mask], traction[mask], target[mask], int(count))
        return q, flow, cum

    def work_parts(target):
        if maximum_m is None:
            return contact_work(system, jump, traction, target, contact_steps)
        counts = frozen_counts
        pd = np.zeros(len(jump)); nd = pd.copy(); work = pd.copy()
        for count in np.unique(counts):
            mask = counts == count
            pd[mask], nd[mask], work[mask] = contact_work(subset(mask), jump[mask], traction[mask], target[mask], int(count))
        return pd, nd, work

    def increment(load, factor, name, progress):
        nonlocal u, jump, traction, plastic, accumulated, maximum_residual, frozen_counts
        frozen_counts = np.ones(len(jump), dtype=int)
        fixed = offset(factor); free = t.T@(u-fixed); scale = max(np.linalg.norm(load), 1.)
        def evaluate(free):
            candidate = fixed+t@free; target = np.einsum('qag,g->qa', system.c, candidate)
            q, slip, cum = project(target)
            force = k@candidate+seed_force+np.einsum('qag,qa,q->g', system.c, q[:, 0], system.weight)
            return candidate, target, q[:, 0], slip[:, 0], cum[:, 0], t.T@(force-load)
        for iteration in range(60):
            candidate, target, q, slip, cum, residual = evaluate(free)
            relative = float(np.linalg.norm(residual)/scale)
            if relative < 1e-10:
                if maximum_m is not None and np.any(partition(target) > frozen_counts):
                    frozen_counts = np.maximum(frozen_counts, partition(target))
                    continue
                pd, nd, work = work_parts(target)
                plastic_work[:] += pd; numerical[:] += nd; endpoint_work[:] += work
                u = candidate; jump = target; traction = q; plastic += slip; accumulated += cum
                maximum_residual = max(maximum_residual, relative)
                path.append(dict(stage=name, progress=progress, factor=factor, relative_force_residual=relative, newton_iterations=iteration, acceptance_force_tolerance=1e-10,
                                 contact_partition_counts=frozen_counts.tolist() if maximum_m is not None else [contact_steps]*len(jump),
                                 open_points=int((np.einsum('qa,qa->q', jump, system.normal) > 0).sum()),
                                 slipping_points=int((cum > 0).sum())))
                if history_sink is not None:
                    history_sink(dict(path[-1]), field_state(name, load))
                return
            h = 2e-9; perturb = np.repeat(target[:, None, :], 7, axis=1)
            for axis in range(3):
                perturb[:, 1+2*axis, axis] += h; perturb[:, 2+2*axis, axis] -= h
            pq, _, _ = project(perturb)
            tangent = np.stack([(pq[:, 1+2*a]-pq[:, 2+2*a])/(2*h) for a in range(3)], axis=2)
            jac = kr+np.einsum('qai,qab,qbj,q->ij', cr, tangent, cr, system.weight, optimize=True)
            correction = np.linalg.solve(jac, -residual)
            for power in range(16):
                trial = free+correction/2**power
                if np.linalg.norm(evaluate(trial)[-1]) < np.linalg.norm(residual):
                    free = trial; break
            else:
                raise AssertionError(('independent line search', name, progress, iteration, relative))
        raise AssertionError(('independent Newton limit', name, progress))

    def field_state(name, load):
        stress = []; strain = []; void = []; pore = []
        for q, ids, d, s0, initial, mat in operators:
            e = np.einsum('qai,i->qa', q['b'], u[ids]); s = s0-e@d.T
            assert np.max(normalized_yield(mat, s)) <= 1e-8, 'bulk plasticity outside reference'
            stress.extend(s); strain.extend(e)
            void.extend((1+initial['void_ratio'])*np.exp(e[:, :3].sum(axis=1))-1)
            pore.extend(gamma*np.maximum(level-q['position'][:, 2], 0))
        normal = np.einsum('qa,qa->q', jump, system.normal)
        qn = np.einsum('qa,qa->q', traction, system.normal)
        shear = traction-qn[:, None]*system.normal
        energy = qn**2/(2*system.kn)+np.einsum('qa,qa->q', shear, shear)/(2*system.kt)
        reaction = k@u+seed_force+np.einsum('qag,qa,q->g', system.c, traction, system.weight)-load
        return dict(name=name, u=u.reshape(-1, 3).copy(), stress=np.asarray(stress),
                           strain=np.asarray(strain), void=np.asarray(void), pore=np.asarray(pore),
                           plastic_strain=np.zeros_like(stress), jump=jump.copy(), traction=traction.copy(),
                           plastic_slip=plastic.copy(), accumulated=accumulated.copy(), elastic_energy=energy,
                           plastic_dissipation=plastic_work.copy(), numerical_dissipation=numerical.copy(),
                           endpoint_work=endpoint_work.copy(),
                           physical_pressure=interface['pressure'].copy(), water_load=interface['load'].copy(),
                           reaction=reaction.reshape(-1, 3), gap=normal)
    grids = progress_grids or {}
    for progress in grids.get('geostatic', [i/steps for i in range(1, steps+1)]):
        increment(base+(1-progress)*holding, 0., 'geostatic', progress)
    frames.append(field_state('initial_equilibrated', base))
    origin = 0.
    for stage in case['stages']:
        for progress in grids.get(stage['name'], [i/steps for i in range(1, steps+1)]):
            factor = origin+(stage['target_factor']-origin)*progress
            increment(base+factor*boundary, factor, stage['name'], progress)
        origin = stage['target_factor']; frames.append(field_state(stage['name'], base+origin*boundary))
    return dict(frames=frames, path=path, steps=steps, contact_steps=contact_steps,
                volume_order=order, contact_maximum_m=maximum_m, maximum_relative_force_residual=maximum_residual,
                volume_rule=('affine_cut_hydrostatic' if cut_volume else
                             'input_gauss_jacobi' if discrete_volume else 'independent_gauss_legendre'),
                reference_uses_native_results=False)


def serializable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {k: serializable(v) for k, v in value.items()}
    if isinstance(value, list):
        return [serializable(v) for v in value]
    return value


def differences(a, b):
    """Measurements only. Do not fit new pass thresholds to these differences."""
    assert len(a['frames']) == len(b['frames'])
    rows = []
    for x, y in zip(a['frames'], b['frames']):
        assert x['name'] == y['name']
        row = dict(frame=x['name'])
        for key in x.keys()-{'name'}:
            xx = np.asarray(x[key]); yy = np.asarray(y[key]); assert xx.shape == yy.shape
            assert np.isfinite(xx).all() and np.isfinite(yy).all()
            row[key] = dict(maximum_absolute=float(np.max(np.abs(xx-yy), initial=0)),
                            relative_l2=float(np.linalg.norm(xx-yy)/max(np.linalg.norm(yy), 1e-12)))
        rows.append(row)
    return rows


def native_fields(run):
    assert run['complete'], 'incomplete native run'
    assert len(run['frames']) == len(run['case']['stages'])+1, 'missing native frames'
    frames = []
    for f in run['frames']:
        s = f['state']; p = s['material_states']; j = s['joints']['points']; result = f['joint_results']
        assert all(x['equivalent_plastic_strain'] == 0 for x in p), 'bulk plasticity outside reference'
        frames.append(dict(name=f['name'], u=np.asarray(s['displacement_m']),
                           stress=np.asarray([x['effective_stress_compression_positive_kpa'] for x in p]),
                           plastic_strain=np.asarray([x['plastic_strain_compression_positive'] for x in p]),
                           void=np.asarray([x['void_ratio'] for x in p]),
                           jump=np.asarray([x['jump_m'] for x in j]), traction=np.asarray([x['traction_kpa'] for x in j]),
                           plastic_slip=np.asarray([x['plastic_slip_m'] for x in j]),
                           accumulated=np.asarray([x['accumulated_plastic_slip_m'] for x in j]),
                           plastic_dissipation=np.asarray([x['plastic_dissipation_kpa_m'] for x in j]),
                           numerical_dissipation=np.asarray([x['numerical_dissipation_kpa_m'] for x in j]),
                           endpoint_work=np.asarray([x['endpoint_work_kpa_m'] for x in j]),
                           gap=np.asarray([x['gap_m'] for x in result['points']]),
                           elastic_energy=np.asarray([x['elastic_energy_kpa_m'] for x in result['points']]),
                           physical_pressure=np.asarray(s['joint_water']['physical_pressure_kpa']),
                           water_load=np.asarray(s['joint_water']['nodal_load_kn']).reshape(-1, 3)))
    return dict(frames=frames)


def compare(run, expected):
    """Full-field measurements; existing dry coupling gates cover u/traction.

    Extra fields remain measurements, not a new acceptance specification.
    Energy/work and status invariants are checked separately by native audits.
    """
    native = native_fields(run)
    subset = copy.deepcopy(expected)
    for f in subset['frames']:
        for key in list(f.keys()-native['frames'][0].keys()):
            del f[key]
    rows = differences(native, subset)
    for row in rows:
        assert row['u']['relative_l2'] < .005, ('existing coupling displacement gate', row)
        assert row['traction']['relative_l2'] < .01, ('existing coupling traction gate', row)
    return rows
