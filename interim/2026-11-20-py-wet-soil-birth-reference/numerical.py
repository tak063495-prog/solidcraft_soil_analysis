"""Independent sparse P1 head solve and Q4 weak-force integration.

P1 uses integrated conductivity, NOT analytic series resistance. Q4 integrals
are split at all coefficient/pressure breaks; no production shape functions.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve


def solve_head(problem, z):
    z = np.asarray(z)
    if not np.all(np.isfinite(z)) or np.any(np.diff(z) <= 1e-12*(z[-1]-z[0])):
        raise ValueError('degenerate/non-increasing mesh; near-duplicate layer nodes')
    if z[0] != problem.edges[0] or z[-1] != problem.edges[-1]:
        raise ValueError('mesh endpoints must match prescribed boundary elevations')
    matrix = lil_matrix((len(z), len(z)))
    conductance = []
    for i, (a, b) in enumerate(zip(z[:-1], z[1:])):
        integral_k = sum(max(0., min(b, d)-max(a, c))*k
                         for c, d, k in zip(problem.edges[:-1], problem.edges[1:], problem.conductivity))
        c = integral_k/(b-a)**2
        conductance.append(c)
        matrix[i:i+2, i:i+2] += c*np.array([[1., -1.], [-1., 1.]])
    matrix = matrix.tocsr()
    head = np.zeros(len(z))
    head[[0, -1]] = [problem.bottom_head, problem.top_head]
    head[1:-1] = spsolve(matrix[1:-1, 1:-1], -(matrix @ head)[1:-1])
    residual = matrix @ head
    # Outward prescribed-head flux = -(K H)_boundary, as in DARCY.ja.md.
    out = -residual[[0, -1]]
    flux = -np.asarray(conductance)*np.diff(head)
    return dict(z=z, head=head, flux=flux, outward=out,
                water_defect=float(sum(out)), free_residual=float(max(abs(residual[1:-1]), default=0.)))


def fitted_mesh(problem, n=8):
    return np.r_[np.concatenate([np.linspace(a, b, max(1, round(n*(b-a)))+1)[:-1]
                                for a, b in zip(problem.edges[:-1], problem.edges[1:])]), problem.edges[-1]]


def head_errors(problem, result):
    """Continuous L2 incl. layer kink; nodes alone could hide interpolation error."""
    squared = 0.
    maximum = 0.
    g, w = leggauss(4)
    raw_flux_max = 0.
    for i, (a, b) in enumerate(zip(result['z'][:-1], result['z'][1:])):
        breaks = sorted({a, b, *(x for x in problem.edges if a < x < b)})
        for lo, hi in zip(breaks[:-1], breaks[1:]):
            zz = lo+(g+1)*(hi-lo)/2
            interpolated = result['head'][i]+(zz-a)/(b-a)*(result['head'][i+1]-result['head'][i])
            error = interpolated-problem.head(zz)
            squared += np.dot(w, error*error)*(hi-lo)/2
            raw_flux = -problem.k(zz)*(result['head'][i+1]-result['head'][i])/(b-a)
            raw_flux_max = max(raw_flux_max, float(max(abs(raw_flux-problem.flux))))
        probes = np.asarray(breaks)
        interp = np.interp(probes, result['z'], result['head'])
        maximum = max(maximum, float(max(abs(interp-problem.head(probes)))))
    return dict(head_l2=float(np.sqrt(squared)), head_max=maximum,
                flux_max=float(max(abs(result['flux']-problem.flux))), raw_material_flux_max=raw_flux_max)


def series(problem, ns=(8, 16, 32, 64), column=None):
    rows = []
    for fitted in (True, False):
        previous = None
        for n in ns:
            # One cut element centered on z=1, width 2/n, fixed cut fraction 1/2.
            # Remaining intervals refined independently on each side.
            if fitted:
                z = np.r_[np.linspace(0, 1, n+1), np.linspace(1, 3, 2*n+1)[1:]]
            else:
                z = np.r_[np.linspace(0, 1-1/n, n), np.linspace(1+1/n, 3, 2*n)]
            result = solve_head(problem, z)
            errors = head_errors(problem, result)
            h = float(max(np.diff(z)))
            order = None
            if previous and min(errors['head_l2'], previous['head_l2']) > 1e-11:
                order = float(np.log(previous['head_l2']/errors['head_l2'])/np.log(previous['h']/h))
            row = dict(fitted=fitted, n=n, cells=len(z)-1, h=h, order=order, **errors,
                       water_defect=result['water_defect'], free_residual=result['free_residual'],
                       bottom_out=float(result['outward'][0]), top_out=float(result['outward'][1]),
                       order_note='roundoff/exact P1 field: no inferred order' if fitted else 'unfitted coefficient jump')
            rows.append(row)
            if column is not None:
                m = integrate_forces(column, problem, numerical_result=result)['metrics']
                row['numerical_force_metrics'] = m
            previous = row
    return rows


def shape(x, z, a, b):
    t = (z-a)/(b-a)
    n = np.array([(1-x)*(1-t), x*(1-t), x*t, (1-x)*t])
    grad = np.array([[-(1-t), -(1-x)/(b-a)], [1-t, -x/(b-a)],
                     [t, x/(b-a)], [-t, (1-x)/(b-a)]])
    return n, grad


def integrate_forces(column, flow=None, subdivisions=4, quadrature=3, numerical_result=None):
    """Unit x width and out-of-plane thickness; full x/z nodal force vectors.

    Distinct mechanical interface nodes are retained, with explicit physical
    pairing only for the resultant audit. No contact law or history is solved.
    """
    g, weights = leggauss(quadrature)
    nodes, elements, bodies = [], [], []
    for body, lo, hi in [('old', 0., column.interface), ('new', column.interface, column.top)]:
        offset = len(nodes)
        zz = np.linspace(lo, hi, subdivisions+1)
        nodes.extend([[x, float(z)] for z in zz for x in [0., 1.]])
        for j, (a, b) in enumerate(zip(zz[:-1], zz[1:])):
            elements.append(([offset+2*j, offset+2*j+1, offset+2*j+3, offset+2*j+2], a, b, body))
        bodies.extend([body]*(2*len(zz)))
    nodes = np.array(nodes)
    arrays = {key: np.zeros_like(nodes) for key in
              ['pore_internal', 'effective_internal', 'body', 'pressure_exterior',
               'pressure_interface', 'effective_boundary', 'gradient_volume', 'seepage']}
    breaks = [column.water] if flow is None else list(flow.edges)
    if numerical_result is not None:
        breaks += list(numerical_result['z'])

    def fields(z):
        p = column.pressure(z) if flow is None else flow.pressure(z)
        f = column.fields(z, p)
        dp = (-column.gamma_w if z < column.water else 0.) if flow is None else float(flow.pressure_gradient(z))
        seep = 0. if flow is None else float(flow.seepage(z))
        if numerical_result is not None:
            zz, hh = numerical_result['z'], numerical_result['head']
            i = int(np.clip(np.searchsorted(zz, z, side='right')-1, 0, len(zz)-2))
            dh = (hh[i+1]-hh[i])/(zz[i+1]-zz[i])
            p = column.gamma_w*(np.interp(z, zz, hh)-z)
            if p < 0:
                raise ValueError('negative numerical saturated pressure')
            f = column.fields(z, p)
            dp, seep = column.gamma_w*(dh-1), -column.gamma_w*dh
        return f, dp, seep

    for ids, a, b, body in elements:
        cuts = sorted({a, b, *(v for v in breaks if a < v < b)})
        for lo, hi in zip(cuts[:-1], cuts[1:]):
            for ix, wx in zip(g, weights):
                x = (ix+1)/2
                for iz, wz in zip(g, weights):
                    z = lo+(iz+1)*(hi-lo)/2
                    dv = wx*wz*(hi-lo)/4
                    n, grad = shape(x, z, a, b)
                    f, dp, seep = fields(z)
                    arrays['pore_internal'][ids] -= grad*float(f['pressure'])*dv
                    arrays['effective_internal'][ids] -= grad*np.array([
                        float(f['effective_horizontal']), float(f['effective_vertical'])])*dv
                    arrays['body'][ids, 1] -= n*float(column.weight(z))*dv
                    arrays['gradient_volume'][ids, 1] += n*dp*dv
                    arrays['seepage'][ids, 1] += n*seep*dv
        # All vertical exterior faces; top/bottom only at body boundary.
        edges = [('left', a, b), ('right', a, b)]
        if a == (0. if body == 'old' else column.interface):
            edges.append(('bottom', a, a))
        if b == (column.interface if body == 'old' else column.top):
            edges.append(('top', b, b))
        for side, lo, hi in edges:
            vertical = side in ('left', 'right')
            cuts = sorted({lo, hi, *(v for v in breaks if lo < v < hi)}) if vertical else [0., 1.]
            normal = {'left': [-1., 0.], 'right': [1., 0.],
                      'bottom': [0., -1.], 'top': [0., 1.]}[side]
            normal = np.array(normal)
            interface = not vertical and lo == column.interface
            for c, d in zip(cuts[:-1], cuts[1:]):
                for gg, ww in zip(g, weights):
                    v = c+(gg+1)*(d-c)/2
                    x, z = ((0. if side == 'left' else 1.), v) if vertical else (v, lo)
                    n, _ = shape(x, z, a, b)
                    f, _, _ = fields(z)
                    ds = ww*(d-c)/2
                    key = 'pressure_interface' if interface else 'pressure_exterior'
                    arrays[key][ids] -= np.outer(n, normal)*float(f['pressure'])*ds
                    arrays['effective_boundary'][ids] -= np.outer(n, normal*np.array([
                        float(f['effective_horizontal']), float(f['effective_vertical'])]))*ds
    pressure_boundary = arrays['pressure_exterior']+arrays['pressure_interface']
    identity = arrays['pore_internal']-pressure_boundary-arrays['gradient_volume']
    residual = (arrays['effective_internal']+arrays['pore_internal']-arrays['body']
                -pressure_boundary-arrays['effective_boundary'])
    external = arrays['body']+pressure_boundary+arrays['effective_boundary']
    def moment(force):
        return float(np.sum(nodes[:, 0]*force[:, 1]-nodes[:, 1]*force[:, 0]))
    pair = arrays['pressure_interface'][:2*(subdivisions+1)][-2:]+arrays['pressure_interface'][2*(subdivisions+1):][:2]
    # Before birth: old top is exterior. After: old top/new bottom are paired
    # interface tractions, never also kept in exterior pressure arrays.
    old_ids = np.flatnonzero(np.array(bodies) == 'old')[-2:]
    new_ids = np.flatnonzero(np.array(bodies) == 'new')[:2]
    old_top = arrays['pressure_interface'][old_ids]
    pre_exposed = np.array([[0., -float(fields(column.interface)[0]['pressure'])/2]]*2)
    metrics = dict(pore_identity_max=float(np.max(abs(identity))),
                   nodal_equilibrium_max=float(np.max(abs(residual))),
                   global_force=external.sum(axis=0).tolist(), global_moment=moment(external),
                   continuous_pair_max=float(np.max(abs(pair))),
                   before_old_exposed_top=pre_exposed.tolist(),
                   exposed_to_internalized_max=float(np.max(abs(pre_exposed-old_top))),
                   after_old_interface=old_top.tolist(),
                   after_new_interface=arrays['pressure_interface'][new_ids].tolist(),
                   pressure_boundary_resultant=pressure_boundary.sum(axis=0).tolist(),
                   seepage_resultant=arrays['seepage'].sum(axis=0).tolist(),
                   water_defect=0. if flow is None else float(-flow.flux+flow.flux),
                   water_note='analytic constant flux; numerical boundary flux is reported separately',
                   topology_note='same explicit pressure trace before/after for topology audit; not a pre-BIRTH hydraulic solve')
    return dict(nodes=nodes, bodies=bodies, forces=arrays, metrics=metrics)
