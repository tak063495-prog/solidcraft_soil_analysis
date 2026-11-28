"""Explicitly run 18-Tet4/48-DOF reassembly; never solve a global path.

Reads SHA-verified saved states. All sums/matrix products use fixed index order.
Numbers read from saved JSON are exact stored decimal values. A separate
binary64-lift control identifies input representation effects.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
from mpmath import mp
from a01hp_contact import point

HERE = Path(__file__).resolve().parent


def read(path, binary=False):
    return json.loads(path.read_text(encoding='utf-8'),
                      parse_float=(lambda s: mp.mpf(float(s))) if binary else mp.mpf)


def write(path, value):
    def encode(x):
        if isinstance(x, mp.matrix):
            return [[encode(x[i, j]) for j in range(x.cols)] for i in range(x.rows)]
        if isinstance(x, mp.mpf):
            return mp.nstr(x, mp.dps, strip_zeros=False)
        if isinstance(x, dict):
            return {k: encode(v) for k, v in x.items()}
        if isinstance(x, (tuple, list)):
            return [encode(v) for v in x]
        return x
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8', newline='\n') as f:
        json.dump(encode(value), f, indent=2, ensure_ascii=False)
        f.write('\n')


def norm(a):
    return mp.sqrt(mp.fsum(a[i, j]**2 for i in range(a.rows) for j in range(a.cols)))


def bulk(mesh, ground, transform, committed, displacement):
    nd = len(mesh['nodes'])*3
    full, force = mp.zeros(nd), mp.zeros(nd, 1)
    old_u = mp.matrix([x for row in committed['displacement_m'] for x in row])
    u = mp.matrix([x for row in displacement for x in row])
    materials = {m['id']: m['mohr_coulomb'] for m in ground['materials']}
    volumes = []
    for e, cell in enumerate(mesh['tetrahedra']):
        x = mp.matrix([[1, *mesh['nodes'][n]] for n in cell])
        inv = x**-1
        volume = abs(mp.det(x))/6
        volumes.append(volume)
        b = mp.zeros(6, 12)
        for a in range(4):
            dx, dy, dz = (inv[i, a] for i in (1, 2, 3))
            b[0, 3*a], b[1, 3*a+1], b[2, 3*a+2] = dx, dy, dz
            b[3, 3*a], b[3, 3*a+1] = dy, dx
            b[4, 3*a+1], b[4, 3*a+2] = dz, dy
            b[5, 3*a], b[5, 3*a+2] = dz, dx
        material = materials[mesh['material_ids'][e]]
        young, nu = material['young_kpa'], material['poisson']
        lam, shear = young*nu/((1+nu)*(1-2*nu)), young/(2*(1+nu))
        d = mp.zeros(6)
        for i in range(3):
            for j in range(3):
                d[i, j] = lam + (2*shear if i == j else 0)
        for i in range(3, 6):
            d[i, i] = shear
        ids = [3*n+a for n in cell for a in range(3)]
        k = b.T*d*b*volume
        du = mp.matrix([u[i]-old_u[i] for i in ids])
        sig = mp.matrix(committed['material_states'][e]['effective_stress_compression_positive_kpa']) - d*b*du
        f = -b.T*sig*volume
        for i, ii in enumerate(ids):
            force[ii] += f[i]
            for j, jj in enumerate(ids):
                full[ii, jj] += k[i, j]
    return transform.T*full*transform, transform.T*force, volumes


def assemble(cr, tangents, weights, indices):
    a = mp.zeros(cr[0].cols)
    for q in indices:
        a += cr[q].T*tangents[q]*cr[q]*weights[q]
    return a


def storage_matrix(event):
    s = event['storage']
    if s['kind'] == 'dense':
        return mp.matrix(s['matrix_columns']).T
    a = mp.zeros(len(event['rhs']))
    for i in range(a.rows):
        for j in range(s['row_ptr'][i], s['row_ptr'][i+1]):
            a[i, s['columns'][j]] = s['values'][j]
    return a


def rigid(mesh, transform):
    # Central face-connected component; never select coincident outer nodes.
    cells = [c for c in mesh['tetrahedra']
             if mp.mpf(1)/3 < sum(mesh['nodes'][n][2] for n in c)/4 < mp.mpf(2)/3]
    nodes = sorted({n for c in cells for n in c})
    center = [mp.fsum(mesh['nodes'][n][a] for n in nodes)/len(nodes) for a in range(3)]
    rb = mp.zeros(transform.rows, 6)
    for n in nodes:
        x, y, z = [mesh['nodes'][n][a]-center[a] for a in range(3)]
        for a in range(3):
            rb[3*n+a, a] = 1
        for a, row in enumerate(((0, z, -y), (-z, 0, x), (y, -x, 0))):
            for j, v in enumerate(row):
                rb[3*n+a, j+3] = v
    rb = transform.T*rb
    # Fixed-order modified Gram-Schmidt, with a second orthogonalization pass.
    q = mp.zeros(rb.rows, 6)
    for j in range(6):
        v = rb[:, j]
        for _ in range(2):
            for k in range(j):
                v -= q[:, k]*(q[:, k].T*v)[0]
        q[:, j] = v/norm(v)
    return q, nodes


def analyze(a, vectors, rb, originals):
    u, s, v = mp.svd_r(a)
    cutoff = a.rows*mp.eps*s[0]
    rank = sum(value > cutoff for value in s)
    # Deterministic sign; degeneracies are compared by projectors, not by signs.
    for k in range(a.rows):
        pivot = max(range(a.rows), key=lambda i: abs(v[k, i]))
        if v[k, pivot] < 0:
            v[k, :] *= -1
            u[:, k] *= -1
    modes = []
    for k in (a.rows-2, a.rows-1):
        right, left = v[k, :].T, u[:, k]
        modes.append(dict(index_descending=k, singular_value=s[k], right=right, left=left,
                          middle_rigid_overlap_squared=norm(rb.T*right)**2,
                          right_action_norm=norm(a*right), left_action_norm=norm(a.T*left),
                          projections={name: dict(left_component=(left.T*x)[0],
                              right_component=(right.T*x)[0]) for name, x in vectors.items()}))
    lp = u[:, a.rows-2:a.rows] * u[:, a.rows-2:a.rows].T
    rp = v[a.rows-2:a.rows, :].T * v[a.rows-2:a.rows, :]
    return dict(spectrum_descending=list(s), observation_cutoff=cutoff, numerical_rank=rank,
        cutoff_is_acceptance_tolerance=False, two_smallest=modes,
        left_projector_two=lp, right_projector_two=rp,
        nullity_observed=a.rows-rank,
        projections_two={name: dict(left_norm=norm(lp*x), right_norm=norm(rp*x),
                                   left_projected_vector=lp*x, right_projected_vector=rp*x)
                         for name, x in vectors.items()},
        projections_true_numerical_null={name: norm(u[:, rank:].T*x) for name, x in vectors.items()},
        svd_reconstruction_error=norm(u*mp.diag(s)*v-a),
        differences={name: dict(frobenius=norm(a-b), relative=norm(a-b)/norm(b))
                     for name, b in originals.items()})


def run(source, output, dps, binary=False, frozen_local=False):
    mp.dps = dps
    started = time.monotonic()
    folder = source/'formal/stacked__native__Pnative__Rnative__g8'
    mesh, ground = (read(folder/(name+'.json'), binary) for name in ('mesh', 'ground'))
    kernel = read(source/'separation/native-kernel-output.json', binary)
    events = [json.loads(line, parse_float=(lambda s: mp.mpf(float(s))) if binary else mp.mpf)
              for line in (folder/'native/progress.jsonl').read_text(encoding='utf-8').splitlines()]
    it = [e for e in events if e['kind'] == 'iterate'][-1]
    tangent = [e for e in events if e['kind'] == 'tangent'][-1]
    linear = [e for e in events if e['kind'] == 'a01_linear_system'][-1]
    t = mp.zeros(len(tangent['dof_rows']), len(linear['rhs']))
    for i, row in enumerate(tangent['dof_rows']):
        for j, c in row:
            t[i, j] = c
    assert len(mesh['tetrahedra']) == 18 and t.cols == 48
    a0 = mp.matrix(tangent['matrix_columns']).T
    apost = storage_matrix(linear)
    k, bulk_force, volumes = bulk(mesh, ground, t, it['committed'], it['displacement_m'])
    cr, weights, ad, fd, qvalues, branches = [], [], [], {v: [] for v in (1, 2, 4)}, [], []
    fd_branch_paths = {v: [] for v in (1, 2, 4)}
    target_errors = []
    for q, loc in enumerate(kernel['local_inputs']):
        c = mp.zeros(3, t.rows)
        layout = kernel['layout'][q]
        if frozen_local:
            coefficients = kernel['jump_coefficients'][q]
            weight = layout['weight_m2']
        else:
            # The existing Tri3 three-point rule: permutations of (2/3,1/6,1/6).
            # Re-evaluate its algebra and the face area, without changing the rule.
            shape = [mp.mpf(2)/3 if x > mp.mpf('0.5') else mp.mpf(1)/6 for x in layout['shape']]
            coefficients = [(n, sign*s) for key, sign in (('nodes_a', -1), ('nodes_b', 1))
                            for n, s in zip(layout[key], shape)]
            vertices = [mp.matrix(mesh['nodes'][n]) for n in layout['nodes_a']]
            a, b = vertices[1]-vertices[0], vertices[2]-vertices[0]
            cross = mp.matrix([a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]])
            weight = norm(cross)/6
        for node, coefficient in coefficients:
            for axis in range(3):
                c[axis, 3*node+axis] = coefficient
        cr.append(c*t)
        weights.append(weight)
        material, old = loc['material'], loc['old']
        rebuilt_target = c*mp.matrix([x for row in it['displacement_m'] for x in row])
        target_errors.append(norm(rebuilt_target-mp.matrix(loc['target'])))
        target = mp.matrix(loc['target']) if frozen_local else rebuilt_target
        args = (old['normal'], material['normal_stiffness_kpa_m'],
                material['tangential_stiffness_kpa_m'], material['friction_coefficient'],
                old['jump_m'], old['traction_kpa'], old['birth']['elastic_offset_m'])
        value, derivative, branch = point(*args, list(target), loc['count'])
        qvalues.append(value)
        ad.append(derivative)
        branches.append(branch)
        for divisor in fd:
            step = loc['h']/divisor
            d = mp.zeros(3)
            probes = []
            for axis in range(3):
                e = mp.eye(3)[:, axis]*step
                plus, _, pb = point(*args, list(target+e), loc['count'])
                minus, _, mb = point(*args, list(target-e), loc['count'])
                probes.append(dict(axis=axis, plus=pb, minus=mb))
                d[:, axis] = (plus-minus)/(2*step)
            fd[divisor].append(d)
            fd_branch_paths[divisor].append(probes)
    indices = {b: [i for i, row in enumerate(branches) if row['branch'] == b]
               for b in ('open', 'stick', 'slip')}
    contact_force = mp.zeros(t.cols, 1)
    for q in range(len(cr)):
        contact_force += cr[q].T*qvalues[q]*weights[q]
    external = t.T*mp.matrix(it['target'])
    native_run = read(folder/'native/run.json', binary)
    holding = t.T*mp.matrix(native_run['transition_audits'][0]['initial_holding_load_kn'])
    # Production reduced_residual has the RHS sign: external minus internal.
    residual = external - bulk_force - contact_force
    vectors = dict(saved_residual=mp.matrix(it['reduced_residual']), reassembled_residual=residual,
                   external_load=external, birth_initial_holding=holding,
                   birth_holding_at_stop=mp.zeros(t.cols, 1),
                   internal_bulk=bulk_force, internal_contact=contact_force,
                   linear_rhs_after_predictor=mp.matrix(linear['rhs']))
    rb, middle = rigid(mesh, t)
    matrices, records, local = {}, {}, {}
    for name, derivatives in [('analytic', ad)]+[('fd_h_div_'+str(v), fd[v]) for v in (1, 2, 4)]:
        parts = {b: assemble(cr, derivatives, weights, indices[b]) for b in indices}
        total = k + sum(parts.values(), mp.zeros(t.cols))
        matrices[name] = total
        originals = dict(saved_pre_predictor=a0, saved_solver_matrix=apost)
        if name != 'analytic':
            originals['hp_analytic'] = matrices['analytic']
        records[name] = analyze(total, vectors, rb, originals)
        records[name]['contact_rigid_spectrum'] = list(mp.svd_r((total-k)*rb, compute_uv=False))
        local[name] = dict(local_tangents=derivatives, branch_matrices=parts, bulk_matrix=k, total_matrix=total)
        print(dps, 'binary' if binary else 'decimal', name,
              [mp.nstr(x, 12) for x in records[name]['spectrum_descending'][-3:]], flush=True)
    for name, original in (('saved_pre_predictor', a0), ('saved_solver_matrix', apost)):
        records[name] = analyze(original, vectors, rb, dict(hp_analytic=matrices['analytic']))
    write(output, dict(dps=dps, input_representation='exact_binary64' if binary else 'exact_saved_decimal',
        frozen_local_inputs=frozen_local,
        passed=False, a01_complete=False, diagnostic_only=True,
        context={key: it[key] for key in ('stage','progress','iteration','relative_force_residual')},
        element_count=18, reduced_dofs=t.cols, middle_nodes=middle, branches=branches,
        branch_gp_indices=indices, counts=[loc['count'] for loc in kernel['local_inputs']],
        h=[loc['h'] for loc in kernel['local_inputs']],
        target_reassembly_differences=target_errors, volumes=volumes,
        vectors=vectors, residual_reassembly_difference=norm(residual-vectors['saved_residual']),
        saved_predictor_matrix_change=norm(apost-a0),
        new_global_solve=False, birth_holding_active_at_stop=False,
        residual_sign='external_minus_internal', fd_branch_paths=fd_branch_paths,
        records=records, raw=local, elapsed_seconds=time.monotonic()-started))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--dps', type=int, choices=(50, 80), required=True)
    p.add_argument('--binary64-inputs', action='store_true')
    p.add_argument('--frozen-local-inputs', action='store_true')
    a = p.parse_args()
    run(a.source, a.output, a.dps, a.binary64_inputs, a.frozen_local_inputs)
