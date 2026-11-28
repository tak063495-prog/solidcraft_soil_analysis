"""Observe range, weighted projection, order and roundoff; never repair a force."""
import argparse
import copy
from pathlib import Path
import numpy as np
from a02cause_binding import study
from a02cause_metrics import metric


def solve(matrix, rhs, weight, normal, mu, row_order=None, gp_order=None):
    rows = np.arange(matrix.shape[0]) if row_order is None else row_order
    gps = np.arange(len(weight)) if gp_order is None else gp_order
    cols = (3*gps[:, None]+np.arange(3)).ravel()
    permuted = matrix[np.ix_(rows, cols)]
    xp, _, rank, singular = np.linalg.lstsq(permuted, rhs[rows], rcond=1e-12)
    x = np.empty_like(xp); x[cols] = xp
    q = x.reshape(-1, 3)/np.sqrt(weight)[:, None]
    qn = np.einsum('qa,qa->q', q, normal)
    excess = np.linalg.norm(q-qn[:, None]*normal, axis=1)-mu*(-qn)
    error = matrix@x-rhs
    return q, dict(rank=int(rank), singular_values=singular.tolist(),
        nodal_error_kn=error.tolist(), maximum_nodal_error_kn=float(max(abs(error))),
        normal_traction_kpa=qn.tolist(), friction_excess_kpa=excess.tolist(),
        nodal_gate=bool(max(abs(error)) <= 1e-10), cone_gate=bool(np.all(qn <= 0) and np.all(excess <= 0)),
        uncorrected_traction_kpa=q.tolist())


def analyze(out, probe):
    out.mkdir()
    mesh, ground, case = study.fixture('transfer_force', 4)
    s = study.Model(mesh, ground, case).system
    matrix = np.einsum('qag,q->gqa', s.c, np.sqrt(s.weight)).reshape(3*len(mesh['nodes']), -1)
    force = lambda q: np.einsum('qag,a,q->g', s.c, q, s.weight)
    good = force([12., 5., -100.]); unrep = good.copy()
    tip = case['joint_system']['joints'][0]['crack_tip']['edges'][0]
    direction = np.array(mesh['nodes'][tip[1]])-mesh['nodes'][tip[0]]
    unrep[3*tip[0]:3*tip[0]+3] += direction; unrep[3*tip[1]:3*tip[1]+3] -= direction
    study.write(out/'mesh.json', mesh); study.write(out/'case.json', case)
    result = dict(matrix_shape=list(matrix.shape), zero_nodal_rows=np.flatnonzero(np.all(matrix == 0, axis=1)).tolist(),
                  tip_nodes=tip, requested_forces_never_modified=True, local_return_on_invalid_initial_state=False,
                  cases={}, native={})
    for name, rhs in [('good', good), ('tension', force([0., 0., 100.])),
                      ('outside-circle', force([40., 0., -100.])), ('unreproducible-tip', unrep),
                      ('cone-boundary', force([30., 0., -100.]))]:
        original = rhs.copy(); q, baseline = solve(matrix, rhs, s.weight, s.normal, s.mu)
        variants = {}
        for variant in ('repeat', 'reverse-nodes', 'reverse-gps', 'nextafter', 'decimal12', 'decimal10'):
            rr = rhs.copy(); rows = gps = None
            if variant == 'reverse-nodes': rows = np.arange(len(rhs))[::-1]
            if variant == 'reverse-gps': gps = np.arange(len(s.weight))[::-1]
            if variant == 'nextafter': rr = np.nextafter(rr, np.inf)
            if variant.startswith('decimal'): rr = np.round(rr, int(variant[7:]))
            qq, obs = solve(matrix, rr, s.weight, s.normal, s.mu, rows, gps)
            variants[variant] = dict(input_change=metric(rr, rhs), traction_change=metric(qq, q),
                bitwise_same_traction=bool(qq.tobytes() == q.tobytes()), observation=obs)
        assert np.array_equal(original, rhs)
        result['cases'][name] = dict(requested_residual_kn=rhs.tolist(), baseline=baseline, variants=variants)
        result['native'][name] = []
        for variant in ('original', 'repeat', 'reversed-face-declaration'):
            cc = copy.deepcopy(case)
            if variant.startswith('reversed'): cc['joint_system']['joints'][0]['crack_tip']['faces'].reverse()
            prefix = out/(name+'-'+variant)
            study.write(prefix.with_suffix('.case.json'), cc)
            study.write(prefix.with_suffix('.request.json'), dict(operation='birth', residual_kn=rhs.tolist(), inactive_owner=None))
            process = study.invoke([probe, out/'mesh.json', prefix.with_suffix('.case.json'), prefix.with_suffix('.request.json'),
                                    prefix.with_suffix('.result.json')], prefix.with_suffix('.log'))
            response = study.read(prefix.with_suffix('.result.json'))
            assert process['exit_code'] == (0 if response['accepted'] else 1) and not process['timed_out']
            result['native'][name].append(dict(variant=variant, accepted=response['accepted'], error=response['error'], process=process))
        a = study.read((out/(name+'-original')).with_suffix('.result.json'))
        b = study.read((out/(name+'-repeat')).with_suffix('.result.json'))
        result['native'][name][1]['exact_same_except_pid'] = study.exact({k:v for k,v in a.items() if k != 'pid'}) == study.exact({k:v for k,v in b.items() if k != 'pid'})
        assert result['native'][name][1]['exact_same_except_pid']
        if a['accepted']:
            actual = np.array([p['birth']['initial_traction_kpa'] for p in a['result']['state']['points']])
            result['native'][name][0]['against_independent_projection'] = metric(actual, q)
            c = study.read((out/(name+'-reversed-face-declaration')).with_suffix('.result.json'))
            def mapped(r):
                values = r['result']['state']['points']; layout = r['layout']
                return {(p['owner_a'],p['owner_b'],tuple(p['shape'])):values[i]['birth']['initial_traction_kpa'] for i,p in enumerate(layout)}
            aa, bb = mapped(a), mapped(c); assert set(aa) == set(bb)
            result['native'][name][2]['mapped_traction_change'] = metric([aa[k] for k in sorted(aa)], [bb[k] for k in sorted(bb)])
    u, sig, _ = np.linalg.svd(matrix, full_matrices=False)
    rank = int(np.count_nonzero(sig > sig[0]*1e-12))
    result.update(numerical_rank=rank, gp_nullity=int(matrix.shape[1]-rank),
                  cutoff_relative=1e-12, rank_is_observation_not_new_acceptance=True)
    outside = unrep-u[:, :rank]@(u[:, :rank].T@unrep)
    result['unreproducible_range_component_kn'] = outside.tolist()
    result['unreproducible_range_component_maximum_kn'] = float(max(abs(outside)))
    study.write(out/'diagnosis.json', result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--probe', type=Path, required=True); a = p.parse_args()
    analyze(a.output, a.probe.resolve())


if __name__ == '__main__': main()
