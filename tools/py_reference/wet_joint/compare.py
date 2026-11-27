"""Field comparisons on explicit common load coordinates; no fitted pass gate."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

FIELDS = ('u', 'jump', 'traction', 'plastic_slip', 'accumulated', 'gap',
          'physical_pressure', 'water_load')


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    with Path(path).open('x', encoding='utf-8') as f:
        json.dump(value, f, indent=2, ensure_ascii=False, allow_nan=False); f.write('\n')


def save(folder, result):
    folder = Path(folder); folder.mkdir(parents=True, exist_ok=False)
    with (folder/'fields.npz').open('xb') as f:
        np.savez_compressed(f, **result['fields'])
    write(folder/'history.json', dict(rows=result['rows'], metadata=result['metadata'],
                                    fields_sha256=digest(folder/'fields.npz')))


def validate(result):
    meta = result['metadata']; n = meta['steps']; rows = result['rows']
    if meta.get('complete') is not True or meta.get('numerical_acceptance') is not False:
        raise ValueError('incomplete or mislabelled diagnostic')
    names = [r['stage'] for r in rows[::n]]
    if len(rows) != 7*n or len(set(names)) != 7:
        raise ValueError('all seven stages required')
    for phase, name in enumerate(names):
        selected = rows[phase*n:(phase+1)*n]
        if [r['stage'] for r in selected] != [name]*n or [r['progress'] for r in selected] != [i/n for i in range(1, n+1)]:
            raise ValueError('missing, duplicate, or nonuniform accepted schedule')
    for row in rows:
        if not np.isfinite([row['factor'], row['progress'], row['relative_force_residual']]).all():
            raise ValueError('nonfinite accepted row')
        if not 0 <= row['relative_force_residual'] < meta['force_tolerance']:
            raise ValueError('accepted force tolerance violated')
    for name, values in result['fields'].items():
        if len(values) != len(rows) or not np.isfinite(values).all():
            raise ValueError(('nonfinite or incomplete field', name))
    return result


def load(folder):
    folder = Path(folder); manifest = read(folder/'history.json')
    if manifest['fields_sha256'] != digest(folder/'fields.npz'):
        raise ValueError('field archive SHA mismatch')
    with np.load(folder/'fields.npz', allow_pickle=False) as f:
        fields = {k: f[k].copy() for k in f.files}
    return validate(dict(rows=manifest['rows'], metadata=manifest['metadata'], fields=fields))


def coordinates(result):
    validate(result)
    return {(r['stage'], Fraction(r['progress']).limit_denominator(1000000)): i
            for i, r in enumerate(result['rows'])}


def common(results):
    maps = [coordinates(r) for r in results]
    keys = [key for key in maps[0] if all(key in m for m in maps[1:])]
    if not keys:
        raise ValueError('no common accepted coordinates')
    for key in keys:
        factors = [r['rows'][m[key]]['factor'] for r, m in zip(results, maps)]
        if any(factor != factors[0] for factor in factors[1:]):
            raise ValueError('boundary factors differ at common coordinates')
    return keys, [np.array([m[k] for k in keys]) for m in maps]


def metric(a, b):
    if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('incompatible/nonfinite comparison')
    difference = a-b; norm = float(np.linalg.norm(difference)); denominator = float(np.linalg.norm(b))
    return dict(absolute_l2=norm, reference_l2=denominator,
                relative_l2=norm/denominator if denominator else (0. if norm == 0 else None),
                relative_l2_undefined=denominator == 0 and norm != 0,
                maximum_absolute=float(np.max(abs(difference), initial=0.)))


def matching(a, b, *, require_all=True):
    a = np.asarray(a); b = np.asarray(b)
    distance, index = cKDTree(b).query(a)
    selected = np.flatnonzero(distance < 1e-10)
    if require_all and len(selected) != len(a):
        raise ValueError('point coordinates differ; use spatial averages')
    if len(set(index[selected])) != len(selected):
        raise ValueError('ambiguous spatial correspondence')
    return selected, index[selected]


def matching_nodes(a, b):
    left = np.array(a['node_position']); right = np.array(b['node_position'])
    ids_a = np.array(a['node_body']); ids_b = np.array(b['node_body'])
    selected = []; mapped = []
    for body in sorted(set(ids_a)):
        ia = np.flatnonzero(ids_a == body); ib = np.flatnonzero(ids_b == body)
        if not len(ib):
            raise ValueError('missing mechanical body')
        aa, bb = matching(left[ia], right[ib])
        selected.extend(ia[aa]); mapped.extend(ib[bb])
    order = np.argsort(selected)
    return np.array(selected)[order], np.array(mapped)[order]


def pair(a, b):
    keys, (ia, ib) = common([a, b]); out = {}
    _, node = matching_nodes(a['metadata'], b['metadata'])
    _, gp = matching(a['metadata']['contact_position'], b['metadata']['contact_position'])
    for key in FIELDS:
        x = a['fields'][key][ia]; y = b['fields'][key][ib]
        y = y[:, node if key in ('u', 'water_load') else gp]
        out[key] = metric(x, y)
    return dict(common_coordinates=len(keys), fields=out, numerical_acceptance=False)


def observed(d1, d2, ratio=2.):
    n1, n2 = float(np.linalg.norm(d1)), float(np.linalg.norm(d2))
    if ratio <= 1:
        raise ValueError('refinement ratio must exceed one')
    return dict(order=float(np.log(n1/n2)/np.log(ratio)) if n1 > 0 and n2 > 0 else None,
                difference_direction_cosine=float(np.sum(d1*d2)/(n1*n2)) if n1 > 0 and n2 > 0 else None,
                coarse_difference_l2=n1, fine_difference_l2=n2,
                asymptotic_convergence_certified=False)


def load_series(results):
    if len(results) < 2:
        raise ValueError('at least two load levels required')
    steps = [r['metadata']['steps'] for r in results]
    controls = [(r['metadata']['local'], r['metadata'].get('native_contact_maximum_m')) for r in results]
    if any(control != controls[0] for control in controls[1:]):
        raise ValueError('load series requires fixed contact integration controls')
    if any(b <= a for a, b in zip(steps, steps[1:])):
        raise ValueError('load series must be strictly refined')
    pairs = [dict(steps=[a, b], **pair(x, y)) for a, b, x, y in zip(steps, steps[1:], results, results[1:])]
    orders = []
    for i in range(len(results)-2):
        group = results[i:i+3]; keys, ids = common(group)
        ratio = steps[i+1]/steps[i]
        if not np.isclose(ratio, steps[i+2]/steps[i+1]):
            raise ValueError('observed order needs a geometric load series')
        fields = {}
        for key in FIELDS:
            values = []
            for r, selected in zip(group, ids):
                axis = 'node_position' if key in ('u', 'water_load') else 'contact_position'
                _, mapping = (matching_nodes(group[0]['metadata'], r['metadata']) if axis == 'node_position'
                              else matching(group[0]['metadata'][axis], r['metadata'][axis]))
                values.append(r['fields'][key][selected][:, mapping])
            fields[key] = observed(values[0]-values[1], values[1]-values[2], ratio)
        orders.append(dict(steps=steps[i:i+3], common_coordinates=len(keys), fields=fields))
    return dict(axis='global_load_increments', pairs=pairs, observed_orders=orders, numerical_acceptance=False)


def patch_means(result, triangles):
    """Area-weighted averages on original parent triangles, never GP identity."""
    positions = np.array(result['metadata']['contact_position']); weights = np.array(result['metadata']['contact_weight'])
    groups = np.full(len(positions), -1)
    for i, triangle in enumerate(np.asarray(triangles)):
        matrix = (triangle[1:]-triangle[0]).T
        uv = np.linalg.lstsq(matrix, (positions-triangle[0]).T, rcond=None)[0].T
        projected = triangle[0]+uv@matrix.T
        inside = ((np.linalg.norm(projected-positions, axis=1) < 1e-10) & (uv.min(axis=1) >= -1e-10)
                  & (uv.sum(axis=1) <= 1+1e-10))
        if np.any(inside & (groups >= 0)):
            raise ValueError('contact GP on ambiguous parent boundary')
        groups[inside] = i
    if np.any(groups < 0):
        raise ValueError('contact point outside parent faces')
    out = {}; areas = []
    for i, triangle in enumerate(np.asarray(triangles)):
        mask = groups == i; area = weights[mask].sum()
        exact = np.linalg.norm(np.cross(triangle[1]-triangle[0], triangle[2]-triangle[0]))/2
        if not np.isclose(area, exact, rtol=1e-11, atol=1e-12):
            raise ValueError('incomplete parent face area')
        areas.append(float(area))
    for key in FIELDS:
        if key in ('u', 'water_load'):
            continue
        values = result['fields'][key]
        out[key] = np.stack([np.einsum('tq...,q->t...', values[:, groups == i], weights[groups == i])/areas[i]
                             for i in range(len(triangles))], axis=1)
    return out, areas


def spatial_series(results, triangles):
    if len(results) < 2:
        raise ValueError('at least two spatial levels required')
    if len({r['metadata']['steps'] for r in results}) != 1 or len({r['metadata']['local'] for r in results}) != 1:
        raise ValueError('space series requires fixed global and contact partitions')
    counts = [r['metadata']['elements'] for r in results]
    if any(b != 8*a for a, b in zip(counts, counts[1:])):
        raise ValueError('space series requires uniform h/2 meshes')
    keys, ids = common(results); means = [patch_means(r, triangles)[0] for r in results]
    pairs = []
    for i in range(len(results)-1):
        na, nb = matching_nodes(results[i]['metadata'], results[i+1]['metadata'])
        pairs.append(dict(elements=counts[i:i+2], common_coordinates=len(keys), common_nodes=len(na),
                          u=metric(results[i]['fields']['u'][ids[i]][:, na], results[i+1]['fields']['u'][ids[i+1]][:, nb]),
                          parent_area_means={key: metric(means[i][key][ids[i]], means[i+1][key][ids[i+1]]) for key in means[i]}))
    orders = []
    for i in range(len(results)-2):
        orders.append(dict(elements=counts[i:i+3], parent_area_means={
            key: observed(means[i][key][ids[i]]-means[i+1][key][ids[i+1]],
                          means[i+1][key][ids[i+1]]-means[i+2][key][ids[i+2]]) for key in means[i]}))
    return dict(axis='spatial_h', pairs=pairs, observed_orders=orders,
                order_unavailable_reason='three spatial levels required' if not orders else None,
                pointwise_spatial_convergence_certified=False, numerical_acceptance=False)


def diagnose(result, model):
    validate(result); f = result['fields']; events = []; stage_error = {}; previous = None
    maximum_residual = 0.; max_reported = 0.
    contact_pressure = -np.einsum('tqa,qa->tq', f['traction'], model.normal)
    pressure_changes = np.diff(contact_pressure, axis=0)
    maximum_pressure_drop = float(np.max(-pressure_changes, initial=0.))
    for i, row in enumerate(result['rows']):
        load = model.base+(1-row['progress'])*(model.seed-model.base) if row['stage'] == 'geostatic' else model.base
        force = model.k@f['u'][i].ravel()+model.seed+np.einsum('qai,qa,q->i', model.c, f['traction'][i], model.weight)
        residual = float(np.linalg.norm(model.t.T@(force-load))/max(np.linalg.norm(load), 1.))
        maximum_residual = max(maximum_residual, residual)
        max_reported = max(max_reported, row['relative_force_residual'])
        if previous is not None:
            delta = f['accumulated'][i]-f['accumulated'][i-1]
            if np.min(delta) < -1e-12:
                raise ValueError('nonmonotone accumulated slip')
            # 1e-14 m is an event display threshold, never an accuracy gate.
            state = np.where(f['gap'][i] > 0, 2, np.where(delta > 1e-14, 1, 0))
            old_state = previous
            for q in np.flatnonzero(state != old_state):
                events.append(dict(stage=row['stage'], progress_interval=[result['rows'][i-1]['progress'] if result['rows'][i-1]['stage'] == row['stage'] else 0., row['progress']],
                                   point=int(q), before=int(old_state[q]), after=int(state[q]),
                                   gap_m=float(f['gap'][i, q]), accumulated_increment_m=float(delta[q])))
            previous = state
        else:
            previous = np.where(f['gap'][i] > 0, 2, 0)
    for key in ('physical_pressure', 'pore', 'water_load'):
        if key in f:
            stage_error[key] = float(np.max(abs(f[key]-f[key][0]), initial=0.))
    jump_error = float(np.max(abs(f['jump']-np.einsum('qai,ti->tqa', model.c, f['u'].reshape(len(result['rows']), -1)))))
    normal_q = np.einsum('tqa,qa->tq', f['traction'], model.normal)
    shear_q = f['traction']-normal_q[..., None]*model.normal
    return dict(accepted_states=len(result['rows']), maximum_reassembled_force_residual=maximum_residual,
                maximum_displacement_jump_error_m=jump_error,
                maximum_normal_law_error_kpa=float(np.max(abs(normal_q-model.kn*np.minimum(f['gap'], 0)))),
                maximum_friction_excess_kpa=float(np.max(np.linalg.norm(shear_q, axis=-1)-model.mu*np.maximum(-normal_q, 0), initial=0.)),
                maximum_physical_pressure_error_kpa=float(np.max(abs(f['physical_pressure']-model.pressure))),
                maximum_water_load_error_kn=float(np.max(abs(f['water_load']-model.water_load.reshape(-1, 3)))),
                maximum_reported_force_residual=max_reported, invariant_water_maximum_changes=stage_error,
                endpoint_and_increment_events=events, event_codes={'0': 'closed/no slip increment', '1': 'closed/slip increment', '2': 'open'},
                event_display_threshold_m=1e-14, maximum_effective_contact_pressure_drop_kpa=maximum_pressure_drop,
                endpoint_events_do_not_resolve_within_increment_path=True, numerical_acceptance=False)


def historical_comparison(current_load_report, path):
    previous = read(path)
    if previous['passed'] is not False or previous['original_16_128_parameter'] != 'global load increments per phase, not mesh subdivisions':
        raise ValueError('historical record identity/acceptance differs')
    rows = []
    for row in current_load_report['pairs']:
        a, b = row['steps']; name = f'{a}__{b}'
        if name not in previous['global_only_cut']:
            continue
        old = previous['global_only_cut'][name]['fields']['plastic_slip']['pooled_relative_l2']
        new = row['fields']['plastic_slip']['relative_l2']
        rows.append(dict(steps=[a, b], recorded_relative_l2=old, new_relative_l2=new, new_minus_recorded=new-old))
    return dict(source_sha256=digest(path), historical_passed=previous['passed'], comparisons=rows,
                recorded_order=previous['richardson']['global_cut']['fields']['plastic_slip']['observed_order'],
                numerical_acceptance=False)


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--axis', choices=('load', 'space', 'pair'), required=True)
    p.add_argument('--paths', nargs='+', type=Path, required=True); p.add_argument('--parents', type=Path)
    p.add_argument('--output', type=Path, required=True); p.add_argument('--historical', type=Path); a = p.parse_args()
    values = [load(path) for path in a.paths]
    if a.axis == 'pair':
        if len(values) != 2:
            p.error('pair needs exactly two paths')
        value = pair(*values)
    elif a.axis == 'space':
        if a.parents is None:
            p.error('--parents JSON triangles is required for space')
        value = spatial_series(values, read(a.parents))
    else:
        value = load_series(values)
    if a.historical:
        if a.axis != 'load':
            p.error('--historical applies to the load series only')
        value['historical_comparison'] = historical_comparison(value, a.historical)
    write(a.output, value)


if __name__ == '__main__':
    main()
