"""Read fresh full-equilibrium evidence; never synthesize accepted solutions."""
import argparse
import itertools
import json
import math
from pathlib import Path
import numpy as np
from a03follow_binding import study, fields, reference
from a03follow_controls import canonical, compare, common, metric, norm, dot, controls
from a03follow_run import write

SELECTED = ('u', 'stress', 'pore_pressure', 'reaction', 'traction', 'endpoint_work',
            'numerical_dissipation', 'load_work_cumulative', 'elastic_energy')


def load(run, spec, engine, common_n=None):
    folder = run/spec['name']
    mapping = study.read(folder/(engine+'-layout.json'))
    rows = [(r, canonical(v, mapping)) for r, v in fields.load_fields(folder/(engine+'-fields.jsonl'))]
    return common(rows, common_n) if common_n is not None else rows


def values(rows, key):
    return np.stack([f[key] for _, f in rows])


def order(x, y, z):
    d1, d2 = np.asarray(y)-x, np.asarray(z)-y
    n1, n2 = norm(d1), norm(d2)
    return dict(coarse_middle_absolute_l2=n1, middle_fine_absolute_l2=n2,
                observed_order=math.log2(n1/n2) if n1 and n2 else None,
                difference_cosine=dot(d1, d2)/(n1*n2) if n1 and n2 else None,
                asymptotic_certification=False, extrapolation_used_for_acceptance=False)


def physical_fingerprint(folder):
    mesh, ground, case = (study.read(p) for p in study.paths(folder))
    elements = sorted(zip(mesh['tetrahedra'], mesh['tetrahedra10'], mesh['material_ids'], mesh['region_ids']))
    for i, key in enumerate(('tetrahedra', 'tetrahedra10', 'material_ids', 'region_ids')):
        mesh[key] = [e[i] for e in elements]
    mesh['face_sets'] = {k: sorted(v) for k, v in mesh['face_sets'].items()}
    import hashlib
    return hashlib.sha256(json.dumps([mesh, ground, case], sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def raw_observations(run, spec, engine):
    folder = run/spec['name']
    rows = load(run, spec, engine)
    with controls(spec['volume'], spec['reverse'], spec['minimum']):
        model = reference.Model(*[study.read(p) for p in study.paths(folder)])
    mapping = study.read(folder/(engine+'-layout.json'))
    weights = model.system.weight[mapping['contact_order']]
    if engine == 'native':
        actual = study.read(folder/'native/run.json')['joint_integration_points']
        weights = np.asarray([q['weight_m2'] for q in actual])[mapping['contact_order']]
        residuals = [a['relative_force_residual'] for a in study.read(folder/'native-audits.json')['rows']]
    else:
        residuals = [a['independent_relative_force_residual'] for a in study.read(folder/'independent-result.json')['rows']]
    result = []
    for (row, f), residual in zip(rows, residuals):
        original_active = np.empty_like(f['active_elements'])
        original_active[mapping['element_order']] = f['active_elements']
        transform, _ = reference.constraint_space(model.mesh, model.case, original_active, model.system, bool(f['joint_tied'][0]))
        projected = transform.T@(f['internal']-f['applied_load']).ravel(order='C')
        fixed_residual = norm(projected)/max(norm(f['applied_load']), 1.)
        integrated = {k+'_kn_m': dot(weights, f[k]) for k in ('endpoint_work', 'numerical_dissipation', 'plastic_dissipation', 'elastic_energy')}
        result.append(dict(row=row, **integrated,
            load_work_cumulative_kn_m=float(f['load_work_cumulative']),
            reaction_work_cumulative_kn_m=float(f['reaction_work_cumulative']),
            free_relative_equilibrium_residual=fixed_residual,
            legacy_numpy_relative_equilibrium_residual=residual,
            raw_residual_including_supports_l2_kn=norm(f['reaction']),
            raw_residual_resultant_kn=[math.fsum(map(float, f['reaction'][:, i])) for i in range(3)],
            maximum_local_work_identity_defect_kpa_m=float(np.max(np.abs(f['initial_energy']+f['endpoint_work']-f['elastic_energy']-f['plastic_dissipation']-f['numerical_dissipation']))),
            physical_plastic_work_kpa=float(np.max(np.abs(f['plastic_work'])))))
    assert len(result) == 4*spec['n']
    return result


def assess(run, output, vector_output=None):
    plan = study.read(run/'plan.json')['specs']
    driver = study.read(run/'driver-result.json')
    result = dict(passed=False, a03_complete=False, diagnostic_only=True, all_field_numerical_acceptance=False,
                  finest_grid_is_exact_reference=False, paths=driver['paths'], geometry={}, ordering={}, cross={}, raw={},
                  water_accounting=dict(flux=None, storage=None, birth_water_transfer=None, balance_acceptance=None),
                  formal_paths=driver['formal_paths'], formal_states=driver['formal_states'],
                  compatibility_paths=driver['compatibility_paths'], compatibility_states=driver['compatibility_states'],
                  standard_cli_paths=driver['standard_cli_paths'], standard_cli_states=driver['standard_cli_states'])
    vectors = {}
    by_name = {s['name']: s for s in plan}
    def find(**kwargs):
        found = [s for s in plan if all(s[k] == v for k, v in kwargs.items())]
        if len(found) != 1:
            raise ValueError(('nonunique specification', kwargs))
        return found[0]
    def difference(name, a, b):
        comparison = compare(a, b)
        vectors[name] = {key: (values(a, key)-values(b, key)).tolist() for key in SELECTED}
        return comparison
    for spec in plan:
        for engine in spec['engines']:
            result['raw'][spec['name']+'/'+engine] = raw_observations(run, spec, engine)
    for shape, stress in itertools.product(('wide', 'tall'), ('k05', 'k08')):
        specs = [find(shape=shape, stress=stress, n=n, minimum=1, volume=5, element=False, face=False, reverse=False) for n in (4, 8, 16)]
        gl = find(shape=shape, stress=stress, n=8, volume=7)
        label = shape+'-'+stress
        entry = {}
        for engine in ('native', 'independent'):
            rows = [load(run, s, engine, 4) for s in specs]
            entry[engine] = dict(
                n4_to_8=difference(label+'/'+engine+'/4-8', rows[0], rows[1]),
                n8_to_16=difference(label+'/'+engine+'/8-16', rows[1], rows[2]),
                orders={key: order(*(values(r, key) for r in rows)) for key in SELECTED})
        entry['gl5_gl7'] = difference(label+'/gl5-gl7', load(run, specs[1], 'independent'), load(run, gl, 'independent'))
        entry['implementation'] = {str(s['n']): difference(label+'/implementation/'+str(s['n']), load(run, s, 'native'), load(run, s, 'independent')) for s in specs}
        result['geometry'][label] = entry
    order_specs = [s for s in plan if 'order' in s['roles']]
    baseline = find(shape='wide', stress='k05', n=8, minimum=1, volume=5, element=False, face=False, reverse=False)
    fingerprints = {s['name']: physical_fingerprint(run/s['name']) for s in order_specs}
    assert len(set(fingerprints.values())) == 1
    result['ordering']['same_physical_fingerprint'] = fingerprints
    result['ordering']['aggregation'] = 'physical connectivity / exact local rule shape / row / component; math.fsum'
    for engine in ('native', 'independent'):
        base_rows = load(run, baseline, engine)
        base_layout = study.read(run/baseline['name']/(engine+'-layout.json'))
        entry = dict(all_combinations={}, independent_effects={}, same_exact_physical_keys=True)
        for s in order_specs:
            mapping = study.read(run/s['name']/(engine+'-layout.json'))
            assert mapping['element_keys'] == base_layout['element_keys']
            assert mapping['contact_keys'] == base_layout['contact_keys']
            entry['all_combinations'][s['name']] = difference('ordering/'+engine+'/'+s['name'], load(run, s, engine), base_rows)
        for factor in ('element', 'face', 'reverse'):
            pairs = []
            for s in order_specs:
                if s[factor]:
                    match = next(t for t in order_specs if all(t[k] == (False if k == factor else s[k]) for k in ('element', 'face', 'reverse')))
                    pairs.append(dict(changed=s['name'], baseline=match['name'],
                        effect=difference('ordering/'+engine+'/'+factor+'/'+s['name'], load(run, s, engine), load(run, match, engine))))
            entry['independent_effects'][factor] = pairs
        result['ordering'][engine] = entry
    for engine in ('native', 'independent'):
        table = {(s['n'], s['minimum']): load(run, s, engine, 4) for s in plan if 'cross' in s['roles']}
        entry = dict(common_coordinates_per_stage=4, total_common_coordinates=16, local={}, global_step={}, interactions={}, orders={})
        for n in (4, 8, 16, 32):
            for a, b in ((1, 2), (2, 4)):
                name = f'N{n}-M{a}-{b}'
                entry['local'][name] = difference('cross/'+engine+'/local/'+name, table[n, b], table[n, a])
        for m in (1, 2, 4):
            for a, b in ((4, 8), (8, 16), (16, 32)):
                name = f'M{m}-N{a}-{b}'
                entry['global_step'][name] = difference('cross/'+engine+'/global/'+name, table[b, m], table[a, m])
            for ns in ((4, 8, 16), (8, 16, 32)):
                entry['orders'][f'M{m}-N{ns[0]}-{ns[1]}-{ns[2]}'] = {
                    key: order(*(values(table[n, m], key) for n in ns)) for key in SELECTED}
        for na, nb in ((4, 8), (8, 16), (16, 32)):
            for ma, mb in ((1, 2), (2, 4)):
                name = f'N{na}-{nb}-M{ma}-{mb}'
                effects, raw = {}, {}
                for key in SELECTED:
                    d = values(table[nb, mb], key)-values(table[na, mb], key)-values(table[nb, ma], key)+values(table[na, ma], key)
                    raw[key] = d.tolist()
                    g = values(table[nb, ma], key)-values(table[na, ma], key)
                    effects[key] = dict(metric=metric(d, g), interaction_dot_global=dot(d, g))
                entry['interactions'][name] = effects
                vectors['cross/'+engine+'/interaction/'+name] = raw
        result['cross'][engine] = entry
    result['work_differences'] = {}
    for spec in plan:
        if len(spec['engines']) == 2:
            left, right = (result['raw'][spec['name']+'/'+e] for e in ('native', 'independent'))
            result['work_differences'][spec['name']] = [dict(row=a['row'], **{k: a[k]-b[k] for k in (
                'endpoint_work_kn_m', 'numerical_dissipation_kn_m', 'elastic_energy_kn_m', 'load_work_cumulative_kn_m')}) for a, b in zip(left, right)]
    result['s1_s5_changed'] = False
    result['new_limited_acceptance'] = False
    result['source_binding_unchanged'] = driver['sources_unchanged_after_solves']
    write(vector_output or output.with_name(output.stem+'-vectors.json'), vectors)
    write(output, result)
    print('formal paths', result['formal_paths'], 'states', result['formal_states'], flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--vectors', type=Path)
    a = p.parse_args()
    assess(a.run, a.output, a.vectors)
