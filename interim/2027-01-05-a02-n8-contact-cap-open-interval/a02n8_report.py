"""Fixed-order observations from new equilibrium solutions and identified old controls."""
import argparse
from collections import Counter
import json
from pathlib import Path
import numpy as np
from a02n8_binding import study, fields, metric, dot, norm, common, switches as legacy_switches, switchdiff
from a02n8_common import write, sha

QUANTITIES = ('endpoint_work', 'traction', 'plastic_slip', 'accumulated', 'plastic_dissipation', 'numerical_dissipation', 'u')


def switches(rows):
    events = legacy_switches(rows)
    for i, ((ra, a), (rb, b)) in enumerate(zip(rows, rows[1:])):
        if not np.all(a['joint_active']) or not np.all(b['joint_active']) or not np.array_equal(a['generation'], b['generation']):
            continue
        aa, bb = np.sign(a['gap']), np.sign(b['gap'])
        for gp in np.flatnonzero(aa != bb):
            events.append(dict(kind='gap_sign', stage=rb['stage'], gp=int(gp), before=int(aa[gp]), after=int(bb[gp]),
                               lower=ra['progress'] if ra['stage'] == rb['stage'] else 0., upper=rb['progress'],
                               accepted_index=i+1, generation=b['generation'].tolist(), exact_internal_switch_root_known=False))
    return events


def load(folder, engine):
    return list(fields.load_fields(folder/(engine+'-fields.jsonl')))


def select(rows, grid):
    if grid == 'common_i_over_2':
        return common(rows, 2)
    if grid == 'stage_endpoints':
        return [v for v in rows if v[0]['progress'] == 1.]
    raise ValueError(grid)


def effect(left, right, label, vectors, grid='common_i_over_2', denominator='right'):
    a, b = select(left, grid), select(right, grid)
    coords = [(r['stage'], r['progress']) for r, _ in a]
    assert coords == [(r['stage'], r['progress']) for r, _ in b]
    result = {}
    for key in QUANTITIES:
        av, bv = np.stack([f[key] for _, f in a]), np.stack([f[key] for _, f in b])
        d = (av-bv).ravel(order='C')
        if denominator not in ('left', 'right'):
            raise ValueError('explicit denominator side required')
        result[key] = metric(av, bv) if denominator == 'right' else metric(bv, av)
        vectors[label+'/'+key] = dict(shape=list(av.shape), coordinates=coords, difference_left_minus_right=d.tolist())
    result['coordinates'] = coords
    result['grid'] = grid
    result['denominator_side'] = denominator
    result['signed_vector_convention'] = 'left_minus_right'
    return result


def branch_summary(events):
    kinds = Counter()
    for e in events:
        if e['kind'] == 'branch':
            kind = 'opening' if e['after'] == 0 else ('closure' if e['before'] == 0 else 'stick_slip_endpoint')
        else:
            kind = e['kind']
        kinds[kind] += 1
    return dict(counts=dict(sorted(kinds.items())),
                unique_slip_activity_gp=len({e['gp'] for e in events if e['kind'] == 'slip_activity'}),
                slip_activity_interval_gp_occurrences=kinds['slip_activity'],
                open_upper_endpoints=sorted({e['upper'] for e in events if e['stage'] == 'open' and e['kind'] == 'branch' and e['after'] == 0}),
                slip_activity_upper_endpoints=sorted({e['upper'] for e in events if e['stage'] == 'open' and e['kind'] == 'slip_activity'}),
                gap_bool_duplicates_opening_closure=True,
                slip_activity_can_overlap_opening_and_is_not_an_extra_physical_switch=True,
                exact_continuous_root=False)


def traces(folder, engine, rows):
    path = folder/('native/events.jsonl' if engine == 'native' else 'events.jsonl')
    events = [json.loads(s) for s in path.read_text(encoding='utf-8').splitlines()]
    names = ['geostatic'] + [s['name'] for s in study.read(folder/'case.json')['stages']]
    groups = {}
    for e in events:
        if e['kind'] not in ('iterate', 'linear', 'accepted', 'repartition', 'a02n8_iterate', 'a02n8_repartition'):
            continue
        stage = names[e['stage']] if engine == 'native' else e['stage']
        groups.setdefault((stage, e['progress']), []).append(e)
    native_attempts = {}
    if engine == 'native':
        run = study.read(folder/'native/run.json')
        for a in run['geostatic']['attempts'] + run['attempts']:
            stage = a.get('stage', 'geostatic')
            if stage == 'geostatic_equilibration':
                stage = 'geostatic'
            native_attempts[stage, a['progress_to']] = a
    output = []
    for row, f in rows:
        key = (row['stage'], row['progress'])
        group = groups[key]
        iters = [v for v in group if v['kind'] in ('iterate', 'a02n8_iterate')]
        if engine == 'native':
            counts = iters[-1]['partitions'][-32:]
        else:
            counts = [v for v in group if v['kind'] == 'accepted'][-1]['counts']
        assert len(counts) == 32
        changes = []
        for before, after in zip(iters, iters[1:]):
            a = before.get('partitions', before.get('counts'))
            b = after.get('partitions', after.get('counts'))
            if a != b:
                changes.append(dict(before_iteration=before['iteration'], after_iteration=after['iteration'],
                                    before_parts=a, after_parts=b,
                                    contact_parts_changed=a[-32:] != b[-32:]))
        emitted = sum(v['kind'] in ('repartition', 'a02n8_repartition') for v in group)
        # The frozen native observer emits candidate-time restart events only.
        # Its deferred R=independent increases are visible in consecutive full
        # partition vectors, so zero emitted events is not zero repartitions.
        actual_repartitions = len(changes) if engine == 'native' else emitted
        output.append(dict(stage=key[0], progress=key[1], contact_parts=counts,
                           joint_active=f['joint_active'].tolist(),
                           actual_contact_parts_if_active=counts if np.all(f['joint_active']) else None,
                           newton_iterate_evaluations_including_converged=len(iters),
                           newton_corrections=(len(native_attempts[key]['linear_solves']) if engine == 'native'
                                               else sum(v['kind'] == 'linear' for v in group)),
                           native_reported_newton_iterations=native_attempts[key]['newton_iterations'] if engine == 'native' else None,
                           repartitions=actual_repartitions,
                           emitted_repartition_events=emitted,
                           partition_changes_between_iterates=changes,
                           repartition_count_source='full consecutive iterate partitions' if engine == 'native' else 'explicit observer events including bulk restarts',
                           last_force_residual=iters[-1].get('relative', iters[-1].get('relative_force_residual')),
                           native_attempt=native_attempts.get(key)))
    return dict(accepted_points=output,
                maximum_contact_parts=max(max(r['actual_contact_parts_if_active']) for r in output if r['actual_contact_parts_if_active'] is not None),
                newton_corrections=sum(r['newton_corrections'] for r in output),
                native_reported_newton_iterations=sum(r['native_reported_newton_iterations'] for r in output) if engine == 'native' else None,
                newton_iterate_evaluations_including_converged=sum(r['newton_iterate_evaluations_including_converged'] for r in output),
                repartitions=sum(r['repartitions'] for r in output),
                maximum_accepted_force_residual=max(r['last_force_residual'] for r in output))


def open_observations(rows):
    events = switches(rows)
    intervals = []
    for (ra, a), (rb, b) in zip(rows, rows[1:]):
        if rb['stage'] != 'open':
            continue
        increment = b['accumulated']-a['accumulated']
        intervals.append(dict(lower=ra['progress'] if ra['stage'] == 'open' else 0., upper=rb['progress'],
                              slip_activity_gp=int(np.count_nonzero(increment > 0.)),
                              accumulated_slip_increment_m=increment.tolist(),
                              endpoint_status=b['contact_status'].tolist(), gap_m=b['gap'].tolist(),
                              opening_gp=int(np.count_nonzero((a['contact_status'] != 0) & (b['contact_status'] == 0))),
                              gap_bool_changes=int(np.count_nonzero(a['contact_open'] != b['contact_open'])),
                              signed_gap_changes=int(np.count_nonzero(np.sign(a['gap']) != np.sign(b['gap']))),
                              endpoint_work_kpa_m=b['endpoint_work'].tolist(), plastic_slip_m=b['plastic_slip'].tolist()))
    return dict(branch_summary=branch_summary(events), intervals=intervals,
                signed_gap_zero_negative_positive_are_distinct=True)


def report(run, old, out):
    out.mkdir(parents=True, exist_ok=False)
    paths = study.read(run/'paths.json')
    data, switches_all, trace, vectors = {}, {}, {}, {}
    complete = [p for p in paths if p.get('complete') and p['exit_code'] == 0]
    for p in complete:
        name = p['path']
        rows = load(run/name, p['engine'])
        data[name] = rows
        switches_all[name] = switches(rows)
        trace[name] = traces(run/name, p['engine'], rows)
    exchange = {}
    for e in ('native', 'independent'):
        for p in ('native', 'independent'):
            for r in ('native', 'independent'):
                suffix = f'{e}-p{p}-r{r}'
                cap, normal = 'exchange/cap/'+suffix, 'exchange/normal/'+suffix
                if cap not in data or normal not in data:
                    continue
                exchange[suffix] = dict(
                    cap_minus_normal=effect(data[cap], data[normal], 'cap_minus_normal/'+suffix, vectors),
                    cap_minus_normal_all_n8_points=effect_all(data[cap], data[normal], 'cap_minus_normal_all_n8_points/'+suffix, vectors),
                    switches_cap_vs_normal=switchdiff(switches_all[cap], switches_all[normal]),
                    cap_trace=trace[cap], normal_trace=trace[normal],
                    cap_branch_summary=branch_summary(switches_all[cap]))
    paired = {}
    for cap in ('cap', 'normal'):
        for p in ('native', 'independent'):
            for r in ('native', 'independent'):
                suffix = f'p{p}-r{r}'
                a, b = f'exchange/{cap}/native-{suffix}', f'exchange/{cap}/independent-{suffix}'
                if a in data and b in data:
                    paired[cap+'/'+suffix] = dict(
                        native_minus_independent=effect(data[a], data[b], 'paired/'+cap+'/'+suffix, vectors),
                        native_minus_independent_all_n8_points=effect_all(data[a], data[b], 'paired_all_n8_points/'+cap+'/'+suffix, vectors),
                        switches=switchdiff(switches_all[a], switches_all[b]))
    within = {}
    for e in ('native', 'independent'):
        base = f'exchange/cap/{e}-p{e}-r{e}'
        if base not in data:
            continue
        for p in ('native', 'independent'):
            for r in ('native', 'independent'):
                key = f'exchange/cap/{e}-p{p}-r{r}'
                if key in data:
                    within[key] = dict(fields=effect(data[key], data[base], 'within/'+key, vectors),
                                       switches=switchdiff(switches_all[key], switches_all[base]))
    baseline, open_results, grid_difference, old_n4_common = {}, {}, {}, {}
    for p in ('native', 'independent'):
        for r in ('native', 'independent'):
            suffix = f'p{p}-r{r}'
            cap_rows = {}
            for e in ('native', 'independent'):
                name = e+'-'+suffix
                cap_rows[e] = load(old/'retained/a02follow-contact-exchange'/name, e)
                normal = load(old/'formal/exchange/stress_free-cycle'/name, e)
                old_n4_common['cap_minus_normal/'+name] = effect(cap_rows[e], normal, 'old_n4/cap_minus_normal/'+name, vectors)
            old_n4_common['paired/'+suffix] = effect(cap_rows['native'], cap_rows['independent'], 'old_n4/paired/'+suffix, vectors)
    for e in ('native', 'independent'):
        for n in (4, 8):
            folder = old/f'formal/cross/n{n}-m4/stress_free-cycle/{e}'
            baseline[n, e] = load(folder, e)
            switches_all[f'old/n{n}/{e}'] = switches(baseline[n, e])
        grid_difference[e+'/old_n8_minus_n4'] = effect(baseline[8, e], baseline[4, e], e+'/old_n8_minus_n4', vectors, denominator='left')
        for b in (2, 4, 8):
            for n in (4, 8):
                key = f'open/n{n}-b{b}/{e}'
                if key not in data:
                    continue
                rows = data[key]
                previous = baseline[n, e] if b == 2 else data.get(f'open/n{n}-b{b//2}/{e}')
                basecase = study.read(old/f'formal/cross/n{n}-m4/stress_free-cycle/{e}/case.json')
                newcase = study.read(run/key/'case.json')
                basecase['settings']['minimum_increment'] = newcase['settings']['minimum_increment']
                assert basecase == newcase
                open_results[key] = dict(observations=open_observations(rows),
                    divided_minus_original=effect(rows, baseline[n, e], key+'/minus_original', vectors),
                    divided_minus_original_stage_endpoints=effect(rows, baseline[n, e], key+'/minus_original_stage_endpoints', vectors, 'stage_endpoints'),
                    divided_minus_previous=effect(rows, previous, key+'/minus_previous', vectors) if previous else None,
                    switches_vs_original=switchdiff(switches(rows), switches(baseline[n, e])))
            a, bb = f'open/n8-b{b}/{e}', f'open/n4-b{b}/{e}'
            if a in data and bb in data:
                grid_difference[f'{e}/b{b}_n8_minus_n4'] = effect(data[a], data[bb], f'{e}/b{b}_n8_minus_n4', vectors, denominator='left')
                original = np.stack([f['plastic_slip'] for _, f in common(baseline[8, e], 2)]) - np.stack([f['plastic_slip'] for _, f in common(baseline[4, e], 2)])
                newdiff = np.stack([f['plastic_slip'] for _, f in common(data[a], 2)]) - np.stack([f['plastic_slip'] for _, f in common(data[bb], 2)])
                grid_difference[f'{e}/b{b}_n8_minus_n4']['plastic_slip_difference_vs_original'] = dict(
                    dot=dot(newdiff, original), original_norm=norm(original), new_norm=norm(newdiff),
                    changed_difference_norm=norm(newdiff-original), causal_percentage=False)
    # Matched first-subinterval width: N4/b4 and N8/b2, N4/b8 and N8/b4.
    matched_width = {}
    for e in ('native', 'independent'):
        for b in (1, 2, 4):
            a, bb = f'open/n8-b{b}/{e}', f'open/n4-b{2*b}/{e}'
            left = baseline[8, e] if b == 1 else data.get(a)
            if left is not None and bb in data:
                matched_width[f'{e}/first_width_1_over_{8*b}'] = effect(left, data[bb], f'matched_width/{e}/{b}', vectors, denominator='left')
    for b in (2, 4, 8):
        for n in (4, 8):
            a, bb = f'open/n{n}-b{b}/native', f'open/n{n}-b{b}/independent'
            if a in data and bb in data:
                paired[f'open/n{n}-b{b}'] = dict(native_minus_independent=effect(data[a], data[bb], f'paired/open/n{n}-b{b}', vectors),
                                               switches=switchdiff(switches_all[a], switches_all[bb]))
    write(out/'difference-vectors.json', vectors)
    write(out/'branch-switches.json', switches_all)
    write(out/'accepted-traces.json', trace)
    write(out/'measurements.json', dict(passed=False, a02_complete=False, diagnostic_only=True,
        new_limited_acceptance=False, paths=paths, complete_principal_paths=len(complete),
        incomplete_paths=[p for p in paths if p not in complete],
        principal_accepted_states=sum(p['accepted_states'] for p in complete),
        observer_free_compatibility_paths=sum(p['engine'] == 'native' for p in complete),
        observer_free_accepted_states=sum(p['accepted_states'] for p in complete if p['engine'] == 'native'),
        exchange=exchange, paired=paired, within_engine=within, open=open_results,
        n4_n8=grid_difference, matched_first_subinterval_width=matched_width,
        old_n4_contact_exchange_summary=study.read(old/'source/contact-policy-exchange.json'),
        old_n4_contact_exchange_on_same_common_i_over_2=old_n4_common,
        input_binding_sha256=sha(old/'input-binding.json'),
        observations_are_bracket_upper_endpoints_not_exact_continuous_roots=True,
        norm_order='stage order, exact common i/2, GP then component, C ravel, math.fsum',
        denominator_floor=1e-12,
        n4_n8_denominator='finer N8 norm, as in the frozen 12-29 MxN effects; pairwise tables name their denominator side',
        tolerance_changed=False, windows_ci_after_push='未確定'))


def effect_all(left, right, label, vectors):
    # At N8 all accepted points are exact binary fractions; no interpolation.
    result = {}
    assert [(r['stage'], r['progress']) for r, _ in left] == [(r['stage'], r['progress']) for r, _ in right]
    coords = [(r['stage'], r['progress']) for r, _ in left]
    for key in QUANTITIES:
        av, bv = np.stack([f[key] for _, f in left]), np.stack([f[key] for _, f in right])
        result[key] = metric(av, bv)
        vectors[label+'/'+key] = dict(shape=list(av.shape), coordinates=coords, difference_left_minus_right=(av-bv).ravel(order='C').tolist())
    result.update(coordinates=coords, grid='all_exact_n8_endpoints')
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('run', 'old', 'output'):
        p.add_argument('--'+key, type=Path, required=True)
    a = p.parse_args()
    report(a.run.resolve(), a.old.resolve(), a.output.resolve())
