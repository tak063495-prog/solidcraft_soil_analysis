"""Read-only reassembly of new diagnostic paths; incomplete paths stay false."""
import argparse
import gzip
import json
from pathlib import Path
import sys
from collections import Counter

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from oblique import fractions, component_bounds
from linear_evidence_codec import read_rows as actual_linear_rows


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def rows(path):
    with gzip.open(path, 'rt', encoding='utf-8') as stream:
        return [json.loads(line) for line in stream]


def stress(states, kind):
    if kind == 'native':
        return np.asarray([s['effective_stress_compression_positive_kpa'] for s in states])
    return np.asarray(states).reshape(-1, 14)[:, :6]


def audit(folder):
    model = read(folder/'model.json'); controls = read(folder/'controls.json'); report = read(folder/'report.json')
    assert not report['passed'] and not report['full_release_acceptance'], 'diagnostic falsely certifies acceptance'
    B, R, w = (np.asarray(model[k]) for k in ('B', 'transform', 'weights'))
    H, final, first, second = (np.asarray(model[k]) for k in ('holding_kn', 'final_load_kn', 'first_group_kn', 'second_group_kn'))
    # Replay the original complement = H - clipped_mean operation exactly.
    # Inverting that subtraction by addition need not recover H's last bit.
    np.testing.assert_array_equal(first, H-second)
    component_bounds(H, first)
    component_bounds(H, second)
    partition_roundoff = float(np.max(np.abs(first+second-H)))
    samples = rows(folder/'samples.jsonl.gz')
    initial = model['initial_states']; pi = np.asarray(model['active_gp']); dead = sorted(set(range(324))-set(pi))
    max_force = 0.; max_strain = 0.; minimum_h = float('inf')
    def frozen(states):
        if model['map'] == 'native':
            return all(states[i] == initial[i] for i in dead)
        return np.array_equal(np.asarray(states).reshape(-1, 14)[dead], np.asarray(initial).reshape(-1, 14)[dead])
    for i, sample in enumerate(samples):
        a, b = fractions(sample['coordinate'], model['waypoint'])
        assert sample['fractions'] == [a, b]
        holding = (1-a)*first+(1-b)*second
        component_bounds(H, holding)
        load = final+holding
        internal = -np.einsum('qai,qa,q->i', B, stress(sample['states'], model['map'])[pi], w)
        residual = float(np.linalg.norm(R.T@load-internal)/max(np.linalg.norm(load), 1.))
        assert residual < 1e-8, 'saved force equilibrium failed'
        assert frozen(sample['states']) and sample['contacts'] == model['contacts'], 'removed history changed'
        assert sample['removed_material_frozen'] and sample['contact_history_frozen']
        max_force = max(max_force, residual); minimum_h = min(minimum_h, float(np.linalg.norm(holding)))
        if i:
            du = np.asarray(sample['reduced_increment_m'])
            maximum = float(np.max(np.abs(np.einsum('qai,i->qa', B, du)).sum(axis=1)))
            assert maximum <= .001
            max_strain = max(max_strain, maximum)
            np.testing.assert_allclose(np.asarray(sample['displacement_m'])-samples[i-1]['displacement_m'],
                                       (R@du).reshape(-1, 3), atol=1e-12, rtol=0.)
    attempts = rows(folder/'attempts.jsonl.gz')
    zero = [a for a in attempts if a.get('label') == 'lambda-zero']
    if not controls['binding'].get('source'):
        raise ValueError('missing input/source binding')
    for attempt in attempts:
        if not attempt['accepted']:
            assert attempt['rejected_committed_state_unchanged']
    rejected = Counter(a['stop_reason'] for a in attempts if not a['accepted'])
    linear = rows(folder/'linear-systems.jsonl.gz') if (folder/'linear-systems.jsonl.gz').exists() else []
    # Every accepted Newton event retains the original linear gate; rejected
    # events may exceed it and are deliberately not counted as accepted solves.
    accepted_linear = [e['true_linear_residual'] for a in attempts if a['accepted'] for e in a['events']]
    if accepted_linear:
        assert max(accepted_linear) <= 1e-10
    else:
        assert report['accepted_increments'] == 0 and len(samples) == 1, 'missing accepted linear evidence'
    full = report['full_removal_reached']
    if full:
        assert samples[-1]['coordinate'] == 1. and np.array_equal(final+(1-1.)*H, final)
    tail_path = folder/'terminal-linear-systems.jsonl.f64.gz'
    tail = list(actual_linear_rows(tail_path)) if tail_path.exists() else []
    terminal = None
    if tail:
        actual = tail[-1]
        matrix = np.asarray(actual['matrix'])
        singular = np.linalg.svd(matrix, compute_uv=False)
        terminal = dict(condition_2=float(singular[0]/singular[-1]), singular_values=singular.tolist(),
                        label=actual['label'], iteration=actual['iteration'],
                        true_relative_residual=float(np.linalg.norm(matrix@actual['solution']-actual['rhs'])/max(np.linalg.norm(actual['rhs']), 1e-30)),
                        may_be_rejected_trial=True, physical_limit_confirmed=False,
                        scope='last recorded actual unsymmetrized matrix, not all trial matrices')
    return dict(passed=False, report=report, saved_fields=len(samples), original_component_bounds_verified=True,
                original_partition_subtraction_exact=True, inverse_addition_roundoff_kn=partition_roundoff,
                maximum_reassembled_force_residual=max_force, maximum_increment_strain_l1=max_strain,
                minimum_remaining_holding_norm_kn=minimum_h, removed_material_and_contacts_frozen=True,
                lambda_zero_actual_perturbed_newton=bool(zero and zero[0]['accepted']),
                lambda_zero_final_residual=zero[0]['events'][-1]['relative_force_residual'] if zero else None,
                rejected_attempt_reasons=dict(rejected), rejected_committed_states_unchanged=True,
                linear_summaries=len(linear), accepted_linear_maximum=max(accepted_linear) if accepted_linear else None,
                saved_actual_tail_matrices=len(tail), terminal_matrix=terminal,
                full_release=full, full_release_freeze_verified=False, completion_restart_verified=False,
                scope='new bounded oblique load-control paths; no production commit or physical failure verdict')


def path_samples(folder):
    observed = rows(folder/'samples.jsonl.gz')
    if folder.name.endswith('-extended'):
        origin = folder.with_name(folder.name.removesuffix('-extended'))
        prior = rows(origin/'samples.jsonl.gz')
        first = dict(prior[-1]); first.pop('reduced_increment_m', None)
        assert observed[0] == first, 'extension did not start at actual terminal field'
        return prior+observed[1:]
    return observed


def compare(a, b):
    x = {p['coordinate']: p for p in path_samples(a)}
    y = {p['coordinate']: p for p in path_samples(b)}
    common = sorted(x.keys() & y.keys())
    if not common:
        return dict(common_coordinates=0, comparison_available=False, complete_history_converged=False)
    result = {}
    for key in ('displacement_m', 'states', 'contacts'):
        if key == 'states' and isinstance(x[common[0]][key][0], dict):
            fields = ('effective_stress_compression_positive_kpa', 'plastic_strain_compression_positive', 'preconsolidation_pressure_kpa', 'void_ratio')
            for field in fields:
                result[field] = metrics([np.asarray([s[field] for s in x[t][key]]) for t in common],
                                        [np.asarray([s[field] for s in y[t][key]]) for t in common])
        elif key != 'contacts':
            result[key] = metrics([np.asarray(x[t][key]) for t in common], [np.asarray(y[t][key]) for t in common])
        else:
            result[key] = dict(exact=all(x[t][key] == y[t][key] for t in common))
    return dict(common_coordinates=len(common), fields=result, complete_history_converged=False)


def metrics(x, y):
    return dict(maximum_absolute=max(float(np.max(np.abs(a-b))) for a, b in zip(x, y)),
                maximum_relative_l2=max(float(np.linalg.norm(a-b)/max(np.linalg.norm(b), 1e-30)) for a, b in zip(x, y)))


def cross_maps(a, b):
    x = {p['coordinate']: p for p in path_samples(a)}
    y = {p['coordinate']: p for p in path_samples(b)}
    common = sorted(x.keys() & y.keys())
    if not common:
        return dict(common_coordinates=0, comparison_available=False, independent_completion_verified=False)
    def fields(sample, kind):
        states = sample['states']
        if kind == 'native':
            return dict(displacement_m=np.asarray(sample['displacement_m']),
                        stress=np.asarray([s['effective_stress_compression_positive_kpa'] for s in states]),
                        plastic_strain=np.asarray([s['plastic_strain_compression_positive'] for s in states]),
                        pc=np.asarray([s['preconsolidation_pressure_kpa'] for s in states]),
                        void=np.asarray([s['void_ratio'] for s in states]))
        state = np.asarray(states).reshape(-1, 14)
        return dict(displacement_m=np.asarray(sample['displacement_m']), stress=state[:, :6],
                    plastic_strain=state[:, 8:14], pc=np.exp(state[:, 6]), void=np.expm1(state[:, 7]))
    nx = [fields(x[t], 'native') for t in common]; iy = [fields(y[t], 'independent') for t in common]
    mx, my = read(a/'model.json'), read(b/'model.json')
    return dict(common_coordinates=len(common), fields={k: metrics([f[k] for f in nx], [f[k] for f in iy]) for k in nx[0]},
                holding_vectors=metrics([np.asarray(mx['holding_kn'])], [np.asarray(my['holding_kn'])]),
                same_coordinate_does_not_guarantee_identical_load=True, independent_completion_verified=False)


def measure(source):
    paths = {p.name: audit(p) for p in sorted(source.glob('*')) if (p/'report.json').is_file()}
    sensitivity = {}
    for kind in ('native', 'independent'):
        for name in ('a80-b20', 'a20-b80'):
            a, b = source/f'{kind}-{name}-128', source/f'{kind}-{name}-256'
            if a.is_dir() and b.is_dir():
                sensitivity[f'{kind}-{name}'] = compare(a, b)
            a, b = source/f'{kind}-{name}-128-extended', source/f'{kind}-{name}-256-extended'
            if a.is_dir() and b.is_dir():
                sensitivity[f'{kind}-{name}-continued'] = compare(a, b)
    cross = {}
    for name in ('a80-b20', 'a20-b80'):
        a, b = source/f'native-{name}-256', source/f'independent-{name}-256'
        if a.is_dir() and b.is_dir():
            cross[name] = cross_maps(a, b)
        a, b = source/f'native-{name}-256-extended', source/f'independent-{name}-256-extended'
        if a.is_dir() and b.is_dir():
            cross[name+'-continued'] = cross_maps(a, b)
    return dict(passed=False, full_removal_accepted=False, paths=paths, step_sensitivity=sensitivity, native_vs_input_only_independent=cross,
                independent_plastic_completion=False, full_release_and_freeze=False,
                p3b_p4_complete=False, structure_water_joint_complete=False, p6_complete=False)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.output.exists():
        p.error('new output required')
    args.output.write_text(json.dumps(measure(args.source), indent=2, allow_nan=False)+'\n', encoding='utf-8')
