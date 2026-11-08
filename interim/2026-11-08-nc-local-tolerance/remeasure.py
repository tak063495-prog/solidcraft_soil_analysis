"""Recheck saved small-case evidence, without solving or altering checkpoints."""
import argparse
import collections
import copy
import json
from pathlib import Path
import sys

ROOT = next(p for p in Path(__file__).resolve().parents
            if (p/'Cargo.toml').exists() and (p/'tools').is_dir())
sys.path.insert(0, str(ROOT/'tools'))
import numpy as np
import fem_reference as fem
import probe_low_preload_nc as low
import verify_wet_structures as wet


def nc(folder):
    low.check_inputs(folder/'input', excavation_increment=1/128,
                     minimum_increment=1/4194304, dense_lu=True,
                     tangent_perturbation=1e-9)
    mesh, ground, case = [low.read(folder/'input'/(n+'.json')) for n in low.NAMES]
    metadata = low.read(folder/'experiment.json')
    assert metadata['local_newton_tolerance'] in ('1e-11', '1e-13')
    run = low.read(folder/'result/run.json')
    origin = run['frames'][2]['state']
    removed = np.flatnonzero(~np.array(run['last_committed']['active_elements']))
    offsets = np.r_[0, np.cumsum(np.bincount([p['element_index'] for p in run['integration_points']],
                                            minlength=len(mesh['tetrahedra'])))]
    removed_points = np.concatenate([np.arange(offsets[i], offsets[i+1]) for i in removed])
    counts = collections.Counter(); progress = []; final = None; observed_commits = set()
    maximum_live = 0
    for line in (folder/'result/trace.jsonl').open(encoding='utf-8'):
        row = json.loads(line); counts[row['kind']] += 1
        if row['kind'] == 'progress':
            progress.append(row['event'])
        if row['kind'] == 'iterate':
            old = row['committed']; observed_commits.add(old['stage_progress'])
            assert old['joints'] == origin['joints']
            live = sum(row['active'][p['owner_a']] and row['active'][p['owner_b']]
                       for p in run['joint_integration_points'])
            maximum_live = max(maximum_live, live)
            assert live == 0, 'this diagnosis requires inactive removal contact'
            final = {k: row[k] for k in ('stage', 'progress', 'iteration', 'relative_force_residual')}
            final['maximum_partition'] = max(row['partitions'])
        if row['kind'] in ('iterate', 'line_search') and 'states' in row:
            for i in removed_points:
                assert row['states'][i] == origin['material_states'][i]
    accepted = [a for a in run['geostatic']['attempts']+run['attempts'] if a['accepted']]
    linear = [s['true_relative_residual'] for a in accepted for s in a['linear_solves']]
    assert max(linear) <= case['settings']['relative_linear_tolerance']
    assert max(a['relative_force_residual'] for a in accepted) <= case['settings']['relative_force_tolerance']
    assert counts['iterate'] > 0, 'no detailed history evidence'
    assert not run['complete'], 'full release requires a different full-path acceptance audit'
    last_newton = next(e for e in reversed(progress) if e['kind'] == 'newton')
    return dict(passed=False, production_adopted=False, complete=run['complete'],
                local_newton_tolerance=metadata['local_newton_tolerance'],
                common_settings=case['settings'],
                failure_reason=run['failure_reason'], last_attempt=run['attempts'][-1],
                remaining_holding_fraction=1-run['last_committed']['stage_progress'],
                observation=low.inspect_run(mesh, ground, case, run),
                maximum_accepted_linear_residual=max(linear),
                maximum_accepted_force_residual=max(a['relative_force_residual'] for a in accepted),
                trace_event_counts=dict(counts), saved_iteration_maximum_active_contact_points=maximum_live,
                saved_committed_contact_history_frozen=True,
                removed_material_history_frozen_in_saved_iterates_and_candidates=True,
                detailed_trace_committed_progresses=sorted(observed_commits),
                last_saved_iterate=final, last_newton_event=last_newton,
                detailed_trace_is_full_path=False,
                input_sha256={p.name: low.sha(p) for p in sorted((folder/'input').iterdir())})


def joint(folder):
    mesh = wet.read(folder/'mesh.json')
    ground = fem.normalize_ground(wet.read(folder/'ground.json'))
    case = fem.normalize_stage_case(wet.read(folder/'stages.json'))
    native = wet.read(folder/'stage-results.json')
    expected = wet.reference(mesh, ground, case, water_modes=('soil',))
    rows = wet.compare(mesh, ground, case, native, expected)
    mutations = []
    for kind in ('open', 'slip', 'negative_slip', 'nan_slip', 'dissipation'):
        altered = copy.deepcopy(native)
        point = altered['frames'][4]['state']['joints']['points'][0]
        if kind in ('open', 'slip'): point['status'] = kind
        elif kind == 'negative_slip': point['accumulated_plastic_slip_m'] = -1e-12
        elif kind == 'nan_slip': point['plastic_slip_m'][0] = float('nan')
        else: point['plastic_dissipation_kpa_m'] = 1e-12
        try:
            wet.compare(mesh, ground, case, altered, expected)
        except AssertionError:
            mutations.append(dict(mutation=kind, accepted=False))
        else:
            raise AssertionError('false history accepted: '+kind)
    before = wet.read(folder/'before-guard.json')
    return dict(unchanged_closed_stick_fields_pass=True, frames=rows,
                before_guard_observations=before, after_guard=mutations,
                opening_sliding_or_general_plasticity_accepted=False,
                p3b_p4_rows_complete=False, p6_complete=False)


def measure(folder):
    return dict(base='c4a2b39d7f420036a8903163a653a6ed943e718c',
                nc={name: nc(folder/'nc'/name) for name in ('baseline', 'local13')},
                joint=joint(folder/'joint'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    low.write(args.output, measure(args.evidence))
