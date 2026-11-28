"""Report only accepted excavation progress from a solve-stages result or checkpoint."""
import argparse
import hashlib
import json
from pathlib import Path


def summarize(native, balance):
    run = native.get('run', native)
    frames = run['frames']
    settings = run['case']['settings']
    state = run['last_committed']
    attempts = run['attempts']
    cut = [a for a in attempts if a['stage'] == 'excavate_upper' and a['accepted']]
    audit = next((a for a in run['excavation_audits'] if a['stage'] == 'excavate_upper'), None)
    reached = cut[-1]['progress_to'] if cut else (0. if audit else None)
    linear = [s['true_relative_residual'] for a in attempts if a['accepted'] for s in a['linear_solves']]
    initial = frames[0] if frames else None
    geo = run.get('geostatic')
    physical_initial = initial is not None and (geo is None or geo['complete'])
    if physical_initial:
        tolerance = run['ground']['initialization']['relative_equilibrium_tolerance']
        scale = sum(abs(f) for f in balance['applied_force_kn'])
        if initial['relative_force_residual'] > tolerance:
            raise ValueError('native initial equilibrium failed')
        for key in ('applied_force_kn', 'support_reaction_kn'):
            if max(abs(a-b) for a, b in zip(initial[key], balance[key], strict=True)) > tolerance*scale:
                raise ValueError('independent initial resultant mismatch: '+key)
    if run['complete']:
        floor = frames[-1]
        active = state['active_elements']
        if (reached != 1. or floor['name'] != 'excavate_upper'
                or any(state['excavation_release_load_kn']) or sum(active)*2 != len(active)
                or audit['removed_element_count'] != len(active)-sum(active)
                or any(f['relative_force_residual'] > settings['relative_force_tolerance'] for f in frames)
                or max(linear, default=0.) > settings['relative_linear_tolerance']):
            raise ValueError('full removal is inconsistent with accepted equilibrium history')
    return {
        'complete': run['complete'], 'failure_reason': run['failure_reason'],
        'last_accepted_lambda': reached, 'initial_physical_equilibrium': physical_initial,
        'initial_free_relative_residual': initial['relative_force_residual'] if physical_initial else None,
        'independent_initial_balance': balance,
        'geostatic_strain_from_seed': geo['maximum_absolute_engineering_strain_from_seed'] if geo else None,
        'frames': [{k: f[k] for k in ('name', 'relative_force_residual', 'applied_force_kn',
                                       'support_reaction_kn', 'force_balance_error_kn')} for f in frames],
        'accepted_path': [[a['stage'], a['progress_to']] for a in attempts if a['accepted']],
        'rejected': [[a['stage'], a['progress_to'], a['reason']] for a in attempts if not a['accepted']],
        'removed_elements': audit['removed_element_count'] if audit else None,
        'zero_release_relative_residual': audit['zero_release_relative_force_residual'] if audit else None,
        'accepted_max_true_linear_residual': max(linear, default=None),
        'floor_mean_heave_during_excavation_m': (frames[-1]['monitors']['floor']['stage_increment_z_m']
                                                  if run['complete'] else None),
        'solver_binary_sha256': run['provenance']['solver_binary_sha256'],
    }


def read(path):
    return json.loads(path.read_bytes())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('result', type=Path)
    parser.add_argument('initial_balance', type=Path)
    args = parser.parse_args()
    data = summarize(read(args.result), read(args.initial_balance))
    data['native_file_sha256'] = hashlib.sha256(args.result.read_bytes()).hexdigest()
    print(json.dumps(data, indent=2, allow_nan=False))
