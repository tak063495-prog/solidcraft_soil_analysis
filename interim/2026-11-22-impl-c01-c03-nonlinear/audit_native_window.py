"""Audit the prefix emitted by joint_path_window without requiring a completed run."""
import argparse
from collections import Counter
import json
from pathlib import Path


def summarize(run):
    case = run['case']
    settings = case['settings']
    steps = round(1 / settings['initial_increment'])
    assert settings['initial_increment'] == 1 / steps
    phases = [('geostatic', run['geostatic']['attempts'])]
    for stage in case['stages']:
        assert stage['initial_increment'] == 1 / steps
        phases.append((stage['name'], [a for a in run['attempts'] if a['stage'] == stage['name']]))
    maximum_force = maximum_linear = 0.
    stages = {}
    for name, attempts in phases:
        if not attempts:
            continue
        coordinates = set()
        progress = 0.
        reasons = Counter()
        for trial in attempts:
            assert trial['progress_from'] == progress
            error = trial.get('joint_path_error')
            if trial['accepted']:
                assert error and not error['reject']
                assert error['maximum_slip_error_m'] <= error['slip_budget_m']
                assert error['maximum_dissipation_error_kpa_m'] <= error['dissipation_budget_kpa_m']
                progress = trial['progress_to']
                coordinates.add(progress)
            else:
                assert trial['committed_state_unchanged_on_rejection']
                reasons[trial['reason']] += 1
            reports = list(trial['linear_solves'])
            if error:
                forces = (trial['relative_force_residual'], error['midpoint_relative_force_residual'],
                          error['fine_relative_force_residual'])
                assert all(0 <= force <= settings['relative_force_tolerance'] for force in forces)
                maximum_force = max(maximum_force, *forces)
                reports.extend(error['midpoint_linear_solves'])
                reports.extend(error['fine_linear_solves'])
            for report in reports:
                residual = report['true_relative_residual']
                assert 0 <= residual <= settings['relative_linear_tolerance']
                maximum_linear = max(maximum_linear, residual)
        assert {i/steps for i in range(1, int(progress*steps)+1)} <= coordinates
        stages[name] = dict(attempts=len(attempts), accepted=len(coordinates), progress=progress,
                            rejected_by_reason=dict(reasons))
    return dict(steps=steps, control=case['joint_system']['global_step_control'],
                stop=run['failure_reason'], stages=stages, maximum_force=maximum_force,
                maximum_linear=maximum_linear, numerical_acceptance=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    args = parser.parse_args()
    print(json.dumps(summarize(json.loads(args.run.read_bytes())), indent=2, allow_nan=False))
