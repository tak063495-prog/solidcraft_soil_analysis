"""Run the independent input-built prefix through the requested cloud observation point."""
import argparse
import collections
import importlib.util
import json
from pathlib import Path
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--steps', type=int, choices=(512, 1024), required=True)
    parser.add_argument('--event-proposal', action='store_true')
    parser.add_argument('--local-substeps', type=int, choices=(1, 4), default=4)
    parser.add_argument('--native-input', action='store_true', help='write the matching native case without solving')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    spec = importlib.util.spec_from_file_location('path_study', here.parent / '2026-12-12-c01-switch-event-localization/study.py')
    study = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(study)
    inputs = here.parents[1] / 'examples'
    mesh = study.read(inputs / 'joint-water-column-tet10-mesh.json')
    ground = study.fem.normalize_ground(study.read(inputs / 'ground-joint-water-mc.json'))
    stages = study.fem.normalize_stage_case(study.read(inputs / 'stages-joint-water-cycle.json'))
    n = args.steps
    case = study.prior.configuration(stages, n, 5e-5)
    case['settings']['minimum_increment'] = min(case['settings']['minimum_increment'], 1/(n*65536))
    control = dict(study.PATH_ERROR)
    if args.event_proposal:
        control.update(event_proposal=True, minimum_contact_substeps=args.local_substeps)
    case['joint_system']['global_step_control'] = dict(path_error=control)
    if args.native_input:
        study.write(args.output, case)
        return
    name, threshold = ('reload', 472/512) if n == 512 else ('open', 44/1024)
    counts = collections.defaultdict(collections.Counter)
    maxima = collections.Counter()
    last = {}
    anchors = collections.defaultdict(set)
    started = time.monotonic()

    class Observed(Exception):
        pass

    def attempt(row):
        phase = row['stage']
        counts[phase]['attempts'] += 1
        counts[phase]['accepted' if row['accepted'] else 'rejected'] += 1
        if row['accepted']:
            audit = row['joint_path_error']
            assert not audit['reject']
            for key in ('maximum_slip_error_m', 'maximum_dissipation_error_kpa_m'):
                budget = audit['slip_budget_m' if key.endswith('error_m') else 'dissipation_budget_kpa_m']
                assert audit[key] <= budget
            for residual in (row['relative_force_residual'], audit['midpoint_relative_force_residual'], audit['fine_relative_force_residual']):
                assert residual < 1e-10
                maxima['force'] = max(maxima['force'], residual)

    def accepted(row, fields):
        phase = row['stage']
        progress = row['progress']
        last[phase] = progress
        if (progress * n).is_integer():
            anchors[phase].add(round(progress * n))
        if phase == name and progress >= threshold:
            raise Observed

    try:
        study.reference._reference(mesh, ground, case, n, 1, 7, False, True, accepted,
                                   maximum_m=5e-5, path_error=control, attempt_sink=attempt)
    except Observed:
        stop = 'observation point'
    except AssertionError as exc:
        if str(exc) != 'independent stage attempt budget':
            raise
        stop = str(exc)
    else:
        raise RuntimeError('observation point was missed')
    for phase, progress in last.items():
        assert anchors[phase] == set(range(1, int(progress*n)+1))
    study.write(args.output, dict(steps=n, event_proposal=args.event_proposal, stop=stop, phase=name, local_substeps=control.get('minimum_contact_substeps', 1),
                                  threshold=threshold, last=last, counts=dict(counts), maximum_force=maxima['force'],
                                  seconds=time.monotonic()-started, numerical_acceptance=False))


if __name__ == '__main__':
    main()
