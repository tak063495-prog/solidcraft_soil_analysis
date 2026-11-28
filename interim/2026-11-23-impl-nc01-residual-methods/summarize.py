"""Bind native saved-system and fresh path results to their inputs and executable."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re

from prepare import METHODS, LEGACY

MODES = ('baseline', 'refined', *METHODS)
ROOT = Path(__file__).resolve().parents[2]
TOL = 1e-10


def sha(data):
    return hashlib.sha256(data).hexdigest()


def read(path):
    return json.loads(path.read_bytes().decode('utf-8-sig'))


def saved_observations(log):
    prior = read(LEGACY.parent / 'measurements.json')
    old = 'tests/fixtures/linear-refinement/nc96.jsonl'
    end = 'interim/2026-11-22-impl-linear-solver/refined-terminal.jsonl'
    systems = {}
    for name, metrics in [(old, prior['nc96_systems']), (end, [prior['refined_terminal']['system']])]:
        for index, (raw, diag) in enumerate(zip((ROOT / name).read_bytes().splitlines(keepends=True), metrics, strict=True)):
            fingerprint = sha(raw)
            if fingerprint != diag['sha256']:
                raise ValueError('saved system differs from pinned diagnosis')
            systems[name, index] = dict(file=name, index=index, sha256=fingerprint,
                                        kappa2=diag['kappa2'], sigma_min=diag['sigma_min'],
                                        sigma_next=diag['sigma_penultimate'], methods={})
    for line in log.read_text(encoding='utf-8-sig').splitlines():
        if 'NC01_REPLAY ' not in line:
            continue
        item = json.loads(line.partition('NC01_REPLAY ')[2])
        mode = item['method']
        table = systems[item['file'], item['index']]['methods']
        result = item['result']
        if mode in table or mode not in MODES or result['passed'] != (result['residual'] <= TOL):
            raise ValueError('invalid saved-system result')
        if not result['passed']:
            steps = re.search(r', steps (\d+),', result['reason'])
            result['iterations'] = 1 + (int(steps[1]) if steps else int(mode == 'refined'))
        table[mode] = dict(seconds=item['seconds'], **result)
    if len(systems) != 8 or any(set(s['methods']) != set(MODES) for s in systems.values()):
        raise ValueError('expected eight distinct systems with six observations each')
    return list(systems.values())


def path_observation(folder, mode):
    process = read(folder / f'{mode}-process.json')
    raw = (folder / mode / 'run.json').read_bytes()
    run = json.loads(raw)
    stage = 'stages.json' if mode == 'baseline' else f'stages-{mode}.json'
    case = read(folder / stage)
    baseline = read(folder / 'stages.json')
    expected = json.loads(json.dumps(baseline))
    if mode == 'refined':
        expected['settings']['dense_lu_refinement'] = True
    elif mode in METHODS:
        expected['settings']['dense_residual_method'] = mode
    if (case != expected or case['settings']['relative_linear_tolerance'] != TOL
            or run['provenance']['solver_binary_sha256'] != process['binary_sha256']
            or process['exit_code'] != 0 or process.get('budget_stop', False)):
        raise ValueError(f'incomparable path {mode}')
    attempts = run['attempts'] + (run.get('geostatic') or {}).get('attempts', [])
    reports = [r for a in attempts if a['accepted'] for r in a.get('linear_solves', [])]
    if not reports or any(not 0 <= r['true_relative_residual'] <= TOL for r in reports):
        raise ValueError('missing or invalid accepted linear reports')
    last = run['attempts'][-1]
    reason = run['failure_reason'] or ''
    residual = re.search(r'true linear residual ([.e+\d-]+)', reason)
    return dict(process=process, case_sha256=sha((folder / stage).read_bytes()),
                run_sha256=sha(raw), complete=run['complete'], accepted=sum(a['accepted'] for a in attempts),
                rejected=sum(not a['accepted'] for a in attempts),
                last_accepted=run['last_committed']['stage_progress'],
                last_stage=last['stage'], target=last['progress_to'], failure_reason=reason,
                terminal_residual=float(residual[1]) if residual else None,
                committed_unchanged=all(a['committed_state_unchanged_on_rejection'] for a in attempts if not a['accepted']),
                accepted_iterations=dict(sorted(Counter(r['iterations'] for r in reports).items())),
                max_accepted_residual=max(r['true_relative_residual'] for r in reports))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder', type=Path)
    parser.add_argument('replay_log', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    systems = saved_observations(args.replay_log)
    paths = {m: path_observation(args.folder, m) for m in MODES}
    if len({p['process']['binary_sha256'] for p in paths.values()}) != 1:
        raise ValueError('six paths must share one executable')
    eligible = [m for m in METHODS if all(s['methods'][m]['passed'] for s in systems)]
    data = dict(linear_tolerance=TOL, input_sha256={f: sha((args.folder / f).read_bytes())
                for f in ('mesh.json', 'ground.json', 'stages.json')},
                saved=systems, paths=paths, eligible_methods=eligible)
    with args.output.open('x', encoding='utf-8', newline='\n') as f:
        json.dump(data, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write('\n')


if __name__ == '__main__':
    main()
