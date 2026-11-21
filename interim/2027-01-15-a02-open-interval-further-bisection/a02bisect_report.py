"""Observe accepted states only. A common coarse interval may span several fine steps."""
import argparse
import csv
import json
from pathlib import Path
import sys

from a02bisect_run import PREVIOUS, fields, grid, sha, write
sys.path.insert(0, str(PREVIOUS))
try:
    from a02n8_report import switches
    from a02n8_binding import metric, norm
finally:
    sys.path.remove(str(PREVIOUS))

import numpy as np

ENGINES = ('native', 'independent')
EVENTS = ('open', 'slip_endpoint', 'slip_activity')
KEYS = ('traction', 'endpoint_work', 'plastic_slip')


def read_path(folder, engine, level):
    data = list(fields.load_fields(folder / (engine+'-fields.jsonl')))
    endpoints = grid(level)['open']
    assert [s['progress'] for s, _ in data if s['stage'] == 'open'] == endpoints
    assert len(data) == 48 + len(endpoints)
    first = [{key: None for key in EVENTS} for _ in range(32)]
    for e in switches(data):
        if e['stage'] != 'open':
            continue
        if e['kind'] == 'branch' and e['after'] in (0, 2):
            key = 'open' if e['after'] == 0 else 'slip_endpoint'
        elif e['kind'] == 'slip_activity':
            key = 'slip_activity'
        else:
            continue
        gp = first[e['gp']]
        if gp[key] is None:
            gp[key] = (e['lower'], e['upper'])
    assert all(gp['open'] and gp['slip_activity'] for gp in first)
    states = {(s['stage'], s['progress']): f for s, f in data}
    assert len(states) == len(data)
    return states, first


def signed(left, right):
    diff = left - right
    return {'left_minus_right': diff.tolist(), 'l2': norm(diff)}


def interval(left, right, lo, hi):
    a, b = left['open', hi], right['open', hi]
    a0, b0 = left['open', lo], right['open', lo]
    return dict(lower=lo, upper=hi,
                traction_kpa=signed(a['traction'], b['traction']),
                increment_work_kpa_m=signed(a['endpoint_work'] - a0['endpoint_work'],
                                             b['endpoint_work'] - b0['endpoint_work']),
                cumulative_work_kpa_m=signed(a['endpoint_work'], b['endpoint_work']))


def path_comparison(left, right, lf, rf):
    points = sorted(p for s, p in left if s == 'open' and (s, p) in right and 2/64 <= p <= 4/64)
    rows = [interval(left, right, a, b) for a, b in zip(points, points[1:])]
    keys = [k for k in left if k[1] in (.5, 1.)]
    assert len(keys) == 14 and all(k in right for k in keys)
    common = {field: metric(np.stack([left[k][field] for k in keys]),
                            np.stack([right[k][field] for k in keys])) for field in KEYS}
    # Condition on equality of the observed upper endpoint, including paths with different step widths.
    groups = {}
    for label, same in [('same_first_open_upper', True), ('different_first_open_upper', False)]:
        ids = [i for i in range(32) if (lf[i]['open'][1] == rf[i]['open'][1]) == same]
        if not ids:
            groups[label] = dict(gp=[], common14_l2={})
            continue
        groups[label] = dict(gp=ids, common14_l2={
            field: norm(np.stack([left[k][field][ids] - right[k][field][ids] for k in keys]))
            for field in KEYS})
    return rows, dict(common14=common, by_open_upper=groups)


def report(root, out):
    out.mkdir(parents=True, exist_ok=False)
    paths, first, gp_rows, summaries, detail, scalar_rows = {}, {}, [], {}, {}, []
    for level in range(4):
        for engine in ENGINES:
            key = f'L{level}/{engine}'
            states, gps = read_path(root/key, engine, level)
            paths[key], first[key] = states, gps
            for i, gp in enumerate(gps):
                entry = dict(level=level, engine=engine, gp=i)
                for name in EVENTS:
                    lo, hi = gp[name] or (None, None)
                    entry[name+'_lower'], entry[name+'_upper'] = lo, hi
                gp_rows.append(entry)
    # Validate the copied control against fresh, input-only runs in both engines.
    for engine in ENGINES:
        old = root/'L0'/engine
        old_case = json.loads((old/'case.json').read_bytes())
        assert json.loads((old/'global-input-endpoints.json').read_bytes()) == grid(0)
        for level in (1, 2, 3):
            new = root/f'L{level}'/engine
            for name in ('mesh.json', 'ground.json'):
                assert json.loads((new/name).read_bytes()) == json.loads((old/name).read_bytes())
            case = json.loads((new/'case.json').read_bytes())
            assert case['settings']['minimum_increment'] == 1/(64*2**level)
            case['settings']['minimum_increment'] = old_case['settings']['minimum_increment']
            assert case == old_case
            assert json.loads((new/'global-input-endpoints.json').read_bytes()) == grid(level)
    with (out/'accepted-contact.jsonl').open('x', encoding='utf-8', newline='\n') as f:
        for name, states in paths.items():
            for (stage, progress), state in states.items():
                if stage == 'open' and 2/64 <= progress <= 4/64:
                    record = dict(path=name, upper=progress, **{field: state[field].tolist()
                        for field in ('contact_status', 'gap', 'accumulated', 'traction', 'endpoint_work')})
                    f.write(json.dumps(record, allow_nan=False, separators=(',', ':'))+'\n')
    def compare(label, left, right):
        rows, summary = path_comparison(paths[left], paths[right], first[left], first[right])
        detail[label] = rows
        summaries[label] = summary
        for row in rows:
            scalar_rows.append(dict(comparison=label, lower=row['lower'], upper=row['upper'],
                **{k+'_l2': v['l2'] for k, v in row.items() if isinstance(v, dict)}))
    for level in range(4):
        compare(f'L{level} native-independent', f'L{level}/native', f'L{level}/independent')
        if level:
            for engine in ENGINES:
                for base in sorted({0, level-1}):
                    compare(f'{engine} L{level}-L{base}', f'L{level}/{engine}', f'L{base}/{engine}')
    for path, rows in [('first-events.csv', gp_rows), ('interval-differences.csv', scalar_rows)]:
        with (out/path).open('x', encoding='utf-8', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
            w.writeheader()
            w.writerows(rows)
    write(out/'signed-differences.json', detail)
    counts = {}
    for path, gps in first.items():
        counts[path] = {}
        for event in EVENTS:
            counted = {}
            for row in gps:
                bracket = row[event]
                label = str(bracket[1]) if bracket else 'never'
                counted[label] = counted.get(label, 0) + 1
            counts[path][event] = counted
    write(out/'measurements.json', dict(passed=False, a02_complete=False, diagnostic_only=True,
          first_counts=counts, comparisons=summaries,
          observed_sha256={name: sha(out/name) for name in
                           ('first-events.csv', 'interval-differences.csv', 'signed-differences.json')}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--paths', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    report(args.paths, args.output)
