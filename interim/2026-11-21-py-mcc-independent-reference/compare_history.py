"""Strict CSV/JSON history comparison; no solver execution or history filling.
The monotone-pc exchange contract intentionally rejects dilative OC softening.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import numpy as np

GP = ('p', 'q', 'pc', 'e', 'eps_v')
NODE = ('ux', 'uy', 'uz')
META = dict(schema='independent-mcc-history-v1', stress_unit='kPa', length_unit='m',
            strain_unit='1', stress_sign='compression_positive',
            strain_sign='compression_positive', displacement_sign='global_xyz',
            stress_kind='effective', strain_measure='log_volume',
            history_origin='explicit_initial_state', history_completion='none',
            pc_policy='nondecreasing', run_status='completed')


def strict_object(pairs):
    result = {}
    for k, v in pairs:
        if k in result:
            raise ValueError('duplicate JSON key: '+k)
        result[k] = v
    return result


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'), object_pairs_hook=strict_object,
                      parse_constant=lambda x: (_ for _ in ()).throw(ValueError('nonfinite JSON '+x)))


def finite(value):
    if isinstance(value, bool):
        raise ValueError('boolean in numeric field')
    x = float(value)
    if not math.isfinite(x):
        raise ValueError('nonfinite numeric field')
    return x


def key(row):
    step = finite(row['step'])
    if step < 0 or int(step) != step or not isinstance(row['id'], str) or not row['id']:
        raise ValueError('invalid step/id')
    return int(step), row['id']


def validate(data):
    meta, rows = data['metadata'], data['rows']
    for k, expected in META.items():
        if meta.get(k) != expected:
            raise ValueError('metadata rejected: '+k)
    for k in ('case_id', 'physics_sha256', 'source_sha256'):
        if not isinstance(meta.get(k), str) or not meta[k]:
            raise ValueError('missing provenance: '+k)
    for k in ('physics_sha256', 'source_sha256'):
        if len(meta[k]) != 64 or any(c not in '0123456789abcdef' for c in meta[k]):
            raise ValueError('invalid SHA256: '+k)
    # A step is a declared load coordinate, not inferred physical time.
    steps = meta['load_coordinates']
    if not isinstance(steps, list) or len(steps) < 2:
        raise ValueError('explicit complete load coordinates required')
    steps = [finite(t) for t in steps]
    if steps[0] != 0 or any(b <= a for a, b in zip(steps, steps[1:])):
        raise ValueError('load coordinates must start at 0 and increase')
    kind = meta['kind']
    fields = GP if kind == 'gp' else NODE if kind == 'node' else None
    if fields is None:
        raise ValueError('unknown record kind')
    ids = meta['ids']
    if not isinstance(ids, list) or not ids or any(not isinstance(i, str) or not i for i in ids) or len(set(ids)) != len(ids):
        raise ValueError('explicit unique IDs required')
    expected_keys = {(s, i) for s in range(len(steps)) for i in ids}
    indexed = {}
    for row in rows:
        k = key(row)
        if k in indexed or k not in expected_keys:
            raise ValueError('duplicate/unexpected history key')
        if finite(row['load_coordinate']) != steps[k[0]]:
            raise ValueError('load coordinate mismatch; no interpolation')
        if set(row) != {'step', 'id', 'load_coordinate', *fields}:
            raise ValueError('missing/unknown history fields; no filling')
        indexed[k] = {f: finite(row[f]) for f in fields}
    if set(indexed) != expected_keys:
        raise ValueError('missing history; no interpolation/forward-fill/zero-fill')
    if kind == 'gp':
        for ident in ids:
            origin = indexed[(0, ident)]
            prev = None
            for step in range(len(steps)):
                row = indexed[(step, ident)]
                if row['p'] <= 0 or row['pc'] <= 0 or row['e'] <= 0 or row['q'] < 0:
                    raise ValueError('p<=0, pc<=0, e<=0 or q<0')
                if prev is not None and row['pc'] < prev['pc']:
                    raise ValueError('pc regression; exchange contract rejects softening')
                # Input consistency guard, NOT an acceptance tolerance (documented).
                ev = math.log((1+origin['e'])/(1+row['e']))
                if abs(ev-(row['eps_v']-origin['eps_v'])) > 1e-7:
                    raise ValueError('e / logarithmic volumetric strain inconsistent')
                prev = row
    return indexed, fields


def load(path, metadata_path=None):
    path = Path(path)
    if path.suffix.lower() == '.json':
        if metadata_path is not None:
            raise ValueError('JSON embeds metadata; sidecar forbidden')
        data = read_json(path)
    elif path.suffix.lower() == '.csv':
        if metadata_path is None:
            raise ValueError('CSV requires explicit --*-metadata JSON sidecar')
        with path.open(encoding='utf-8', newline='') as f:
            reader = csv.DictReader(f)
            if reader.fieldnames is None or len(reader.fieldnames) != len(set(reader.fieldnames)):
                raise ValueError('missing/duplicate CSV headers')
            rows = list(reader)
        data = dict(metadata=read_json(metadata_path), rows=rows)
    else:
        raise ValueError('only JSON/CSV supported')
    validate(data)
    return data


def compare(reference, candidate):
    a, fields = validate(reference)
    b, fields_b = validate(candidate)
    # Sources may differ; physical problem, IDs and load grid may not.
    am, bm = reference['metadata'], candidate['metadata']
    for k in set(am) | set(bm):
        if k != 'source_sha256' and am.get(k) != bm.get(k):
            raise ValueError('incomparable metadata: '+k)
    if fields != fields_b or a.keys() != b.keys():
        raise ValueError('incompatible histories')
    keys = sorted(a)
    metrics = {}
    for f in fields:
        ar = np.array([a[k][f] for k in keys])
        br = np.array([b[k][f] for k in keys])
        err, scale = np.linalg.norm(br-ar), np.linalg.norm(ar)
        metrics[f] = dict(max_abs=float(np.max(np.abs(br-ar))),
                          relative_l2=float(err/scale) if scale else None,
                          reference_l2=float(scale), accepted=None, tolerance='未定義')
    return dict(accepted=None, tolerance='未定義', status='diagnostic_only', metrics=metrics,
                reference_source_sha256=am['source_sha256'], candidate_source_sha256=bm['source_sha256'])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--candidate', type=Path, required=True)
    p.add_argument('--reference-metadata', type=Path)
    p.add_argument('--candidate-metadata', type=Path)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    try:
        a = load(args.reference, args.reference_metadata)
        b = load(args.candidate, args.candidate_metadata)
        result = compare(a, b)
        result['input_sha256'] = {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                                 for path in [args.reference, args.candidate,
                                              args.reference_metadata, args.candidate_metadata] if path}
        with args.output.open('x', encoding='utf-8') as f:
            json.dump(result, f, ensure_ascii=False, indent=2, allow_nan=False)
            f.write('\n')
    except (ValueError, KeyError, TypeError, OverflowError) as exc:
        p.exit(2, 'REJECTED: '+str(exc)+'\n')


if __name__ == '__main__':
    main()
