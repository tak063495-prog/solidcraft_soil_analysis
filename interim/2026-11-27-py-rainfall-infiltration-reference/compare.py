"""Strict external CSV/JSON adapter. No interpolation, conversion or acceptance.

This is a proposed external interchange contract, not a production schema.
"""
import argparse
import csv
import json
from pathlib import Path
import math

QUANTITIES = {
    'head': {'unit': 'm', 'meaning': 'total_head', 'sign': 'elevation_up'},
    'saturation': {'unit': '1', 'meaning': 'theta/theta_s', 'sign': 'nonnegative'},
    'flux': {'unit': 'm/s', 'meaning': 'instantaneous_Darcy_flux', 'sign': 'down'},
    'infiltration': {'unit': 'm', 'meaning': 'cumulative_net_top_infiltration', 'sign': 'inward'},
    'effective_stress': {'unit': 'kPa', 'meaning': 'Terzaghi_vertical', 'sign': 'compression_positive'},
}
FIELDS = {'quantity', 'id', 'time_s', 'coordinate_m', 'value'}


def read(path):
    path = Path(path)
    if path.suffix == '.csv':
        meta = json.loads(path.with_suffix('.metadata.json').read_text(encoding='utf-8'))
        with path.open(encoding='utf-8', newline='') as f:
            reader = csv.DictReader(f)
            if set(reader.fieldnames or []) != FIELDS:
                raise ValueError('CSV columns must match exactly')
            rows = list(reader)
        for row in rows:
            for key in ('time_s', 'coordinate_m', 'value'):
                row[key] = float(row[key])
        return {'metadata': meta, 'rows': rows}
    if path.suffix != '.json':
        raise ValueError('only CSV or JSON is supported')
    return json.loads(path.read_text(encoding='utf-8'))


def validate(packet):
    if set(packet) != {'metadata', 'rows'}:
        raise ValueError('packet keys differ')
    m = packet['metadata']
    required = {'contract', 'case_id', 'parameters', 'source_kind', 'complete',
                'clock', 'event_trace', 'interpolation', 'quantities', 'times_s', 'locations',
                'hydraulic_model', 'mechanical_model'}
    if set(m) != required:
        raise ValueError('metadata keys differ')
    if m['contract'] != 'external-rainfall-reference/1' or m['complete'] is not True:
        raise ValueError('unsupported contract or incomplete computation')
    if m['source_kind'] not in ('independent_reference', 'production_output', 'synthetic'):
        raise ValueError('source_kind required')
    if m['clock'] != 'physical_seconds' or m['interpolation'] != 'none':
        raise ValueError('physical seconds and no interpolation required')
    if m['event_trace'] != 'right_limit_boundary_continuous_interior':
        raise ValueError('switch-time boundary trace must use the right limit')
    if m['quantities'] != QUANTITIES:
        raise ValueError('units, signs, head/saturation/stress/flow meaning differ')
    times = m['times_s']
    if not times or any(type(t) not in (float, int) or not math.isfinite(t) or t < 0 for t in times):
        raise ValueError('invalid times')
    if times != sorted(set(times)):
        raise ValueError('times must be unique and ordered')
    if set(m['locations']) != set(QUANTITIES):
        raise ValueError('missing quantity locations')
    expected = {}
    for quantity, locations in m['locations'].items():
        if not locations:
            raise ValueError('empty location set')
        ids = []
        for location in locations:
            if set(location) != {'id', 'coordinate_m'} or not isinstance(location['id'], str) or not location['id']:
                raise ValueError('explicit ID and coordinate required')
            if type(location['coordinate_m']) not in (float, int) or not math.isfinite(location['coordinate_m']):
                raise ValueError('invalid coordinate')
            ids.append(location['id'])
            for t in times:
                expected[quantity, location['id'], t] = location['coordinate_m']
        if len(ids) != len(set(ids)):
            raise ValueError('duplicate location ID')
    values = {}
    for row in packet['rows']:
        if set(row) != FIELDS:
            raise ValueError('row keys differ')
        for name in ('time_s', 'coordinate_m', 'value'):
            if type(row[name]) not in (float, int) or not math.isfinite(row[name]):
                raise ValueError('nonfinite or nonnumeric data')
        key = row['quantity'], row['id'], row['time_s']
        if key not in expected or key in values:
            raise ValueError('unknown/duplicate ID, quantity or exact time')
        if row['coordinate_m'] != expected[key]:
            raise ValueError('coordinate differs; no mapping or interpolation')
        if row['quantity'] == 'saturation' and not 0 <= row['value'] <= 1:
            raise ValueError('saturation outside [0,1]')
        values[key] = row['value']
    if values.keys() != expected.keys():
        raise ValueError('missing ID/time/quantity')
    return values


def compare(reference, candidate):
    a, b = validate(reference), validate(candidate)
    ma, mb = reference['metadata'].copy(), candidate['metadata'].copy()
    ma.pop('source_kind')
    mb.pop('source_kind')
    if ma != mb:
        raise ValueError('case, parameters, models, times or locations differ')
    if a.keys() != b.keys():
        raise ValueError('key sets differ')
    metrics = {}
    for q in QUANTITIES:
        differences = [b[k]-a[k] for k in a if k[0] == q]
        metrics[q] = {'count': len(differences),
                      'rms': math.sqrt(sum(d*d for d in differences)/len(differences)),
                      'maximum_absolute': max(abs(d) for d in differences),
                      'unit': QUANTITIES[q]['unit'], 'tolerance': None, 'accepted': None}
    return {'accepted': None, 'tolerance': '未定義', 'metrics': metrics,
            'production_comparison': candidate['metadata']['source_kind'] == 'production_output',
            'reason': '旧3Dゲートを本参照の5量へ適用する定義なし。差ゼロも受入ではない。'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference')
    parser.add_argument('candidate')
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    try:
        report = compare(read(args.reference), read(args.candidate))
        code = 2
    except (ValueError, KeyError, TypeError) as error:
        report = {'accepted': None, 'rejected': True, 'reason': str(error)}
        code = 1
    Path(args.out).write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    return code


if __name__ == '__main__':
    raise SystemExit(main())
