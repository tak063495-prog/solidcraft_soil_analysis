"""Strict same-ID CSV/JSON diagnostic; no production schema/exporter changes."""
import argparse
import copy
import csv
import json
import math
from pathlib import Path
import sys
from reference import export

TABLES = ('hydraulic', 'samples', 'forces', 'budget')
STRINGS = {'id', 'body'}


def pairs(items):
    out = {}
    for key, value in items:
        if key in out:
            raise ValueError('duplicate JSON key: '+key)
        out[key] = value
    return out


def read_json(path):
    def reject(value):
        raise ValueError('non-finite JSON: '+value)
    return json.loads(Path(path).read_text(encoding='utf-8-sig'), object_pairs_hook=pairs,
                      parse_constant=reject)


def finite(value):
    if type(value) not in (float, int) or not math.isfinite(value):
        raise ValueError('finite numeric value required')
    return float(value)


def read(path):
    path = Path(path)
    if path.is_file():
        return read_json(path)
    data = read_json(path/'metadata.json')
    for table in TABLES:
        with (path/(table+'.csv')).open(encoding='utf-8-sig', newline='') as f:
            reader = csv.DictReader(f)
            if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
                raise ValueError('duplicate/missing CSV columns')
            rows = []
            for row in reader:
                if None in row or any(v is None or v == '' for v in row.values()):
                    raise ValueError('ragged/missing CSV value')
                rows.append({k: v if k in STRINGS else float(v) for k, v in row.items()})
            data[table] = rows
    return data


def write_csv(data, destination):
    path = Path(destination)
    path.mkdir(parents=True, exist_ok=False)
    metadata = {k: copy.deepcopy(v) for k, v in data.items() if k not in TABLES}
    (path/'metadata.json').write_text(json.dumps(metadata, indent=2, allow_nan=False)+'\n')
    for table in TABLES:
        with (path/(table+'.csv')).open('w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(data[table][0]))
            writer.writeheader()
            writer.writerows(data[table])


def indexed(rows):
    out = {}
    if type(rows) is not list or not rows:
        raise ValueError('nonempty table required')
    for row in rows:
        key = row['id']
        if type(key) is not str or not key or key in out:
            raise ValueError('invalid/duplicate ID')
        out[key] = row
    return out


def diagnose(data, reference):
    expected = reference['metadata']
    meta = data['metadata']
    if set(meta) != set(expected):
        raise ValueError('metadata keys must match explicit contract')
    for key, value in expected.items():
        if key == 'producer':
            if type(meta[key]) is not str or not meta[key]:
                raise ValueError('producer required')
        elif type(meta[key]) is not type(value) or meta[key] != value:
            raise ValueError('incompatible metadata: '+key)
    differences = dict(surface_m=dict(max_abs=abs(finite(data['surface_m'])-reference['surface_m']),
                                      tolerance='未定義', accepted=None))
    tables = {}
    for table in TABLES:
        actual, target = indexed(data[table]), indexed(reference[table])
        if set(actual) != set(target):
            raise ValueError(table+': missing/extra ID; no interpolation')
        tables[table] = actual
        errors = {}
        for key, row in target.items():
            got = actual[key]
            if set(got) != set(row):
                raise ValueError(table+': missing/extra columns')
            for field, value in row.items():
                if type(value) is str or field == 'z_m':
                    if field == 'z_m':
                        finite(got[field])
                    if got[field] != value:
                        raise ValueError('ID/body/coordinate mismatch; no projection')
                else:
                    errors.setdefault(field, []).append(finite(got[field])-value)
        differences[table] = {key: dict(max_abs=max(map(abs, values)),
                                        rms=math.sqrt(sum(v*v for v in values)/len(values)),
                                        tolerance='未定義', accepted=None)
                              for key, values in errors.items()}
    # Representation identity, allowing only binary64 arithmetic roundoff;
    # not a field-error acceptance tolerance. See README.
    for row in tables['samples'].values():
        scale = max(1., abs(row['total_kpa']), abs(row['effective_kpa']), abs(row['pressure_kpa']))
        if abs(row['total_kpa']-row['effective_kpa']-row['pressure_kpa']) > 64*math.ulp(scale):
            raise ValueError('total/effective/pore-pressure identity inconsistent')
    for row in tables['hydraulic'].values():
        pressure = 9.81*(row['head_m']-row['z_m'])
        if abs(row['pressure_kpa']-pressure) > 64*math.ulp(max(1.,abs(pressure),abs(row['pressure_kpa']))):
            raise ValueError('pressure head/total head identity inconsistent')
    closure = max(abs(finite(r['total_kpa'])-finite(r['effective_kpa'])-finite(r['pressure_kpa']))
                  for r in tables['samples'].values())
    residual = max(abs(r['effective_internal_kn']+r['pore_internal_kn']-r['weight_kn']
                       -r['pore_boundary_kn']-r['effective_boundary_kn'])
                   for r in tables['forces'].values())
    return dict(valid_contract=True, producer=meta['producer'], differences=differences,
                total_effective_identity_max_kpa=closure, nodal_equilibrium_max_kn=residual,
                outward_water_defect_m3_s=sum(r['outward_m3_s'] for r in tables['budget'].values()),
                accepted=None, a03_complete=False,
                pending=['A-03拡張の許容値は未定義', '本体実出力の真正性は別途検証が必要',
                         'MCC・旧履歴の施工時初期状態は未定義'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('output must be new')
    try:
        saved = read_json(args.reference)
        ref = export(saved['metadata']['n'], saved['metadata']['alpha'])
        # Never accept a silently edited reference file as a new truth.
        if saved != ref:
            raise ValueError('reference must match current independent calculation exactly')
        report = diagnose(read(args.input), ref)
        status = 2
    except (ValueError, KeyError, TypeError, OSError, OverflowError) as error:
        report = dict(valid_contract=False, reason=str(error), accepted=False, a03_complete=False)
        status = 1
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return status


if __name__ == '__main__':
    sys.exit(main())
