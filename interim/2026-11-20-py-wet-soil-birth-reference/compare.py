"""Strict future CSV/JSON export diagnostic. Never grants A-03 acceptance.

Exit 1 = invalid/rejected, 2 = valid diagnostic, acceptance undefined.
The synthetic roundtrip is NOT a comparison against a production solver.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from reference import cases
from numerical import integrate_forces
from policies import mcc_invariants


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError('duplicate JSON key: '+key)
        result[key] = value
    return result


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'), object_pairs_hook=unique_object)


def read_export(path):
    path = Path(path)
    if path.is_file():
        return read_json(path)
    data = read_json(path/'metadata.json')
    data['materials'] = read_json(path/'materials.json')
    for table in ('nodes', 'forces', 'fluxes'):
        with (path/(table+'.csv')).open(encoding='utf-8', newline='') as stream:
            data[table] = list(csv.DictReader(stream))
    return data


def finite(value):
    if isinstance(value, bool) or value is None or value == '':
        raise ValueError('missing/non-numeric quantity')
    result = float(value)
    if not np.isfinite(result):
        raise ValueError('non-finite quantity')
    return result


def table(rows, key, expected_ids):
    result = {}
    for row in rows:
        identity = str(row[key])
        if identity in result:
            raise ValueError('duplicate '+key)
        result[identity] = row
    if set(result) != set(expected_ids):
        raise ValueError('missing/extra '+key)
    return result


def reject_fill(obj):
    if isinstance(obj, dict):
        for key, value in obj.items():
            if key in ('pressure_filled', 'pressure_interpolated', 'history_reused', 'pressure_clipped') and value is not False:
                raise ValueError('forbidden inferred/filled pressure or reused history: '+key)
            reject_fill(value)
    elif isinstance(obj, list):
        for value in obj:
            reject_fill(value)


def audit_materials(materials, before):
    """Trust boundary: before is a separately retained input snapshot, not output."""
    old = materials['old']
    if not before['old'] or old != before['old']:
        raise ValueError('old material history changed/missing')
    new = materials['new']
    if not new or len({x['id'] for x in new}) != len(new):
        raise ValueError('new material IDs missing/duplicate')
    old_ids = {x['id'] for x in before['old']}
    if len(old_ids) != len(before['old']) or set(x['id'] for x in new) != set(before['new_ids']):
        raise ValueError('material point inventory differs from independent input snapshot')
    for state in new:
        if state['id'] in old_ids or state['origin'] != 'explicit_fresh_input' or state['source_id'] is not None:
            raise ValueError('existing history transferred to new soil')
        if state['model'] not in ('MC', 'MCC'):
            raise ValueError('unknown model')
        for key in ('plastic_strain', 'equivalent_plastic_strain', 'plastic_work_kpa'):
            values = np.asarray(state[key], dtype=float)
            if not np.all(np.isfinite(values)) or np.any(values != 0):
                raise ValueError('aged plastic history in new soil')
        if np.asarray(state['plastic_strain']).shape != (6,) or state['last_return'] != 'elastic':
            raise ValueError('invalid fresh state')
        if state['model'] == 'MCC':
            mcc_invariants(state)  # fail closed on absent explicit pc/e/stress
    return ['new MCC BIRTH pc/e/initialization: 仕様未定義'] if any(x['model'] == 'MCC' for x in new) else []


def expected(case):
    column, flow = cases()[case]
    force = integrate_forces(column, flow)
    return column, flow, force


def diagnose(data, before, case):
    reject_fill(data)
    for flag in ('pressure_filled', 'pressure_interpolated', 'pressure_clipped', 'history_reused'):
        if data[flag] is not False:
            raise ValueError('explicit false provenance flag required: '+flag)
    if data['schema'] != 'wet-soil-birth-export/1' or data['case'] != case:
        raise ValueError('wrong schema/case')
    required = dict(length='m', pressure='kPa', force='kN', flux='m3/s')
    if data['units'] != required or data['head_kind'] != 'total_z_up' or data['stress_sign'] != 'compression_positive':
        raise ValueError('units/sign/head convention mismatch')
    if data['complete'] is not True or data['status'] != 'completed':
        raise ValueError('incomplete/timeout/cancelled output')
    if data['generation'] != 1 or data['stage'] != 'birth_initial_state':
        raise ValueError('wrong birth generation/stage')
    pending = audit_materials(data['materials'], before)
    c, flow, ref = expected(case)
    ids = [str(i) for i in range(len(ref['nodes']))]
    nodes = table(data['nodes'], 'id', ids)
    forces = table(data['forces'], 'id', ids)
    errors = {key: [] for key in ['head_m', 'pressure_kpa', 'effective_vertical_kpa', 'effective_horizontal_kpa']}
    origin = 'explicit_hydrostatic' if flow is None else 'steady_darcy_solution'
    for identity, (x, z) in zip(ids, ref['nodes']):
        row = nodes[identity]
        if finite(row['x']) != x or finite(row['z']) != z or row['body'] != ref['bodies'][int(identity)]:
            raise ValueError('coordinates/body mapping mismatch')
        if row['pressure_origin'] != origin:
            raise ValueError('undefined/filled pressure provenance')
        p = float(c.pressure(z) if flow is None else flow.pressure(z))
        # Above the hydrostatic water table, head is explicitly inapplicable,
        # not filled from water level. Pressure=0 is a specified dry-zone model.
        href = None if flow is None and z > c.water else float(c.water if flow is None else flow.head(z))
        if href is None:
            if row['head_m'] not in (None, '') or row['head_applicable'] not in (False, 'false'):
                raise ValueError('unsaturated head was invented')
        else:
            if row['head_applicable'] not in (True, 'true'):
                raise ValueError('missing applicable head')
            errors['head_m'].append(finite(row['head_m'])-href)
        f = c.fields(z, p)
        for key, wanted in [('pressure_kpa', p), ('effective_vertical_kpa', f['effective_vertical']),
                            ('effective_horizontal_kpa', f['effective_horizontal'])]:
            errors[key].append(finite(row[key])-float(wanted))
    observed_forces = {}
    for key, wanted in ref['forces'].items():
        observed = np.array([[finite(forces[i][key+'_'+axis]) for axis in ('x', 'z')] for i in ids])
        errors[key+'_kn'] = (observed-wanted).ravel().tolist()
        observed_forces[key] = observed
    boundary = observed_forces['pressure_exterior']+observed_forces['pressure_interface']
    residual = (observed_forces['effective_internal']+observed_forces['pore_internal']
                -observed_forces['body']-boundary-observed_forces['effective_boundary'])
    external = observed_forces['body']+boundary+observed_forces['effective_boundary']
    fluxes = table(data['fluxes'], 'boundary', ['bottom', 'top'])
    outward = [finite(fluxes[key]['outward_m3_s']) for key in ('bottom', 'top')]
    q = 0. if flow is None else flow.flux
    errors['boundary_flux_m3_s'] = [outward[0]+q, outward[1]-q]
    metrics = {key: dict(max_abs=float(max(np.abs(value), default=0.)),
                         rms=float(np.sqrt(np.mean(np.square(value)))) if value else None,
                         tolerance='未定義', accepted=None) for key, value in errors.items()}
    return dict(schema='wet-soil-birth-comparison/1', case=case, valid_export=True,
                accepted=None, a03_complete=False, tolerance='未定義', metrics=metrics,
                water_defect_m3_s=float(sum(outward)),
                nodal_equilibrium_max_kn=float(np.max(abs(residual))),
                global_force_kn=external.sum(axis=0).tolist(),
                global_moment_kn_m=float(np.sum(ref['nodes'][:, 0]*external[:, 1]-ref['nodes'][:, 1]*external[:, 0])),
                mcc_material_diagnostics=[mcc_invariants(s) for s in data['materials']['new'] if s['model'] == 'MCC'],
                pending=['A-03 production BIRTH acceptance', 'old material constitutive evolution: 仕様未定義', *pending],
                note='valid_export is schema/provenance validation, not numerical acceptance')


def synthetic(case):
    c, flow, ref = expected(case)
    old = [dict(id='old-gp-0', model='MC', plastic_strain=[.001, 0., 0., 0., 0., 0.],
                equivalent_plastic_strain=.001, plastic_work_kpa=.02)]
    before = dict(old=old, new_ids=['new-gp-0'])
    data = dict(schema='wet-soil-birth-export/1', case=case, producer='SYNTHETIC_REFERENCE_ONLY',
                units=dict(length='m', pressure='kPa', force='kN', flux='m3/s'),
                head_kind='total_z_up', stress_sign='compression_positive',
                complete=True, status='completed', generation=1, stage='birth_initial_state',
                pressure_filled=False, pressure_interpolated=False, pressure_clipped=False, history_reused=False,
                materials=dict(old=old, new=[dict(id='new-gp-0', model='MC', origin='explicit_fresh_input',
                    source_id=None, plastic_strain=[0.]*6, equivalent_plastic_strain=0.,
                    plastic_work_kpa=0., last_return='elastic')]), nodes=[], forces=[], fluxes=[])
    for i, ((x, z), body) in enumerate(zip(ref['nodes'], ref['bodies'])):
        p = float(c.pressure(z) if flow is None else flow.pressure(z))
        f = c.fields(z, p)
        applicable = bool(flow is not None or z <= c.water)
        data['nodes'].append(dict(id=str(i), x=float(x), z=float(z), body=body,
            head_applicable=applicable, head_m=float(c.water if flow is None else flow.head(z)) if applicable else None,
            pressure_kpa=p, effective_vertical_kpa=float(f['effective_vertical']),
            effective_horizontal_kpa=float(f['effective_horizontal']),
            pressure_origin='explicit_hydrostatic' if flow is None else 'steady_darcy_solution'))
        row = dict(id=str(i))
        for key, value in ref['forces'].items():
            for j, axis in enumerate(('x', 'z')):
                row[key+'_'+axis] = float(value[i, j])
        data['forces'].append(row)
    q = 0. if flow is None else flow.flux
    data['fluxes'] = [dict(boundary='bottom', outward_m3_s=-q), dict(boundary='top', outward_m3_s=q)]
    return data, before


def write_csv_bundle(data, path):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    meta = {k: v for k, v in data.items() if k not in ('nodes', 'forces', 'fluxes', 'materials')}
    (path/'metadata.json').write_text(json.dumps(meta, indent=2)+'\n', encoding='utf-8')
    (path/'materials.json').write_text(json.dumps(data['materials'], indent=2)+'\n', encoding='utf-8')
    for key in ('nodes', 'forces', 'fluxes'):
        with (path/(key+'.csv')).open('w', newline='', encoding='utf-8') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(data[key][0]))
            writer.writeheader()
            for row in data[key]:
                writer.writerow({k: str(v).lower() if isinstance(v, bool) else v for k, v in row.items()})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=cases(), required=True)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--history-before', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        result = diagnose(read_export(args.input), read_json(args.history_before), args.case)
        code = 2
    except (ValueError, KeyError, TypeError, OSError, IndexError) as error:
        result = dict(valid_export=False, accepted=False, a03_complete=False, rejection=str(error))
        code = 1
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    return code


if __name__ == '__main__':
    raise SystemExit(main())
