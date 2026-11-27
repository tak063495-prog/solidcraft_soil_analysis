"""Strict future-output comparison; no undefined tolerance is turned into pass.

JSON packet or CSV directory (metadata.json, nodes.csv, sections.csv).
No automatic units, frame conversion, interpolation, clipping or fallback.
"""
import argparse
import copy
import csv
import json
from pathlib import Path
import sys
import numpy as np

NODE_FIELDS = ['x_m', 'y_m', 'z_m', 'ux_m', 'uy_m', 'uz_m',
               'rx_rad', 'ry_rad', 'rz_rad', 'fx_kn', 'fy_kn', 'fz_kn',
               'mx_kn_m', 'my_kn_m', 'mz_kn_m']
SECTION_FIELDS = ['r', 's', 'Nxx_kn_m', 'Nyy_kn_m', 'Nxy_kn_m',
                  'Mpxx_kn', 'Mpyy_kn', 'Mpxy_kn', 'Msxx_kn', 'Msyy_kn',
                  'Msxy_kn', 'Qx_kn_m', 'Qy_kn_m']
UNITS = dict(length='m', rotation='rad', force='kN', moment='kN*m',
             membrane='kN/m', bending='kN', pressure='kPa')


def metadata(e=.06, case='eccentric-cantilever', n=2):
    return dict(schema='tri6-offset-reference-output/1',
                requested_formulation='mindlin_tri6_offset', actual_formulation='mindlin_tri6_offset',
                completed=True, cancelled=False, source_kind='independent-reference',
                units=UNITS.copy(), coordinate_system='local-right-handed-XYZ',
                basis_global=[[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], normal=[0., 0., 1.],
                displacement_surface='attachment', section_moment_surface='neutral-and-attachment',
                offset_definition='x_p=x_s-e*n;u_p=u_s+e*[n]cross*theta;Ms=Mp-e*N',
                stress_sign='tension-positive', reaction_sign='support-on-plate',
                install_reference='zero-displacement-zero-rotation', e_m=e, thickness_m=.2,
                E_kpa=2.e5, nu=0., case=case, subdivisions=n,
                load=dict(N_kn_m=3., Ms_kn=0.), boundary='x=0 all six DOF fixed',
                geometry='L=2m,H=1m,straight-planar-tri6')


def _finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and np.isfinite(value)


def validate(packet):
    if set(packet) != {'metadata', 'nodes', 'sections'}:
        raise ValueError('packet keys differ')
    m = packet['metadata']
    required = set(metadata())
    if set(m) != required:
        raise ValueError('explicit metadata missing/unknown: '+str(required ^ set(m)))
    if m['schema'] != metadata()['schema'] or m['units'] != UNITS:
        raise ValueError('schema or units mismatch')
    if m['requested_formulation'] != 'mindlin_tri6_offset' or m['actual_formulation'] != 'mindlin_tri6_offset':
        raise ValueError('legacy/default formulation fallback refused, including e=0')
    if m['completed'] is not True or m['cancelled'] is not False:
        raise ValueError('unfinished/cancelled result refused')
    if m['source_kind'] not in ('independent-reference', 'synthetic-test', 'production-export'):
        raise ValueError('explicit provenance required')
    if not all(_finite(m[k]) for k in ('e_m', 'thickness_m', 'E_kpa', 'nu')):
        raise ValueError('nonfinite material/offset')
    if m['thickness_m'] <= 0 or abs(m['e_m']) > m['thickness_m']/2:
        raise ValueError('|e|<=t/2 violated; no clipping')
    if m['E_kpa'] <= 0 or not -1 < m['nu'] < .5:
        raise ValueError('invalid material')
    for key in ('coordinate_system', 'basis_global', 'normal', 'displacement_surface',
                'section_moment_surface', 'offset_definition', 'stress_sign', 'reaction_sign', 'install_reference'):
        if m[key] != metadata()[key]:
            raise ValueError('offset sign/frame/reference mismatch: '+key)
    for table, fields in [('nodes', NODE_FIELDS), ('sections', SECTION_FIELDS)]:
        rows = packet[table]
        if not isinstance(rows, list) or not rows:
            raise ValueError('nonempty complete '+table+' required')
        seen = set()
        for row in rows:
            if set(row) != {'id', *fields} or not isinstance(row['id'], str) or not row['id'] or row['id'] in seen:
                raise ValueError('missing/duplicate/unknown point or field in '+table)
            seen.add(row['id'])
            if not all(_finite(row[k]) for k in fields):
                raise ValueError('nonfinite/non-numeric field in '+table)
    return packet


def read(path):
    path = Path(path)
    def unique(pairs):
        d = {}
        for k, v in pairs:
            if k in d: raise ValueError('duplicate JSON key: '+k)
            d[k] = v
        return d
    if path.is_dir():
        packet = dict(metadata=json.loads((path/'metadata.json').read_text(encoding='utf-8'), object_pairs_hook=unique))
        for table, fields in [('nodes', NODE_FIELDS), ('sections', SECTION_FIELDS)]:
            with (path/(table+'.csv')).open(encoding='utf-8-sig', newline='') as f:
                reader = csv.DictReader(f)
                if reader.fieldnames != ['id', *fields]:
                    raise ValueError('CSV headers differ: '+table)
                rows = list(reader)
            if any(None in row or any(v is None for v in row.values()) for row in rows):
                raise ValueError('CSV missing/extra cells')
            packet[table] = [dict(id=row['id'], **{k: float(row[k]) for k in fields}) for row in rows]
    else:
        packet = json.loads(path.read_text(encoding='utf-8'), object_pairs_hook=unique)
    return validate(packet)


def write_csv(packet, path):
    validate(packet); path = Path(path); path.mkdir(parents=True, exist_ok=False)
    (path/'metadata.json').write_text(json.dumps(packet['metadata'], indent=2, allow_nan=False)+'\n', encoding='utf-8')
    for table, fields in [('nodes', NODE_FIELDS), ('sections', SECTION_FIELDS)]:
        with (path/(table+'.csv')).open('w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['id', *fields]); writer.writeheader(); writer.writerows(packet[table])


def compare(reference, candidate):
    validate(reference); validate(candidate)
    a, b = copy.deepcopy(reference['metadata']), copy.deepcopy(candidate['metadata'])
    a.pop('source_kind'); b.pop('source_kind')
    if a != b:
        raise ValueError('input/offset/material/mesh/boundary contract mismatch')
    result = dict(accepted=None, tolerance='未定義',
                  scope='B-02 future diagnostic; not production acceptance',
                  production_export=candidate['metadata']['source_kind'] == 'production-export', tables={})
    for table, fields in [('nodes', NODE_FIELDS), ('sections', SECTION_FIELDS)]:
        ra, rb = ({row['id']: row for row in packet[table]} for packet in (reference, candidate))
        if ra.keys() != rb.keys():
            raise ValueError('point inventory mismatch; interpolation forbidden: '+table)
        ids = sorted(ra); report = {}
        for k in fields:
            x, y = (np.array([rows[i][k] for i in ids]) for rows in (ra, rb))
            if k in ('x_m', 'y_m', 'z_m', 'r', 's') and not np.array_equal(x, y):
                raise ValueError('coordinate/sample mismatch: '+k)
            error = y-x
            report[k] = dict(max_absolute=float(np.max(abs(error))), rms=float(np.sqrt(np.mean(error**2))),
                             relative_L2=None if np.linalg.norm(x) == 0 else float(np.linalg.norm(error)/np.linalg.norm(x)),
                             signed_differences=error.tolist())
        # A packet may not use the opposite section reference and merely relabel it.
        if table == 'sections':
            e = b['e_m']
            defect = np.array([[row[m]-row[p]+e*row[n] for n, p, m in
                               [('Nxx_kn_m', 'Mpxx_kn', 'Msxx_kn'), ('Nyy_kn_m', 'Mpyy_kn', 'Msyy_kn'),
                                ('Nxy_kn_m', 'Mpxy_kn', 'Msxy_kn')]] for row in rb.values()])
            # Roundoff-only identity check, not an engineering acceptance tolerance.
            scale = max(1., max(abs(row[k]) for row in rb.values() for k in fields))
            if np.max(abs(defect)) > 128*np.finfo(float).eps*scale:
                raise ValueError('section offset sign/reference identity violated: Ms != Mp-e*N')
        result['tables'][table] = dict(point_ids=ids, fields=report)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--candidate', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    try:
        report = compare(read(args.reference), read(args.candidate))
        status = 2  # A valid diagnostic cannot assert an undefined acceptance.
    except (ValueError, OSError, TypeError, KeyError) as exc:
        report, status = dict(accepted=False, rejected=True, reason=str(exc)), 1
    args.output.parent.mkdir(exist_ok=True, parents=True)
    with args.output.open('x', encoding='utf-8') as f:
        f.write(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False)+'\n')
    return status


if __name__ == '__main__':
    sys.exit(main())
