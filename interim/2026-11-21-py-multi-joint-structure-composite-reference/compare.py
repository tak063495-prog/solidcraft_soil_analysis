"""Fail-closed raw CSV/JSON diagnostic contract for a future production exporter.

No exporter or production solver is invoked. Undefined physical tolerances never
become acceptance. Data-only self comparison does not certify physical truth.
"""
import argparse
import copy
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from reference import SCHEMA, UNITS, SIGNS, dump

CITATION = '../2026-11-13-joint-ext-combos/acceptance.json:tolerances'
TOLERANCES = {'nodal_force_kn': 1e-10, 'work_kn_m': 1e-10,
              'independent_elastic_relative': 1e-7}
KINDS = ('nodes', 'contacts', 'members', 'elements')


class Refusal(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise Refusal(message)


def key_for(kind, row):
    return (row['joint_id'], row['gp_id'], row['generation']) if kind == 'contacts' else row['id']


def finite(value):
    if isinstance(value, dict):
        return all(finite(v) for v in value.values())
    if isinstance(value, list):
        return all(finite(v) for v in value)
    return not isinstance(value, float) or np.isfinite(value)


def validate(doc):
    require(doc.get('schema') == SCHEMA, 'SCHEMA')
    require(isinstance(doc.get('producer'), str) and bool(doc['producer']), 'PRODUCER_REQUIRED')
    require(doc.get('complete') is True and doc.get('frames'), 'INCOMPLETE')
    require(doc.get('units') == UNITS, 'UNITS; conversion must be explicit upstream')
    require(doc.get('signs') == SIGNS, 'SIGNS; no inferred reversal')
    require(doc.get('raw') is True and doc.get('transformations') == [],
            'CORRECTION_OR_UNATTESTED_RAW; smoothing/clipping/averaging forbidden')
    require(finite(doc), 'NONFINITE')
    require(doc['metadata']['water'] == 'absent' and doc['metadata']['generation'] == 0,
            'WET_OR_LIFECYCLE_OUT_OF_SCOPE')
    for kind in ('node_ids', 'joint_ids', 'member_ids'):
        require(len(doc[kind]) == len(set(doc[kind])), 'DUPLICATE_'+kind)
    previous = {}
    identities = None
    for step, f in enumerate(doc['frames'], 1):
        require(f['step'] == step and set(KINDS) <= set(f), 'STEP_OR_FIELDS')
        keys = {}
        for kind in KINDS:
            try:
                keys[kind] = [key_for(kind, row) for row in f[kind]]
            except (KeyError, TypeError) as exc:
                raise Refusal('MISSING_ID_OR_GENERATION') from exc
            require(len(keys[kind]) == len(set(keys[kind])), 'DUPLICATE_'+kind)
        require(set(keys['nodes']) == set(doc['node_ids']), 'MISSING_NODE_ID')
        require(set(keys['members']) == set(doc['member_ids']), 'MISSING_MEMBER_ID')
        require({k[0] for k in keys['contacts']} == set(doc['joint_ids']), 'MISSING_JOINT_ID')
        expected_order_ids = {jid+'/'+gid for jid, gid, _ in keys['contacts']}
        require(isinstance(doc.get('evaluation_order'), list) and
                len(doc['evaluation_order']) == len(expected_order_ids) and
                set(doc['evaluation_order']) == expected_order_ids, 'EVALUATION_ORDER_IDS')
        if identities is None:
            identities = {k: set(v) for k, v in keys.items()}
        require(identities == {k: set(v) for k, v in keys.items()}, 'MISSING_GP_OR_ENTITY')
        for row in f['contacts']:
            needed = {'normal', 'position_m', 'weight_m2', 'jump_m', 'traction_kpa', 'force_kn',
                      'plastic_slip_m', 'accumulated_slip_m', 'state', 'open', 'gap_m',
                      'elastic_energy_kn_m', 'friction_increment_kn_m', 'work_increment_kn_m',
                      'work_defect_increment_kn_m'}
            require(needed <= set(row), 'MISSING_CONTACT_FIELD')
            require(row['generation'] == 0, 'GENERATION')
            require(row['state'] in ('open', 'stick', 'slip') and
                    row['open'] == (row['state'] == 'open'), 'OPEN_STATE')
            require(row['weight_m2'] > 0 and row['accumulated_slip_m'] >= 0, 'CONTACT_MEASURE')
            n = np.array(row['normal']); q = np.array(row['traction_kpa'])
            require(abs(n@n-1) <= TOLERANCES['independent_elastic_relative'], 'NORMAL')
            # Force/stress/measure check uses the quoted reassembly diagnostic.
            require(np.max(abs(q*row['weight_m2']-row['force_kn'])) <= TOLERANCES['nodal_force_kn'],
                    'FORCE_TRACTION_UNITS_OR_SIGN')
            key = key_for('contacts', row)
            old = previous.get(key, 0.)
            require(row['accumulated_slip_m'] >= old, 'SLIP_HISTORY_DECREASE')
            previous[key] = row['accumulated_slip_m']
        for row in f['nodes']:
            needed = {'position_m', 'displacement_m', 'bulk_internal_kn', 'contact_internal_kn',
                      'structure_internal_kn', 'applied_kn', 'reaction_kn', 'residual_kn'}
            require(needed <= set(row), 'MISSING_NODE_FIELD')
            r = (np.array(row['bulk_internal_kn'])+row['contact_internal_kn']+
                 np.array(row['structure_internal_kn'])-row['applied_kn']-np.array(row['reaction_kn']))
            require(np.max(abs(r-row['residual_kn'])) <= TOLERANCES['nodal_force_kn'], 'NODAL_FORCE_ACCOUNTING')
            require(np.max(abs(r)) <= TOLERANCES['nodal_force_kn'], 'NODAL_EQUILIBRIUM')
        for row in f['members']:
            require({'kind', 'node_ids', 'axial_force_kn', 'end_internal_kn',
                     'end_moment_kn_m', 'elastic_energy_kn_m'} <= set(row), 'MISSING_SECTION_FORCE')
        require({'external_work_kn_m', 'bulk_energy_kn_m', 'contact_energy_kn_m',
                 'member_energy_kn_m', 'friction_dissipation_kn_m', 'diagnostics'} <= set(f), 'MISSING_WORK')
    return doc


def export_csv(doc, path):
    validate(doc)
    path = Path(path)
    manifest = {k: v for k, v in doc.items() if k != 'frames'}
    manifest['frame_headers'] = [{k: v for k, v in f.items() if k not in KINDS} for f in doc['frames']]
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=['step', 'kind', 'identity', 'payload_json'])
        writer.writeheader()
        for f in doc['frames']:
            for kind in KINDS:
                for row in f[kind]:
                    writer.writerow(dict(step=f['step'], kind=kind,
                                         identity=json.dumps(key_for(kind, row), ensure_ascii=False),
                                         payload_json=json.dumps(row, ensure_ascii=False, allow_nan=False)))
    manifest['csv_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    dump(path.with_suffix('.manifest.json'), manifest)


def load(path, manifest=None):
    path = Path(path)
    if path.suffix == '.json':
        require(manifest is None, 'UNEXPECTED_MANIFEST')
        return validate(json.loads(path.read_text(encoding='utf-8-sig')))
    require(path.suffix == '.csv' and manifest is not None, 'CSV_MANIFEST_REQUIRED')
    header = json.loads(Path(manifest).read_text(encoding='utf-8-sig'))
    require(header.pop('csv_sha256') == hashlib.sha256(path.read_bytes()).hexdigest(), 'CSV_HASH')
    frames = header.pop('frame_headers')
    for f in frames:
        for k in KINDS:
            f[k] = []
    with path.open(newline='', encoding='utf-8') as stream:
        reader = csv.DictReader(stream)
        require(reader.fieldnames == ['step', 'kind', 'identity', 'payload_json'], 'CSV_COLUMNS')
        for row in reader:
            step = int(row['step']); kind = row['kind']
            require(1 <= step <= len(frames) and kind in KINDS, 'CSV_STEP_KIND')
            payload = json.loads(row['payload_json'])
            require(json.loads(row['identity']) == json.loads(json.dumps(key_for(kind, payload))), 'CSV_ID')
            frames[step-1][kind].append(payload)
    return validate(dict(header, frames=frames))


def remap(doc, joint_map):
    doc = copy.deepcopy(doc)
    require(set(joint_map) == set(doc['joint_ids']), 'JOINT_MAP_MISSING_SOURCE_ID')
    require(len(set(joint_map.values())) == len(joint_map), 'JOINT_MAP_NOT_BIJECTIVE')
    doc['joint_ids'] = [joint_map[v] for v in doc['joint_ids']]
    # Keep GP evaluation identities coherent with the explicit JOINT mapping.
    old_prefixes = sorted(joint_map, key=len, reverse=True)
    doc['evaluation_order'] = [next(joint_map[jid]+entry[len(jid):]
                                    for jid in old_prefixes if entry.startswith(jid+'/'))
                               for entry in doc['evaluation_order']]
    for f in doc['frames']:
        for r in f['contacts']:
            r['joint_id'] = joint_map[r['joint_id']]
    return doc


def raw_differences(a, b, path='', output=None):
    """All numeric differences retained; no mean/clip/interpolation/nearest GP."""
    output = {} if output is None else output
    if isinstance(a, dict):
        require(isinstance(b, dict) and a.keys() == b.keys(), 'FIELD_SET:'+path)
        for key in a:
            raw_differences(a[key], b[key], path+'/'+key, output)
    elif isinstance(a, list):
        require(isinstance(b, list) and len(a) == len(b), 'SHAPE:'+path)
        for i, (aa, bb) in enumerate(zip(a, b)):
            raw_differences(aa, bb, path+f'/{i}', output)
    elif isinstance(a, (float, int)) and not isinstance(a, bool):
        require(isinstance(b, (float, int)) and not isinstance(b, bool), 'NUMERIC_TYPE:'+path)
        if a != b:
            output[path] = dict(reference=a, candidate=b, absolute=abs(a-b))
    else:
        require(type(a) == type(b) and a == b, 'DISCRETE:'+path)
    return output


def field_differences(a, b):
    require(a['metadata'] == b['metadata'], 'PHYSICS_OR_DISCRETIZATION; no automatic resampling')
    require(set(a['joint_ids']) == set(b['joint_ids']), 'JOINT_MAP_TARGET_MISMATCH')
    require(set(a['node_ids']) == set(b['node_ids']) and set(a['member_ids']) == set(b['member_ids']), 'ENTITY_IDS')
    require(len(a['frames']) == len(b['frames']), 'HISTORY_COUNT')
    differences = {}
    for fa, fb in zip(a['frames'], b['frames']):
        require((fa['step'], fa['phase'], fa['progress']) == (fb['step'], fb['phase'], fb['progress']), 'LOAD_PATH')
        for kind in KINDS:
            ma = {key_for(kind, r): r for r in fa[kind]}
            mb = {key_for(kind, r): r for r in fb[kind]}
            require(ma.keys() == mb.keys(), 'MISSING_GP_OR_ENTITY:'+kind)
            for key in ma:
                raw_differences(ma[key], mb[key], f"/{fa['step']}/{kind}/{key}", differences)
        for key in ('external_work_kn_m', 'bulk_energy_kn_m', 'contact_energy_kn_m',
                    'member_energy_kn_m', 'friction_dissipation_kn_m'):
            raw_differences(fa[key], fb[key], f"/{fa['step']}/{key}", differences)
    return differences


def compare(reference, candidate, joint_map, order_run):
    validate(reference); validate(candidate)
    require(order_run is not None, 'ORDER_WITNESS_REQUIRED; one output cannot reveal evaluation order dependence')
    validate(order_run)
    require(candidate['evaluation_order'] != order_run['evaluation_order'], 'ORDER_WITNESS_NOT_REORDERED')
    c = remap(candidate, joint_map); o = remap(order_run, joint_map)
    differences = field_differences(reference, c)
    order_differences = field_differences(c, o)
    # An exact raw-data contract is intentionally stricter than an engineering
    # tolerance. Undefined slip/opening/section/order criteria cannot be invented.
    # Nonidentical data is returned in full and comparison acceptance is refused.
    return dict(schema='dry-multi-joint-comparison/1', accepted=None,
                validation_complete=True, exact_raw_match=not differences,
                exact_order_match=not order_differences,
                comparison_refused=bool(differences or order_differences),
                reason='RAW_DIFFERENCE_REQUIRES_DEFINED_POLICY' if differences or order_differences else 'diagnostic identity only; A-01 remains unmet',
                quoted_diagnostic_tolerances=TOLERANCES, tolerance_source=CITATION,
                slip_opening_section_order_tolerance='未定義',
                differences=differences, order_differences=order_differences,
                caveat='Undisclosed correction cannot be identified from final values alone; require trusted raw export/provenance and independent reassembly.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--candidate', type=Path, required=True)
    p.add_argument('--manifest', type=Path)
    p.add_argument('--order-run', type=Path)
    p.add_argument('--order-manifest', type=Path)
    p.add_argument('--joint-map', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--export-csv', type=Path)
    a = p.parse_args()
    try:
        reference = load(a.reference)
        candidate = load(a.candidate, a.manifest)
        if a.export_csv:
            export_csv(candidate, a.export_csv)
        result = compare(reference, candidate, json.loads(a.joint_map.read_text(encoding='utf-8-sig')),
                         load(a.order_run, a.order_manifest) if a.order_run else None)
        result['input_files_sha256'] = {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                                        for path in (a.reference, a.candidate, a.joint_map,
                                                     a.manifest, a.order_run, a.order_manifest) if path is not None}
        status = 2 if result['comparison_refused'] else 0
    except (Refusal, KeyError, TypeError, ValueError) as exc:
        result = dict(accepted=None, comparison_refused=True, reason=str(exc)); status = 2
    dump(a.output, result)
    print(json.dumps({k: v for k, v in result.items() if k not in ('differences', 'order_differences')}))
    return status


if __name__ == '__main__':
    raise SystemExit(main())
