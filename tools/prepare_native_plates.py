"""Expand compact CAD plate declarations after meshing, with reproducible linkage.

The captured CAD attributes remain unchanged. Generated node/member records live
in a separate prepared geometry inventory and are rederived before native export.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import numpy as np

import fem_reference
import plate_surfaces
from import_solidcraft_capture import equal, read
from prepare_model_record import record


KEY = 'plate_surfaces'
RESERVED = {'plate-expansion.json', 'plate-surface-binding.json', 'unexpanded-stages.json'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fingerprint(value):
    data = json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode('utf-8')
    return hashlib.sha256(data).hexdigest()


def derive(mesh, captured_case):
    if KEY not in captured_case:
        return copy.deepcopy(mesh), copy.deepcopy(captured_case), None
    specification = captured_case[KEY]
    if not isinstance(specification, dict):
        raise ValueError('embedded plate_surfaces must be a specification object')
    if captured_case.get('joint_system') is not None or captured_case.get('anchor_system') is not None:
        raise ValueError('embedded plate surfaces with contact/anchor coupling require separate verification')
    case = copy.deepcopy(captured_case)
    del case[KEY]
    result_mesh, result_case, binding = plate_surfaces.attach(mesh, case, specification)
    # These singleton sets are only an attachment serialization detail. Use the
    # existing coordinate contract after independently checking the same unique
    # soil node at the solver's tolerance; retain every user-authored node set.
    xyz = np.asarray(mesh['nodes'], dtype=float)
    structural_nodes = {node['id']: node for node in result_case['beam_system']['nodes']}
    for soil in binding['generated_soil_attachment_nodes']:
        name = f'plate:mesh-node:{soil}'
        node = structural_nodes[binding['soil_node_to_structure_id'][str(soil)]]
        matches = np.flatnonzero(np.linalg.norm(xyz - node['position_m'], axis=1) <= 1e-8)
        if matches.tolist() != [soil]:
            raise ValueError('generated native plate coordinate attachment is ambiguous')
        if node.get('soil_node_set') != name or result_mesh['node_sets'].get(name) != [soil]:
            raise ValueError('generated native plate attachment differs from its soil-node binding')
        del node['soil_node_set']
        node['match_soil_node'] = True
        del result_mesh['node_sets'][name]
    if not equal({k: v for k, v in mesh.items() if k != 'node_sets'},
                 {k: v for k, v in result_mesh.items() if k != 'node_sets'}):
        raise ValueError('plate expansion changed soil geometry or connectivity')
    return result_mesh, result_case, binding


def receipt(source, mesh, case, binding):
    return dict(schema='soil-fem-native-plate-expansion/1',
                source_sha256={name: sha(source/name) for name in ('mesh.json', 'stages.json', 'native-capture.json')},
                implementation_sha256={Path(path).name: sha(Path(path)) for path in
                                       (__file__, plate_surfaces.__file__, fem_reference.__file__)},
                generated_fingerprints={name: fingerprint(value) for name, value in
                                        [('mesh', mesh), ('case', case), ('binding', binding)]},
                geometry_connectivity_unchanged=True,
                node_binding='Unique coordinate attachments after meshing; exact shared soil-node translations',
                coordinate_matching_tolerance_m=1e-8,
                generated_singleton_sets_serialized=False,
                persistent_cad_face_identity_claimed=False, p6_complete=False)


def prepare(source, output):
    source, output = source.resolve(), output.resolve()
    if output.exists() or output.is_relative_to(source):
        raise ValueError('plate expansion requires a new folder outside its source inventory')
    inventory = read(source/'prepared-model.json')
    if inventory.get('schema') != 'soil-fem-prepared-model/1' or inventory.get('equilibrium_verified') is not False:
        raise ValueError('geometry-only prepared inventory required')
    files = inventory.get('files_sha256', {})
    if not {'mesh.json', 'stages.json', 'native-capture.json'} <= set(files) or RESERVED & set(files):
        raise ValueError('missing source files or already expanded plate inventory')
    for name, digest in files.items():
        path = source/name
        if Path(name).name != name or '/' in name or '\\' in name or ':' in name or not path.resolve().is_relative_to(source):
            raise ValueError('prepared inventory must contain safe plain filenames')
        if sha(path) != digest:
            raise ValueError('prepared plate input changed: '+name)
    original_mesh = read(source/'mesh.json', 256*1024**2)
    original_case = read(source/'stages.json')
    if KEY not in original_case:
        raise ValueError('compact native plate declaration required')
    mesh, case, binding = derive(original_mesh, original_case)
    bound = receipt(source, mesh, case, binding)
    # Validate everything before creating an output or copying an immutable file.
    if any(sha(source/name) != digest for name, digest in files.items()):
        raise ValueError('plate source inventory changed during expansion')
    output.mkdir(parents=True)
    for name, digest in files.items():
        if name in ('mesh.json', 'stages.json'):
            continue
        with (output/name).open('xb') as stream:
            stream.write((source/name).read_bytes())
        if sha(output/name) != digest:
            raise ValueError('plate input changed during copy')
    for name, value in [('mesh.json', mesh), ('stages.json', case),
                        ('unexpanded-stages.json', original_case), ('plate-surface-binding.json', binding),
                        ('plate-expansion.json', bound)]:
        with (output/name).open('x', encoding='utf-8') as stream:
            json.dump(value, stream, indent=2, allow_nan=False)
    record(output)
    verify(source, output, mesh, case, original_case)
    return bound


def verify(conformal, prepared, mesh, case, captured_case):
    if KEY not in captured_case:
        if not equal(case, captured_case) or (prepared/'plate-expansion.json').exists():
            raise ValueError('job stages differ from native capture')
        return None
    if not equal(read(conformal/'stages.json'), captured_case):
        raise ValueError('compact plate declaration differs from its captured source')
    original_mesh = read(conformal/'mesh.json', 256*1024**2)
    expected_mesh, expected_case, expected_binding = derive(original_mesh, captured_case)
    if not equal(mesh, expected_mesh) or not equal(case, expected_case):
        raise ValueError('native plate mesh/stages differ from regenerated CAD declarations')
    if not equal(read(prepared/'unexpanded-stages.json'), captured_case):
        raise ValueError('native compact plate declaration changed')
    if not equal(read(prepared/'plate-surface-binding.json'), expected_binding):
        raise ValueError('native plate owner/trace binding differs')
    expected_receipt = receipt(conformal, expected_mesh, expected_case, expected_binding)
    if not equal(read(prepared/'plate-expansion.json'), expected_receipt):
        raise ValueError('native plate expansion source/implementation linkage differs')
    return expected_receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    value = prepare(args.source, args.output)
    print(json.dumps(dict(expanded=True,generated_fingerprints=value['generated_fingerprints'])))
