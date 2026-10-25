"""Nested Tet10 geometry, C0 field, trace, material and split-contact checks."""
import argparse
import copy
import json
from pathlib import Path
import subprocess
import numpy as np
import fem_reference as ref
import refine_tet10 as nested
import freeze_sources


def write(path, value):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, allow_nan=False)


def two_tetrahedra(split):
    xyz, cells, lookup = [], [], {}
    vertices = (np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]),
                np.array([[0., 0., 0.], [0., 1., 0.], [1., 0., 0.], [0., 0., -1.]]))
    for body, p in enumerate(vertices):
        nodes = []
        positions = list(p)+[(p[a]+p[b])/2 for a, b in ref.EDGES]
        for position in positions:
            key = (body if split else 0, *position)
            if key not in lookup:
                lookup[key] = len(xyz)
                xyz.append(position.tolist())
            nodes.append(lookup[key])
        cells.append(nodes)
    mesh = dict(schema='soil-fem-mesh/3', units='m-kN-kPa', axis='z_up', crs='SYNTHETIC_TEST', origin_m=[0., 0., 0.],
                nodes=xyz, tetrahedra=[c[:4] for c in cells], tetrahedra10=cells,
                quadrature_order=3, material_ids=[1, 2], region_ids=[1, 2], node_sets={'all': list(range(len(xyz)))},
                regions=[dict(id=i+1, body=name, material_id=i+1, bounds_m=[p.min(axis=0).tolist(), p.max(axis=0).tolist()],
                              volume_m3=1/6, source_sha256=str(i)*64) for i, (name, p) in enumerate(zip(('upper', 'lower'), vertices))],
                excavation_sets={'cut_upper': [0], 'cut_lower': [1]},
                face_sets={'upper_trace': cells[0][:3], 'lower_trace': cells[1][:3]})
    mesh['face_sets'] = {k: [v] for k, v in mesh['face_sets'].items()}
    for name, e in (('upper_trace', 0), ('lower_trace', 1)):
        mesh['node_sets'][name] = sorted(ref.face_nodes(mesh, e, mesh['face_sets'][name][0]).tolist())
    return mesh


def field_patch(old, new, audit):
    rng = np.random.default_rng(20261026)
    values = rng.normal(size=(len(old['nodes']), 3))
    transferred = np.full((len(new['nodes']), 3), np.nan)
    probes = np.array([[.13, .21, .27], [.05, .07, .73], [.24, .23, .22]])
    shapes, _ = ref.basis(probes, True)
    maximum_error, maximum_trace_error = 0., 0.
    for child, parent in zip(new['tetrahedra10'], audit['parent_of_element']):
        coarse = old['tetrahedra10'][parent]
        positions = [coarse.index(node) for node in child[:4]]
        natural = nested.BARYCENTRIC[positions, 1:]/4
        child_node_positions = np.vstack((natural, [(natural[a]+natural[b])/2 for a, b in ref.EDGES]))
        prolongation, _ = ref.basis(child_node_positions, True)
        local_values = prolongation@values[coarse]
        for node, value in zip(child, local_values):
            if np.isfinite(transferred[node]).all():
                maximum_trace_error = max(maximum_trace_error, float(abs(transferred[node]-value).max()))
            transferred[node] = value
        mapped = natural[0]+probes@(natural[1:]-natural[0])
        parent_shapes, _ = ref.basis(mapped, True)
        maximum_error = max(maximum_error, float(abs(shapes@local_values-parent_shapes@values[coarse]).max()))
    assert np.isfinite(transferred).all()
    assert max(maximum_error, maximum_trace_error) < 1e-12
    assert np.array_equal(transferred[:len(values)], values)
    return dict(maximum_quadratic_field_reproduction_error=maximum_error,
                maximum_shared_trace_field_error=maximum_trace_error)


def verify(source, solver, output):
    output.mkdir(parents=True, exist_ok=False)
    freeze_sources.freeze(output)
    bindings = {str(p): nested.sha(p) for p in (source, solver)}
    original = json.loads(source.read_text(encoding='utf-8-sig'))
    curved = copy.deepcopy(original)
    xyz = np.array(curved['nodes'])
    xyz[:, 2] += .003*(xyz[:, 0]-6)**2+.002*xyz[:, 0]*xyz[:, 1]
    curved['nodes'] = xyz.tolist()
    curved['crs'] = 'SYNTHETIC_QUADRATIC_SHEAR_TEST'
    for region in curved['regions']:
        nodes = np.unique(np.asarray(curved['tetrahedra10'])[np.asarray(curved['region_ids']) == region['id']])
        region['bounds_m'] = [xyz[nodes].min(axis=0).tolist(), xyz[nodes].max(axis=0).tolist()]
        region['source_sha256'] = '0'*64
    reports = []
    for label, mesh in (('actual_solidcraft', original), ('manufactured_curved', curved),
                        ('two_material_shared', two_tetrahedra(False)), ('split_contact_trace', two_tetrahedra(True))):
        folder = output/label
        folder.mkdir()
        write(folder/'original.json', mesh)
        refined, audit = nested.refine(mesh)
        write(folder/'refined.json', refined)
        write(folder/'refinement.json', audit)
        assert all(refined['material_ids'][e] == mesh['material_ids'][parent] and
                   refined['region_ids'][e] == mesh['region_ids'][parent]
                   for e, parent in enumerate(audit['parent_of_element']))
        for name, ids in mesh.get('excavation_sets', {}).items():
            assert refined['excavation_sets'][name] == [8*e+i for e in ids for i in range(8)]
        assert all(len(owners) in (1, 2) for owners in ref.topology(refined).values())
        patch = field_patch(mesh, refined, audit)
        trace = None
        if label in ('two_material_shared', 'split_contact_trace'):
            a, b = [set(refined['node_sets'][name]) for name in ('upper_trace', 'lower_trace')]
            assert len(a) == len(b) == 15
            if label == 'split_contact_trace':
                assert not a & b
                assert sorted(map(tuple, np.array(refined['nodes'])[sorted(a)])) == sorted(map(tuple, np.array(refined['nodes'])[sorted(b)]))
            else:
                assert a == b
            trace = dict(upper_nodes=len(a), lower_nodes=len(b), shared_node_ids=len(a & b))
        check = subprocess.run([str(solver), 'validate', str(folder/'refined.json')], capture_output=True, text=True, encoding='utf-8')
        (folder/'rust-validation.log').write_text(check.stdout+check.stderr, encoding='utf-8')
        assert check.returncode == 0, check.stderr
        reports.append(dict(case=label, parent_elements=audit['parent_elements'], child_elements=audit['child_elements'],
                            nodes=audit['child_nodes'], geometry_and_volume=audit,
                            field_patch=patch, interface_trace=trace, rust_validation_exit_code=check.returncode))
        print(json.dumps(dict(case=label, children=audit['child_elements'], patch=patch, rust_validation_exit_code=0)), flush=True)
    assert bindings == {str(p): nested.sha(p) for p in (source, solver)}
    assert json.loads((output/'source-snapshot.json').read_text()) == freeze_sources.manifest()
    write(output/'verification.json', dict(schema='soil-fem-nested-tet10-verification/1', passed=True,
          cases=reports, sha256=bindings, goal_complete=False, p6_complete=False, spatial_convergence_verified=False,
          limitations='Geometry/topology/field prolongation only; no material history transfer or nonlinear convergence acceptance'))
    print('Nested Tet10 refinement checks passed', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ('source', 'solver', 'output'):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    verify(*(getattr(args, k).resolve() for k in ('source', 'solver', 'output')))
