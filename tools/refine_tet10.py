"""Uniform nested Tet10 subdivision preserving the isoparametric geometry.

Each parent yields eight quadratic children. Node identity uses rational parent
corner topology, never coordinate merging: coincident split/contact traces stay
split. Material/region/excavation ancestry and named face traces are inherited.
This refines an input mesh, not solved material histories or initial GP seeds.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import numpy as np
import fem_reference as ref


CHILDREN = ((0, 4, 6, 7), (4, 1, 5, 8), (6, 5, 2, 9), (7, 8, 9, 3),
            (4, 9, 5, 6), (4, 9, 6, 7), (4, 9, 7, 8), (4, 9, 8, 5))
BARYCENTRIC = np.asarray(((4, 0, 0, 0), (0, 4, 0, 0), (0, 0, 4, 0), (0, 0, 0, 4),
                          (2, 2, 0, 0), (0, 2, 2, 0), (2, 0, 2, 0),
                          (2, 0, 0, 2), (0, 2, 0, 2), (0, 0, 2, 2)))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def refine(mesh):
    if mesh.get('schema') != 'soil-fem-mesh/3' or not mesh.get('tetrahedra10'):
        raise ValueError('quadratic soil-fem-mesh/3 required')
    old_xyz = np.asarray(mesh['nodes'])
    old_cells = np.asarray(mesh['tetrahedra10'])
    if old_cells.shape != (len(mesh['tetrahedra']), 10):
        raise ValueError('Tet10 inventory differs from corner cells')
    xyz = old_xyz.tolist()
    keys, reverse = {}, {}

    def key(corners, weights):
        return tuple(sorted((int(node), int(weight)) for node, weight in zip(corners, weights) if weight))

    for cell in old_cells:
        for node, weights in zip(cell, BARYCENTRIC):
            label = key(cell[:4], weights)
            if label in keys and keys[label] != node:
                raise ValueError('coarse edge has inconsistent quadratic node identity')
            if node in reverse and reverse[node] != label:
                raise ValueError('coarse mesh mixes corner/midpoint node roles')
            keys[label] = int(node)
            reverse[int(node)] = label
    topology = ref.topology(mesh)
    new_cells, parents, origins = [], [], []
    maximum_shared_geometry_error = 0.

    def node_id(cell, weights):
        nonlocal maximum_shared_geometry_error
        label = key(cell[:4], weights)
        natural = np.asarray(weights[1:], dtype=float)/4
        shape, _ = ref.basis(natural[None, :], True)
        position = shape[0]@old_xyz[cell]
        if label not in keys:
            keys[label] = len(xyz)
            xyz.append(position.tolist())
        else:
            maximum_shared_geometry_error = max(maximum_shared_geometry_error,
                                                float(abs(np.asarray(xyz[keys[label]])-position).max()))
        return keys[label]

    for parent, cell in enumerate(old_cells):
        for template in CHILDREN:
            weights = BARYCENTRIC[list(template)].copy()
            natural = weights[:, 1:]/4
            determinant = float(np.linalg.det((natural[1:]-natural[0]).T))
            assert abs(abs(determinant)-.125) < 1e-14
            if determinant < 0:
                weights[[1, 2]] = weights[[2, 1]]
                natural[[1, 2]] = natural[[2, 1]]
            nodes = [node_id(cell, w) for w in weights]
            for a, b in ref.EDGES:
                total = weights[a]+weights[b]
                assert not np.any(total % 2)
                nodes.append(node_id(cell, total//2))
            new_cells.append(nodes)
            parents.append(parent)
            origins.append(natural)
    assert maximum_shared_geometry_error < 1e-10
    result = copy.deepcopy(mesh)
    result['nodes'] = xyz
    result['tetrahedra10'] = new_cells
    result['tetrahedra'] = [c[:4] for c in new_cells]
    for name in ('material_ids', 'region_ids'):
        result[name] = [mesh[name][e] for e in parents]
    result['excavation_sets'] = {name: [8*e+i for e in selected for i in range(8)]
                                 for name, selected in mesh.get('excavation_sets', {}).items()}
    new_topology = ref.topology(result)
    face_sets, original_traces = {}, {}
    for name, faces in mesh.get('face_sets', {}).items():
        children, trace = [], set()
        for face in faces:
            owner = topology[tuple(sorted(face))][0]
            nodes = ref.face_nodes(mesh, owner, face).tolist()
            assert len(nodes) == 6
            a, b, c, ab, bc, ca = nodes
            children.extend(((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)))
            trace.update(nodes)
        assert all(tuple(sorted(face)) in new_topology for face in children)
        face_sets[name] = [list(f) for f in children]
        original_traces[name] = trace
    result['face_sets'] = face_sets
    policies = {}
    for name, nodes in mesh.get('node_sets', {}).items():
        if set(nodes) == set(range(len(old_xyz))):
            result['node_sets'][name] = list(range(len(xyz)))
            policies[name] = 'entire original node inventory -> entire refined inventory'
        elif name in original_traces and set(nodes) == original_traces[name]:
            trace = set()
            for face in face_sets[name]:
                owner = new_topology[tuple(sorted(face))][0]
                trace.update(ref.face_nodes(result, owner, face).tolist())
            result['node_sets'][name] = sorted(trace)
            policies[name] = 'exact named face trace -> inherited refined face trace'
        else:
            policies[name] = 'literal original node IDs preserved; no new geometric selection inferred'
    # A quadratic parent map composed with each affine child map is represented
    # exactly by the child Tet10 geometry, also for curved parent edges/faces.
    probes = np.asarray(((.11, .17, .23), (.4, .1, .2), (.07, .31, .29)))
    child_shape, _ = ref.basis(probes, True)
    maximum_geometry_error = 0.
    new_xyz = np.asarray(xyz)
    for cell, parent, natural in zip(new_cells, parents, origins):
        parent_points = natural[0]+probes@(natural[1:]-natural[0])
        parent_shape, _ = ref.basis(parent_points, True)
        expected = parent_shape@old_xyz[old_cells[parent]]
        actual = child_shape@new_xyz[cell]
        maximum_geometry_error = max(maximum_geometry_error, float(abs(actual-expected).max()))
    assert maximum_geometry_error < 1e-10
    old_volumes = np.asarray([ref.volume(mesh, e)['weight'].sum() for e in range(len(old_cells))])
    volumes = np.asarray([ref.volume(result, e)['weight'].sum() for e in range(len(new_cells))])
    volume_error = float(abs(volumes.reshape(-1, 8).sum(axis=1)-old_volumes).max())
    assert volume_error < 1e-9
    audit = dict(schema='soil-fem-tet10-nested-refinement/1', refined=True, parent_elements=len(old_cells),
                 child_elements=len(new_cells), parent_nodes=len(old_xyz), child_nodes=len(xyz),
                 parent_of_element=parents, node_set_policies=policies,
                 maximum_parent_geometry_reproduction_error_m=maximum_geometry_error,
                 maximum_shared_geometry_error_m=maximum_shared_geometry_error,
                 maximum_parent_volume_conservation_error_m3=volume_error,
                 original_volume_m3=float(old_volumes.sum()), refined_volume_m3=float(volumes.sum()),
                 excavation_volumes_m3={name: float(volumes[ids].sum()) for name, ids in result['excavation_sets'].items()},
                 goal_complete=False, p6_complete=False, spatial_convergence_verified=False,
                 limitations='Mesh geometry/topology only; regenerate initial GP seeds; no solved history transfer or arbitrary geometric interpretation of literal node sets')
    return result, audit


def run(source, output):
    audit_path = output.with_suffix('.refinement.json')
    if output.exists() or audit_path.exists():
        raise ValueError('refined mesh or audit output already exists')
    digest = sha(source)
    mesh = json.loads(source.read_text(encoding='utf-8-sig'))
    refined, audit = refine(mesh)
    assert sha(source) == digest
    with output.open('x', encoding='utf-8') as stream:
        json.dump(refined, stream, allow_nan=False)
    audit['sha256'] = {str(source): digest, str(output): sha(output), str(Path(__file__)): sha(Path(__file__)), str(Path(ref.__file__)): sha(Path(ref.__file__))}
    with audit_path.open('x', encoding='utf-8') as stream:
        json.dump(audit, stream, indent=2, allow_nan=False)
    print(json.dumps({k: audit[k] for k in ('parent_elements', 'child_elements', 'child_nodes', 'maximum_parent_geometry_reproduction_error_m', 'maximum_parent_volume_conservation_error_m3', 'excavation_volumes_m3')}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('source', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    run(args.source.resolve(), args.output.resolve())
