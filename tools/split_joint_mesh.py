"""Split declared interfaces, retaining only explicitly declared Tet10 tips.

Legacy body-pair selection remains fully separated. Explicit faces name their
corner IDs and both owning element indices, including within a single body.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import numpy as np
import fem_reference as ref


def split(mesh, spec):
    if spec['schema'] != 'soil-fem-joint-split/1':
        raise ValueError('invalid joint split schema')
    result = copy.deepcopy(mesh)
    cells = ref.cells(mesh)
    owners = ref.topology(mesh)
    if any(len(os) > 2 for os in owners.values()):
        raise ValueError('nonmanifold volume face')
    regions = {r['id']: r['body'] for r in mesh['regions']}
    bodies = [regions[r] for r in mesh['region_ids']]

    def face_nodes(face, owner):
        return [int(n) for n in ref.face_nodes(mesh, owner, face)]

    # Check the combinatorial edge interpolation before modifying any IDs.
    edges = {}
    midpoints = {}
    corners = set(cells[:, :4].flat)
    if cells.shape[1] == 10:
        for cell in cells:
            if len(set(cell)) != 10:
                raise ValueError('invalid Tet10 connectivity')
            for k, (i, j) in enumerate(ref.EDGES):
                edge = tuple(sorted((int(cell[i]), int(cell[j]))))
                mid = int(cell[k+4])
                if (mid in corners or edges.get(edge, mid) != mid
                        or midpoints.get(mid, edge) != edge):
                    raise ValueError('wrong/nonconforming Tet10 intermediate node')
                edges[edge] = mid
                midpoints[mid] = edge

    declared, groups, tips = {}, {}, {}
    for interface in spec['interfaces']:
        name = interface['id']
        if not name or name in groups:
            raise ValueError('invalid/duplicate interface')
        selected = []
        explicit = 'faces' in interface
        if explicit:
            if not interface['faces']:
                raise ValueError('empty explicit face selection')
            for item in interface['faces']:
                face = tuple(sorted(item['nodes']))
                a, b = item['owner_a'], item['owner_b']
                if (len(face) != 3 or len(set(face)) != 3 or a == b
                        or set(owners.get(face, [])) != {a, b}):
                    raise ValueError('explicit face/owner pair mismatch')
                for side, owner in (('a', a), ('b', b)):
                    body = interface.get('body_'+side)
                    if body is not None and body != bodies[owner]:
                        raise ValueError('explicit owner/body mismatch')
                selected.append((face, a, b))
        else:
            a, b = interface['body_a'], interface['body_b']
            if a == b or a not in bodies or b not in bodies:
                raise ValueError('invalid interface body')
            for face, os in owners.items():
                if len(os) == 2 and {bodies[o] for o in os} == {a, b}:
                    oa = next(o for o in os if bodies[o] == a)
                    ob = next(o for o in os if bodies[o] == b)
                    selected.append((face, oa, ob))
        if not selected:
            raise ValueError('no shared volume faces for declared interface')
        surface_edges = {}
        for face, a, b in selected:
            if face in declared:
                raise ValueError('interface face selected twice')
            if set(face_nodes(face, a)) != set(face_nodes(face, b)):
                raise ValueError('nonconforming opposing face interpolation')
            declared[face] = (a, b)
            for i, j in ((0, 1), (1, 2), (2, 0)):
                edge = tuple(sorted((face[i], face[j])))
                surface_edges[edge] = surface_edges.get(edge, 0) + 1
        if any(count > 2 for count in surface_edges.values()):
            raise ValueError('nonmanifold/branching crack surface')
        tip_edges, degree = {}, {}
        for entry in interface.get('tip_edges', []):
            if not explicit or cells.shape[1] != 10 or len(entry) != 3:
                raise ValueError('tip requires explicit faces and Tet10 edges')
            v, w, mid = entry
            edge = tuple(sorted((v, w)))
            if (v == w or edges.get(edge) != mid or edge in tip_edges
                    or surface_edges.get(edge) != 1):
                raise ValueError('invalid/duplicate tip edge or intermediate node')
            tip_edges[edge] = int(mid)
            for n in edge:
                degree[n] = degree.get(n, 0) + 1
            # A volume edge star must be one fan, never disconnected or branched.
            star = {e for e, cell in enumerate(cells) if set(edge) <= set(cell[:4])}
            links = {e: set() for e in star}
            for face, os in owners.items():
                if set(edge) <= set(face) and len(os) == 2:
                    x, y = os
                    links[x].add(y)
                    links[y].add(x)
            pending, visited = [min(star)], set()
            while pending:
                e = pending.pop()
                if e not in visited:
                    visited.add(e)
                    pending.extend(links[e])
            if visited != star or any(len(s) > 2 for s in links.values()):
                raise ValueError('nonmanifold tip edge star')
        if any(n > 2 for n in degree.values()):
            raise ValueError('branching crack tip')
        groups[name] = selected
        tips[name] = tip_edges
    if not groups:
        raise ValueError('nonempty joint split specification required')
    network_edges, network_degree = set(), {}
    for es in tips.values():
        for edge in es:
            if edge in network_edges:
                raise ValueError('duplicate crack-tip edge across interfaces')
            network_edges.add(edge)
            for node in edge:
                network_degree[node] = network_degree.get(node, 0) + 1
    if any(degree > 2 for degree in network_degree.values()):
        raise ValueError('branching crack tip across interfaces')
    tip_nodes = {n for es in tips.values() for (v, w), mid in es.items() for n in (v, w, mid)}
    affected = {n for face, (a, b) in declared.items() for n in face_nodes(face, a)}
    stars = {n: [] for n in affected}
    for ei, cell in enumerate(cells):
        for n in cell:
            if n in stars:
                stars[n].append(ei)
    glued = {n: [] for n in affected}
    for face, os in owners.items():
        if len(os) != 2 or face in declared:
            continue
        for n in face_nodes(face, os[0]):
            if n in glued:
                glued[n].append(os)
    remap, coordinates = {}, copy.deepcopy(mesh['nodes'])
    for node, incident in stars.items():
        parent = {e: e for e in incident}

        def find(e):
            while parent[e] != e:
                parent[e] = parent[parent[e]]
                e = parent[e]
            return e

        for a, b in glued[node]:
            parent[find(b)] = find(a)
        if node in tip_nodes:
            for e in incident:
                parent[find(e)] = find(incident[0])
        components = {}
        for e in incident:
            components.setdefault(find(e), []).append(e)
        for ci, component in enumerate(sorted(components.values(), key=min)):
            new = int(node) if ci == 0 else len(coordinates)
            if ci:
                coordinates.append(mesh['nodes'][node].copy())
            for e in component:
                remap[e, node] = new

    def mapped(owner, node):
        return remap.get((owner, node), int(node))

    new_cells = [[mapped(e, n) for n in cell] for e, cell in enumerate(cells)]
    result['nodes'] = coordinates
    result['tetrahedra'] = [c[:4] for c in new_cells]
    if mesh.get('tetrahedra10'):
        result['tetrahedra10'] = new_cells
    face_sets, declarations = {}, {}
    for name, faces in mesh.get('face_sets', {}).items():
        transformed = {tuple(sorted(mapped(o, n) for n in face))
                       for face in faces for o in owners[tuple(sorted(face))]}
        face_sets[name] = [list(f) for f in sorted(transformed)]
    for name, selected in groups.items():
        allowed = {n for (v, w), mid in tips[name].items() for n in (v, w, mid)}
        for side, index in (('a', 1), ('b', 2)):
            key = f'{name}:{side}'
            if key in face_sets:
                raise ValueError('joint side name conflicts with existing set')
            face_sets[key] = [sorted(mapped(item[index], n) for n in item[0]) for item in selected]
        faces = []
        for face, a, b in selected:
            nodes = face_nodes(face, a)
            shared = {n for n in nodes if mapped(a, n) == mapped(b, n)}
            if shared != set(nodes) & allowed:
                raise ValueError('tied crack fronts require explicit legal tip edges; undeclared shared node')
            if len(shared) == len(nodes):
                raise ValueError('fully tied face has no jump degrees of freedom')
            faces.append(dict(face_a=sorted(mapped(a, n) for n in face),
                              face_b=sorted(mapped(b, n) for n in face), owner_a=a, owner_b=b))
        if tips[name]:
            # Tip IDs stay unchanged; bind the complete surface and owner sides.
            declarations[name] = dict(edges=[[v, w, mid] for (v, w), mid in sorted(tips[name].items())], faces=faces)
    result['face_sets'] = face_sets
    for name, nodes in result['node_sets'].items():
        selected = set(nodes)
        selected.update(new for (e, old), new in remap.items() if old in nodes)
        result['node_sets'][name] = sorted(selected)
    result['node_sets']['all'] = list(range(len(coordinates)))
    new_owners = ref.topology(result)
    for name in groups:
        for side in ('a', 'b'):
            key = f'{name}:{side}'
            result['node_sets'][key] = sorted({int(n) for f in face_sets[key]
                for n in ref.face_nodes(result, new_owners[tuple(f)][0], tuple(f))})
    assert result['material_ids'] == mesh['material_ids'] and result['region_ids'] == mesh['region_ids']
    assert np.array_equal(np.asarray(result['nodes'])[np.asarray(new_cells)], np.asarray(mesh['nodes'])[cells])
    summary = {'added_displacement_nodes': len(coordinates)-len(mesh['nodes']),
               'interfaces': {n: len(f) for n, f in groups.items()}}
    if declarations:
        summary['crack_tip_declarations'] = declarations
    return result, summary

def prepare(source,spec_path,output):
    source=source.resolve();data=json.loads((source/'prepared-model.json').read_text(encoding='utf-8'))
    assert data['schema']=='soil-fem-prepared-model/1' and not data['equilibrium_verified']
    output.mkdir(parents=True,exist_ok=False)
    for name,expected in data['files_sha256'].items():
        if Path(name).name!=name:raise ValueError('invalid prepared filename')
        p=source/name;payload=p.read_bytes();assert hashlib.sha256(payload).hexdigest()==expected
        with (output/name).open('xb') as f:f.write(payload)
    original=output/'mesh.json';payload=original.read_bytes()
    with (output/'mesh-conformal.json').open('xb') as f:f.write(payload)
    mesh,summary=split(json.loads(payload),json.loads(spec_path.read_text(encoding='utf-8')))
    original.write_text(json.dumps(mesh,allow_nan=False),encoding='utf-8')
    provenance=json.loads((output/'mesh.provenance.json').read_text(encoding='utf-8'))
    provenance['joint_split']={**summary,'source_mesh_sha256':hashlib.sha256(payload).hexdigest(),'spec_sha256':hashlib.sha256(spec_path.read_bytes()).hexdigest()}
    (output/'mesh.provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    with (output/'joint-split-spec.json').open('xb') as f:f.write(spec_path.read_bytes())
    record={'schema':'soil-fem-prepared-model/1','equilibrium_verified':False,
        'files_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in output.iterdir() if p.is_file()}}
    with (output/'prepared-model.json').open('x',encoding='utf-8') as f:json.dump(record,f,indent=2)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--prepared',type=Path,required=True);p.add_argument('--spec',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();prepare(a.prepared,a.spec,a.output)
