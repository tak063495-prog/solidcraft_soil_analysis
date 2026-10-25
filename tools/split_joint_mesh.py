"""Split declared material interfaces into independent displacement nodes.

Coordinates, element geometry, material IDs and CAD volumes remain unchanged.
Node stars stay connected across all undeclared/glued faces; unsupported tied
crack fronts are rejected instead of silently disconnecting an adjacent body.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import numpy as np
import fem_reference as ref


def split(mesh,spec):
    assert spec['schema']=='soil-fem-joint-split/1'
    result=copy.deepcopy(mesh);cells=ref.cells(mesh);owners=ref.topology(mesh)
    regions={r['id']:r['body'] for r in mesh['regions']};bodies=[regions[r] for r in mesh['region_ids']]
    declared={};groups={}
    for interface in spec['interfaces']:
        name=interface['id'];a=interface['body_a'];b=interface['body_b']
        if not name or name in groups or a==b or a not in bodies or b not in bodies:raise ValueError('invalid/duplicate interface or body')
        selected=[]
        for face,os in owners.items():
            if len(os)==2 and {bodies[o] for o in os}=={a,b}:
                if face in declared:raise ValueError('interface face selected twice')
                oa=next(o for o in os if bodies[o]==a);ob=next(o for o in os if bodies[o]==b)
                declared[face]=(oa,ob);selected.append((face,oa,ob))
        if not selected:raise ValueError('no shared volume faces for declared interface')
        groups[name]=selected
    if not groups:raise ValueError('nonempty joint split specification required')
    def face_nodes(face,owner):
        cell=cells[owner]
        local=next(f for f in ref.FACES if set(cell[list(f)])==set(face));nodes=list(cell[list(local)])
        if len(cell)==10:
            for i,j in ((0,1),(1,2),(2,0)):
                e=next(k for k,p in enumerate(ref.EDGES) if set(p)=={local[i],local[j]});nodes.append(cell[e+4])
        return nodes
    affected={n for face,(a,b) in declared.items() for n in face_nodes(face,a)}
    stars={n:[] for n in affected}
    for ei,cell in enumerate(cells):
        for n in cell:
            if n in stars:stars[n].append(ei)
    glued={n:[] for n in affected}
    for face,os in owners.items():
        if len(os)!=2 or face in declared:continue
        for n in face_nodes(face,os[0]):
            if n in glued:glued[n].append(os)
    remap={};coordinates=copy.deepcopy(mesh['nodes'])
    for node,incident in stars.items():
        parent={e:e for e in incident}
        def find(e):
            while parent[e]!=e:parent[e]=parent[parent[e]];e=parent[e]
            return e
        for a,b in glued[node]:parent[find(b)]=find(a)
        components={}
        for e in incident:components.setdefault(find(e),[]).append(e)
        components=sorted(components.values(),key=min)
        for ci,component in enumerate(components):
            new=int(node) if ci==0 else len(coordinates)
            if ci:coordinates.append(mesh['nodes'][node].copy())
            for e in component:remap[e,node]=new
    def mapped(owner,node):return remap.get((owner,node),int(node))
    new_cells=[[mapped(e,n) for n in cell] for e,cell in enumerate(cells)]
    result['nodes']=coordinates;result['tetrahedra']=[c[:4] for c in new_cells]
    if mesh.get('tetrahedra10'):result['tetrahedra10']=new_cells
    face_sets={}
    for name,faces in mesh.get('face_sets',{}).items():
        transformed={tuple(sorted(mapped(o,n) for n in face)) for face in faces for o in owners[tuple(sorted(face))]}
        face_sets[name]=[list(f) for f in sorted(transformed)]
    for name,selected in groups.items():
        for side,index in (('a',1),('b',2)):
            key=f'{name}:{side}'
            if key in face_sets:raise ValueError('joint side name conflicts with existing set')
            face_sets[key]=[sorted(mapped(item[index],n) for n in item[0]) for item in selected]
        for face,a,b in selected:
            if any(mapped(a,n)==mapped(b,n) for n in face_nodes(face,a)):
                raise ValueError('tied crack fronts are unsupported; use a fully separated body interface')
    result['face_sets']=face_sets
    # Preserve geometry-selected monitors on both displacement sides, then add
    # explicit side sets for contact loads, diagnostics and boundary controls.
    for name,nodes in result['node_sets'].items():
        selected=set(nodes);selected.update(new for (e,old),new in remap.items() if old in nodes)
        result['node_sets'][name]=sorted(selected)
    result['node_sets']['all']=list(range(len(coordinates)))
    new_owners=ref.topology(result)
    for name in groups:
        for side in ('a','b'):
            key=f'{name}:{side}'
            result['node_sets'][key]=sorted({int(n) for f in face_sets[key] for n in ref.face_nodes(result,new_owners[tuple(f)][0],tuple(f))})
    assert result['material_ids']==mesh['material_ids'] and result['region_ids']==mesh['region_ids']
    assert np.array_equal(np.asarray(result['nodes'])[np.asarray(new_cells)],np.asarray(mesh['nodes'])[cells])
    return result,{'added_displacement_nodes':len(coordinates)-len(mesh['nodes']),'interfaces':{n:len(f) for n,f in groups.items()}}


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
