"""Add body-owned reference boundary sets to an already generated mesh.

Only named node/face sets are added. Geometry, connectivity and material IDs
remain the original objects. Coordinate planes require every quadratic node.
"""
import numpy as np

EDGES=((0,1),(1,2),(2,0),(0,3),(1,3),(2,3))
FACES=((1,2,3),(0,3,2),(0,1,3),(0,2,1))

def add_body_boundary_sets(mesh):
    xyz=np.asarray(mesh['nodes'])
    cells=mesh.get('tetrahedra10') or mesh['tetrahedra']
    region_ids=mesh['region_ids']
    faces={}
    for e,cell in enumerate(mesh['tetrahedra']):
        for face in FACES:
            key=tuple(sorted(cell[i] for i in face))
            faces.setdefault(key,[]).append(e)
    def face_nodes(face,owner):
        cell=cells[owner]
        for local in FACES:
            if {cell[i] for i in local}!=set(face):continue
            ids=[cell[i] for i in local]
            if mesh.get('tetrahedra10'):
                for a,b in ((0,1),(1,2),(2,0)):
                    edge=next(i for i,pair in enumerate(EDGES) if set(pair)=={local[a],local[b]})
                    ids.append(cell[edge+4])
            return ids
        raise ValueError('body boundary face lost its owner')
    proposed_faces={};proposed_nodes={}
    for region in mesh['regions']:
        boundary=[];exterior=[]
        for face,owners in faces.items():
            own=[e for e in owners if region_ids[e]==region['id']]
            if len(own)!=1:continue
            if len(owners)==2 and region_ids[owners[0]]==region_ids[owners[1]]:continue
            boundary.append((face,own[0]))
            if len(owners)==1:exterior.append((face,own[0]))
        named=[('boundary',boundary),('exterior',exterior)]
        for axis,label in enumerate('xyz'):
            for side,suffix in enumerate(('min','max')):
                coordinate=region['bounds_m'][side][axis]
                selected=[(face,owner) for face,owner in boundary
                    if np.allclose(xyz[face_nodes(face,owner),axis],coordinate,atol=1e-8,rtol=0)]
                if selected:named.append((label+suffix,selected))
        for suffix,selected in named:
            if not selected:continue
            name='body:'+region['body']+':'+suffix
            face_set=[list(face) for face,_ in selected]
            node_set=sorted({node for face,owner in selected for node in face_nodes(face,owner)})
            for mapping,value in ((mesh['face_sets'],face_set),(mesh['node_sets'],node_set)):
                if name in mapping and mapping[name]!=value:
                    raise ValueError('body boundary set collides with an existing set: '+name)
            proposed_faces[name]=face_set
            proposed_nodes[name]=node_set
    mesh['face_sets'].update(proposed_faces)
    mesh['node_sets'].update(proposed_nodes)
