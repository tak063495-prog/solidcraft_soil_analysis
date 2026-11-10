"""Small deterministic conformal meshes for the declared Tet10 crack tip."""
import copy
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import fem_reference as ref
from split_joint_mesh import split


def conformal(x=(0., 1., 2.), y=(0., 1.), z=(0., .5, 1.), same_body=False):
    xyz = [[a,b,c] for c in z for b in y for a in x]
    nx, ny = len(x), len(y)
    def index(i,j,k): return i+nx*(j+ny*k)
    tets, regions = [], []
    pattern = ((0,1,3,7),(0,3,2,7),(0,2,6,7),(0,6,4,7),(0,4,5,7),(0,5,1,7))
    for k in range(len(z)-1):
        for j in range(len(y)-1):
            for i in range(len(x)-1):
                cube = [index(i+a,j+b,k+c) for c in (0,1) for b in (0,1) for a in (0,1)]
                for p in pattern:
                    tets.append([cube[n] for n in p])
                    regions.append(1 if same_body or z[k] < .5 else 2)
    nodes = copy.deepcopy(xyz)
    edges, cells = {}, []
    for tet in tets:
        cell = tet.copy()
        for a,b in ref.EDGES:
            edge = tuple(sorted((tet[a],tet[b])))
            if edge not in edges:
                edges[edge] = len(nodes)
                nodes.append(((np.asarray(nodes[edge[0]])+nodes[edge[1]])/2).tolist())
            cell.append(edges[edge])
        cells.append(cell)
    mesh = dict(schema='soil-fem-mesh/3', units='m-kN-kPa', axis='z_up', crs='local', origin_m=[0.,0.,0.],
                nodes=nodes, tetrahedra=tets, tetrahedra10=cells, quadrature_order=3,
                material_ids=regions.copy(), region_ids=regions, node_sets={}, face_sets={}, regions=[])
    for rid in sorted(set(regions)):
        ei = [e for e,r in enumerate(regions) if r == rid]
        pts = np.asarray(nodes)[np.unique(np.asarray(cells)[ei])]
        volume = sum(abs(np.linalg.det((np.asarray(nodes)[tets[e]][1:]-np.asarray(nodes)[tets[e]][0]).T))/6 for e in ei)
        mesh['regions'].append(dict(id=rid,body='LowerClay' if rid==1 else 'UpperClay',material_id=rid,
            bounds_m=[pts.min(axis=0).tolist(),pts.max(axis=0).tolist()],volume_m3=volume,source_sha256=hashlib.sha256(b'joint-tip-fixture').hexdigest()))
    xyz = np.asarray(nodes)
    for name,a,value in [('bottom',2,z[0]),('top',2,z[-1]),('left',0,x[0]),('right',0,x[-1])]:
        mesh['node_sets'][name] = np.flatnonzero(np.abs(xyz[:,a]-value)<1e-12).tolist()
    mesh['node_sets']['all'] = list(range(len(nodes)))
    owners = ref.topology(mesh)
    mesh['face_sets']['exterior'] = [list(f) for f,o in owners.items() if len(o)==1]
    mesh['face_sets']['top'] = [list(f) for f,o in owners.items() if len(o)==1 and np.all(xyz[list(f),2]==z[-1])]
    mesh['face_sets']['bottom'] = [list(f) for f,o in owners.items() if len(o)==1 and np.all(xyz[list(f),2]==z[0])]
    return mesh


def specification(mesh, full=False):
    xyz = np.asarray(mesh['nodes'])
    selected = []
    for face,os in ref.topology(mesh).items():
        p = xyz[list(face)]
        if len(os)==2 and np.all(np.abs(p[:,2]-.5)<1e-12) and (full or p[:,0].max()<=1):
            a,b = sorted(os,key=lambda e:xyz[mesh['tetrahedra'][e],2].mean())
            selected.append(dict(nodes=list(face),owner_a=a,owner_b=b))
    tip_x = 0. if full else 1.
    edges = {}
    for item in selected:
        ns = ref.face_nodes(mesh,item['owner_a'],item['nodes'])
        for k,(i,j) in enumerate(((0,1),(1,2),(2,0))):
            if xyz[ns[i],0]==tip_x and xyz[ns[j],0]==tip_x:
                edge = tuple(sorted((int(ns[i]),int(ns[j]))))
                edges[edge] = int(ns[3+k])
    return dict(schema='soil-fem-joint-split/1',interfaces=[dict(id='Seam',faces=selected,tip_edges=[[v,w,m] for (v,w),m in sorted(edges.items())])])


def fixture(full=False, same_body=False, **kwargs):
    source = conformal(same_body=same_body, **kwargs)
    spec = specification(source, full)
    mesh, summary = split(source, spec)
    configuration = dict(schema='soil-fem-joints/1',geometric_tolerance_m=1e-10,maximum_increment_m=5e-4,tangent_step_m=1e-8,
        joints=[dict(id='Seam',side_a='Seam:a',side_b='Seam:b',material=dict(normal_stiffness_kpa_m=1e6,tangential_stiffness_kpa_m=2e5,friction_coefficient=.3),
                     crack_tip=summary['crack_tip_declarations']['Seam'])])
    return source, spec, mesh, configuration
