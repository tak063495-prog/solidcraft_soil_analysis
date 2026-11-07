"""Prepare only: pinned Griffiths geometry, fixed-depth in-plane Tet10 h series.

No 2D displacements are inferred. No slope solve runs on import or preparation.
Default material is unreduced F=1. An explicit pinned reference file prepares
a separate series at the archived fixed reduction factor, never an Fs claim.
"""
import argparse
from collections import Counter
from pathlib import Path
import numpy as np
from plastic_spatial_convergence import ROOT, inputs, lift, read, sha, write

GEOMETRY = ROOT/'tests/data/griffiths_corner_geometry.json'


def slope_mesh(geometry, divisions, depth_divisions=1):
    if divisions < 1 or divisions > 16 or divisions & (divisions-1) or depth_divisions < 1:
        raise ValueError('power-of-two in-plane subdivisions 1..16 and positive depth divisions required')
    xz=np.asarray(geometry['points_xz_m'],float)
    quads=geometry['quads']; fixed=np.asarray(geometry['fixed_x_z'],bool)
    if geometry['schema'] != 'soil-fem-slope-corner-geometry/1':
        raise ValueError('corner geometry schema required')
    owners=Counter(tuple(sorted((q[a],q[(a+1)%4]))) for q in quads for a in range(4))
    if set(owners.values()) != {1,2}:
        raise ValueError('invalid quad boundary topology')
    segments=[[xz[[a,b]] for (a,b),count in owners.items() if count==1 and fixed[a,axis] and fixed[b,axis]]
              for axis in range(2)]
    points=[]; lookup={}; triangles=[]; parents=[]
    def point(value):
        key=tuple(np.round(value,12))
        if key not in lookup:
            lookup[key]=len(points); points.append(np.asarray(value))
        return lookup[key]
    for parent,quad in enumerate(quads):
        corners=xz[quad]
        def mapped(i,j):
            r,s=i/divisions,j/divisions
            return point((1-r)*(1-s)*corners[0]+r*(1-s)*corners[1]+r*s*corners[2]+(1-r)*s*corners[3])
        for i in range(divisions):
            for j in range(divisions):
                a,b,c,d=[mapped(*ij) for ij in ((i,j),(i+1,j),(i+1,j+1),(i,j+1))]
                triangles.extend([[a,b,c],[a,c,d]]); parents.extend([parent,parent])
    count=len(points)
    nodes=[[float(x),layer/depth_divisions,float(z)] for layer in range(depth_divisions+1) for x,z in points]
    cells=[]; parent_quads=[]
    for layer in range(depth_divisions):
        for tri,parent in zip(triangles,parents,strict=True):
            a,b,c=sorted(tri)
            a0,b0,c0=[v+layer*count for v in (a,b,c)]
            a1,b1,c1=[v+(layer+1)*count for v in (a,b,c)]
            for tet in [[a0,b0,c0,c1],[a0,b0,b1,c1],[a0,a1,b1,c1]]:
                xyz=np.asarray([nodes[i] for i in tet])
                if np.linalg.det((xyz[1:]-xyz[0]).T)<0:
                    tet[1],tet[2]=tet[2],tet[1]
                cells.append(tet); parent_quads.append(parent)
    mesh=dict(schema='soil-fem-mesh/2',units='m-kN-kPa',axis='z_up',crs='LOCAL',origin_m=[0.,0.,0.],
              nodes=nodes,tetrahedra=cells,material_ids=[1]*len(cells),region_ids=[1]*len(cells),
              regions=[dict(id=1,body='Slope',material_id=1,bounds_m=[[0.,0.,0.],[32.,1.,10.]],
                  volume_m3=geometry['area_m2'],source_sha256=sha(GEOMETRY))],
              node_sets={},face_sets={},excavation_sets={})
    mesh,_=lift(mesh,[[0.,0.,0.]]*len(nodes))
    xyz=np.asarray(mesh['nodes']); sample=xyz[:,[0,2]]
    sets={'all':list(range(len(xyz)))}
    for axis,name in enumerate(('fixed_x','fixed_z')):
        mask=np.zeros(len(sample),bool)
        for a,b in segments[axis]:
            v=b-a; t=(sample-a)@v/(v@v)
            mask|=(t>=-1e-10)&(t<=1+1e-10)&(np.linalg.norm(sample-a-t[:,None]*v,axis=1)<1e-10)
        sets[name]=np.flatnonzero(mask).tolist()
    mesh['node_sets']=sets
    return mesh,parent_quads


def prepare(output, depth_divisions=1, reference_path=None):
    reference = None
    if reference_path is not None:
        from griffiths_2d import load_reference
        reference = load_reference(reference_path)
    geometry=read(GEOMETRY); output=Path(output); output.mkdir(parents=True,exist_ok=False)
    if reference is not None:
        from griffiths_2d import digest, GEOMETRY_DIGEST
        if digest(geometry) != GEOMETRY_DIGEST:
            raise ValueError('pinned slope corner geometry changed')
    _,ground,case=inputs(2)
    zero=ground['initialization']['geostatic']['seed']['states']['Soil']
    ground['initialization']['geostatic']['seed']['states']={'Slope':zero}
    ground['bodies']=[dict(body='Slope',material_id=1)]
    material=ground['materials'][0]['mohr_coulomb']
    material.update(young_kpa=100000.,poisson=.3,cohesion_kpa=10.,friction_deg=20.)
    ground['materials'][0]['name']='Griffiths unreduced F=1; no critical 2D displacement comparison'
    if reference is not None:
        material.update({k:v for k,v in reference['material'].items() if k!='unit_weight_kn_m3'})
        ground['materials'][0]['unit_weight_kn_m3']=reference['material']['unit_weight_kn_m3']
        ground['materials'][0]['name']='Pinned Griffiths 2D reference at archived fixed reduction factor'
        zero['void_ratio']=.8
    case['boundary']=dict(constraints=[dict(node_set=name,direction=np.eye(3)[axis].tolist(),value_m=0.)
        for name,axis in [('all',1),('fixed_x',0),('fixed_z',2)]])
    case['settings']['initial_increment']=1/32
    case['stages']=[dict(kind='boundary',name='retain_full_gravity_state',target_factor=0.,initial_increment=1.)]
    records=[]
    for n in (1,2,4):
        folder=output/f'xy{n}-ny{depth_divisions}'; folder.mkdir()
        mesh,parents=slope_mesh(geometry,n,depth_divisions)
        for name,value in [('mesh',mesh),('ground',ground),('case',case),('parent-quads',parents)]:
            write(folder/f'{name}.json',value)
        names=['mesh.json','ground.json','case.json','parent-quads.json']
        if reference is not None:
            from griffiths_2d import sample_nodes
            write(folder/'reference-map.json',sample_nodes(mesh,reference))
            # Existing full mechanics verifier reads this established filename.
            write(folder/'stages.json',case)
            names+=['reference-map.json','stages.json']
        records.append(dict(folder=folder.name,elements=len(mesh['tetrahedra']),nodes=len(mesh['nodes']),
                            relative_h_xz=1/n,depth_divisions=depth_divisions,
                            files_sha256={name:sha(folder/name) for name in names}))
    report=dict(schema='soil-fem-griffiths-in-plane-preparation/1',prepared=True,solved=False,
        p6_complete=False,griffiths_slope_five_percent_closed=False,material_reduction_factor=1.,
        source_geometry_sha256=sha(GEOMETRY),levels=records,
        note='No reference supplied: F=1 only. All ny are identical; this is not the archived F~1.35 material.')
    if reference is not None:
        from griffiths_2d import REFERENCE_DIGEST
        report.update(material_reduction_factor=reference['material_reduction_factor'],
            reference_2d_content_sha256=REFERENCE_DIGEST,
            note='Pinned original 661-node data supplied. Matching fixed reduced material; no solve or 5% verdict yet.')
    write(output/'preparation.json',report)
    return report


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--depth-divisions',type=int,default=1)
    parser.add_argument('--reference-2d',type=Path)
    args=parser.parse_args()
    prepare(args.output,args.depth_divisions,args.reference_2d)
