"""Compare physical locations/area averages; node averages change with a mesh."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import fem_reference as ref


def displacement(mesh,state,query):
    xyz=np.asarray(mesh['nodes']);u=np.asarray(state['displacement_m']);r0=np.asarray(query)
    for cell in ref.cells(mesh)[np.asarray(state['active_elements'])]:
        p=xyz[cell]
        if np.any(r0<p.min(axis=0)-1e-8) or np.any(r0>p.max(axis=0)+1e-8):continue
        r=np.linalg.solve((p[1:4]-p[0]).T,r0-p[0])
        for _ in range(20):
            shape,gradient=ref.basis(r[None],len(cell)==10);position=p[0]+shape[0]@(p-p[0]);residual=r0-position
            if np.linalg.norm(residual)<1e-10:break
            r+=np.linalg.solve((p-p[0]).T@gradient[0],residual)
        if min(1-r.sum(),*r)>=-1e-8 and np.linalg.norm(residual)<1e-9:return shape[0]@u[cell]
    raise ValueError(f'query outside remaining curved mesh: {query}')


def data(run):
    mesh=json.loads((run/'ground/mesh.json').read_text());result=json.loads((run/'stage-results.json').read_text());state=result['last_committed']
    assert result['complete'];u=np.asarray(state['displacement_m']);active=np.asarray(state['active_elements']);integral=0.;area=0.
    for key,owners in ref.topology(mesh).items():
        live=[i for i in owners if active[i]]
        if len(owners)!=2 or len(live)!=1:continue
        q=ref.face(mesh,live[0],key);vector=q['area_vector'].sum(axis=0)
        if vector[2]<=.5*np.linalg.norm(vector):continue
        weights=np.linalg.norm(q['area_vector'],axis=1);area+=weights.sum();integral+=weights@(q['shape']@u[q['nodes'],2])
    return {'elements':len(mesh['tetrahedra']),'floor_area_m2':float(area),'floor_area_mean_u_z_m':float(integral/area),
        'controls_u_m':[displacement(mesh,state,q).tolist() for q in [[5,5,6],[2,5,8],[0,5,8]]],
        'nodal_monitors':result['frames'][-1]['monitors']}


def compare(coarse,fine,output):
    if output.exists():raise ValueError('comparison exists')
    a=data(coarse);b=data(fine)
    report={'schema':'soil-fem-quadratic-spatial-sensitivity/1','control_locations_m':[[5,5,6],[2,5,8],[0,5,8]],
        'coarse':a,'fine':b,'floor_area_mean_relative_change':float(abs(a['floor_area_mean_u_z_m']/b['floor_area_mean_u_z_m']-1)),
        'interpretation':'Two meshes quantify remaining sensitivity; no spatial convergence or physical calibration claim. Nodal means are mesh-dependent.',
        'inputs_sha256':{str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for run in (coarse,fine) for p in (run/'ground/mesh.json',run/'stage-results.json')}}
    output.write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('coarse',type=Path);p.add_argument('fine',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();compare(a.coarse,a.fine,a.output)
