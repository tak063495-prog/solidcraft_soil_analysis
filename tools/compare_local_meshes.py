"""Record spatial discretization sensitivity separately from integration accuracy."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def data(run):
    mesh=json.loads((run/'ground/mesh.json').read_text());result=json.loads((run/'reference-refined-result.json').read_text())
    state=result['last_committed'];nodes=np.array(mesh['nodes']);tets=np.array(mesh['tetrahedra']);active=np.array(state['active_elements'])
    return mesh,nodes,tets,active,np.array(state['displacement_m'])


def interpolate(xyz,tets,u,queries):
    p=xyz[tets];inv=np.linalg.inv((p[:,1:]-p[:,:1]).transpose(0,2,1));lo=p.min(axis=1);hi=p.max(axis=1);values=[]
    for q in queries:
        candidates=np.flatnonzero(np.all((q>=lo-1e-8)&(q<=hi+1e-8),axis=1))
        local=np.einsum('eij,ej->ei',inv[candidates],q-p[candidates,0]);w=np.column_stack((1-local.sum(axis=1),local))
        inside=np.flatnonzero(np.all(w>=-1e-8,axis=1))
        if not len(inside):raise ValueError(f'query outside remaining soil: {q}')
        i=inside[0];values.append(w[i]@u[tets[candidates[i]]])
    return np.array(values)


def floor_mean(mesh,xyz,tets,active,u,name='Pit'):
    selected=set(mesh['excavation_sets'][name]);faces={};integral=0.;area=0.
    for ei,tet in enumerate(tets):
        for f in ((0,1,2),(0,1,3),(0,2,3),(1,2,3)):
            key=tuple(sorted(int(tet[i]) for i in f));faces.setdefault(key,[]).append(ei)
    for face,owners in faces.items():
        if len(owners)!=2 or sum(i in selected for i in owners)!=1:continue
        retained=next(i for i in owners if i not in selected)
        if not active[retained]:continue
        p=xyz[list(face)];n=np.cross(p[1]-p[0],p[2]-p[0]);center=p.mean(axis=0)
        if n@(center-xyz[tets[retained]].mean(axis=0))<0:n=-n
        if n[2]<=.5*np.linalg.norm(n):continue
        a=np.linalg.norm(n)/2;area+=a;integral+=a*u[list(face),2].mean()
    assert abs(area-16)<1e-8
    return integral/area


def compare(coarse_run,fine_run,output):
    if output.exists():raise ValueError('comparison output exists')
    cm,cxyz,ct,ca,cu=data(coarse_run);fm,fxyz,ft,fa,fu=data(fine_run)
    used=np.unique(ft[fa]);mapped=interpolate(cxyz,ct[ca],cu,fxyz[used])
    error=np.linalg.norm(mapped-fu[used])/np.linalg.norm(fu[used])
    controls=np.array([[5.,5.,6.],[3.,5.,8.],[1.,5.,8.]])
    ccontrol=interpolate(cxyz,ct[ca],cu,controls);fcontrol=interpolate(fxyz,ft[fa],fu,controls)
    cmean=floor_mean(cm,cxyz,ct,ca,cu);fmean=floor_mean(fm,fxyz,ft,fa,fu)
    assert cmean>0 and fmean>0
    report={'schema':'soil-fem-local-spatial-sensitivity/1','coarse_elements':len(ct),'fine_elements':len(ft),
            'coarse_floor_area_mean_u_z_m':float(cmean),'fine_floor_area_mean_u_z_m':float(fmean),
            'relative_change_floor_area_mean':float(abs(cmean/fmean-1)),
            'relative_change_interpolated_3d_displacement':float(error),
            'control_points_m':controls.tolist(),'coarse_control_u_z_m':ccontrol[:,2].tolist(),'fine_control_u_z_m':fcontrol[:,2].tolist(),
            'interpretation':'Two mesh levels quantify remaining spatial sensitivity; they do not establish mesh convergence or field accuracy.',
            'inputs_sha256':{str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for run in (coarse_run,fine_run) for p in (run/'ground/mesh.json',run/'reference-refined-result.json')}}
    output.write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('coarse',type=Path);p.add_argument('fine',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();compare(a.coarse,a.fine,a.output)
