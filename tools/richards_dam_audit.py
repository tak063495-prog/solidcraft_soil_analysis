"""Batch reassembly for the straight, homogeneous isotropic dam fixture.

No independent solution is generated here. FV is a separate solver. Arrays
are cached/batched to avoid the small-fixture helper's repeated mesh conversion.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import fem_reference as fem
from richards_column_reference import retention


def audit(mesh,result):
    assert result['complete'] and not mesh.get('tetrahedra10') and result['case']['analysis']['kind']=='steady'
    assert all(result['active_elements']) and result['cell_element_indices']==list(range(len(mesh['tetrahedra'])))
    state=result['last_committed'];cells=np.asarray(mesh['tetrahedra']);xyz=np.asarray(mesh['nodes']);count=len(cells)
    heads=np.asarray(state['relative_cell_head_m'])+result['case']['head_datum_m'];flows=np.asarray(state['cell_outflow_m3_s']);lambdas=np.zeros((count,4))
    law=result['case']['materials'][0]['retention'];ks=result['case']['materials'][0]['saturated_conductivity_m_s'][0][0]
    for m in result['case']['materials']:assert m['retention']==law and np.array_equal(m['saturated_conductivity_m_s'],np.eye(3)*ks) and not m.get('volume_source_per_s',0)
    for fi,face in enumerate(result['faces']):
        for e,j in face['owners']:lambdas[e,j]=state['relative_face_head_m'][fi]+result['case']['head_datum_m']
    natural,weights=fem.volume_rule({'tetrahedra10':[True],'quadrature_order':3});rv=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    ref=2*(natural[:,None,:]-rv[None,:,:]);diagonals=np.zeros((count,4));maximum_constitutive=0.
    for first in range(0,count,2048):
        last=min(count,first+2048);p=xyz[cells[first:last]];j=np.swapaxes(p[:,1:,:]-p[:,0,None,:],1,2);det=np.linalg.det(j);assert np.all(det>0)
        position=p[:,0,None,:]+np.einsum('eij,qj->eqi',j,natural)
        basis=np.einsum('eij,qfj->eqfi',j,ref)/det[:,None,None,None]
        kr=retention(heads[first:last,None]-position[:,:,2],law)[1]
        mass=np.einsum('eqfi,eqgi,eq->efg',basis,basis,weights[None,:]*det[:,None]/(ks*kr))
        residual=np.einsum('eij,ej->ei',mass,flows[first:last])+lambdas[first:last]-heads[first:last,None]
        maximum_constitutive=max(maximum_constitutive,float(np.max(np.abs(residual))))
        diagonals[first:last]=np.diagonal(np.linalg.inv(mass),axis1=1,axis2=2)
    assert maximum_constitutive<1e-9,maximum_constitutive
    settings=result['case']['settings'];tolerance=settings['absolute_flow_tolerance_m3_s']+settings['relative_flow_tolerance']*np.abs(flows).sum(1)
    mass_ratio=float(np.max(np.abs(flows.sum(1))/tolerance));assert mass_ratio<=1.02,mass_ratio
    max_jump=0.;head_conditions={tuple(sorted(f)) for f in mesh['face_sets']['reservoir']};air_conditions={tuple(sorted(f)) for f in mesh['face_sets']['air']};q_in=q_out=right=seep_area=0.
    for fi,face in enumerate(result['faces']):
        q=np.asarray([flows[e,j] for e,j in face['owners']]);flux=q.sum();tol=settings['absolute_flow_tolerance_m3_s']+settings['relative_flow_tolerance']*np.abs(q).sum()
        if len(q)==2:
            max_jump=max(max_jump,abs(flux));assert abs(flux)<=1.02*tol;continue
        corners=xyz[np.asarray(face['corners'])];mean=corners.mean(0);area=np.linalg.norm(np.cross(corners[1]-corners[0],corners[2]-corners[0]))/2
        assert abs(area-face['area_m2'])<1e-12
        h=state['relative_face_head_m'][fi]+result['case']['head_datum_m']
        if tuple(face['corners']) in head_conditions:assert abs(h-.5)<1e-8;q_in-=flux
        elif tuple(face['corners']) in air_conditions:
            psi=h-mean[2];e,j=face['owners'][0];t=diagonals[e,j]
            assert psi<=1.01*settings['head_tolerance_m'] and flux>=-1.02*tol
            assert abs(min(-psi,flux/t))<=1.01*settings['head_tolerance_m']+1e-14
            q_out+=flux
            if abs(mean[0]-1)<1e-12:
                right+=flux
                if abs(psi)<=settings['head_tolerance_m'] and flux>tol:seep_area+=area
        else:assert abs(flux)<=1.02*tol
    assert abs(q_out-q_in)<1e-12
    return dict(schema='soil-fem-dam-discrete-audit/1',passed=True,cells=count,
        maximum_constitutive_head_residual_m=maximum_constitutive,maximum_scaled_cell_mass_residual=mass_ratio,maximum_internal_flux_jump_m3_s=max_jump,
        right_outflow_m3_s=right,seepage_area_m2=seep_area,independent_boundary_mass_error_m3_s=abs(q_out-q_in),
        scope='All straight RT0 cells reassembled in batches from primary DOFs; discrete audit only. Independent nonlinear solution is supplied separately by FV.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mesh',type=Path);p.add_argument('result',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    result=audit(json.loads(a.mesh.read_text(encoding='utf-8')),json.loads(a.result.read_text(encoding='utf-8')))
    result['inputs_sha256']={str(f.resolve()):hashlib.sha256(f.read_bytes()).hexdigest() for f in (a.mesh,a.result,Path(__file__))}
    with a.output.open('x',encoding='utf-8') as stream:json.dump(result,stream,indent=2)
    print(result)
