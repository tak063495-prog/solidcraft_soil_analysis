"""Independent NumPy geometry/constitutive reassembly of stored mixed DOFs.

An audit/reconstruction, not an independent nonlinear PDE solution. The
separate finite-volume reference supplies that comparison for 1D columns.
"""
import numpy as np
import fem_reference as fem
from richards_column_reference import retention


def geometry(mesh,e):
    cell=fem.cells(mesh)[e];p=np.asarray(mesh['nodes'])[cell]
    # RT0 mass and nonlinear storage need a positive volume rule also for Tet4.
    r,w=fem.volume_rule({'tetrahedra10':[True],'quadrature_order':mesh.get('quadrature_order',3)})
    shape,gradient=fem.basis(r,len(cell)==10);j=np.einsum('ai,qaj->qij',p-p[0],gradient);det=np.linalg.det(j)
    assert np.all(det>0)
    rv=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    basis=np.einsum('qij,qfj->qfi',j,2*(r[:,None,:]-rv[None,:,:]))/det[:,None,None]
    return dict(position=p[0]+shape@(p-p[0]),weight=w*det,basis=basis)


def fields(mesh,result,frame,geometries=None):
    active=np.flatnonzero(result['active_elements']);state=frame['state'];datum=result['case']['head_datum_m']
    assert active.tolist()==result['cell_element_indices']
    mats={m['id']:m for m in result['case']['materials']};ids=result['material_ids']
    lambdas=np.zeros((len(active),4))
    for fi,face in enumerate(result['faces']):
        for ci,local in face['owners']:lambdas[ci,local]=state['relative_face_head_m'][fi]
    blocks=[];max_flux_error=max_storage_error=0.
    for ci,e in enumerate(active):
        geo=geometries[e] if geometries is not None else geometry(mesh,e)
        h=state['relative_cell_head_m'][e];psi=(datum-geo['position'][:,2])+h;mat=mats[ids[e]]
        theta,kr,offset=retention(psi,mat['retention']);assert np.isfinite(kr).all() and np.all(kr>0)
        flux=np.asarray(state['cell_outflow_m3_s'][e]);q=np.einsum('qfi,f->qi',geo['basis'],flux)
        inverse=np.linalg.inv(mat['saturated_conductivity_m_s'])[None,:,:]/kr[:,None,None]
        mass=np.einsum('qfi,qij,qgj,q->fg',geo['basis'],inverse,geo['basis'],geo['weight'])
        expected=np.linalg.solve(mass,h-lambdas[ci]);conductance=np.diag(np.linalg.inv(mass))
        error=float(np.max(np.abs(expected-flux)));scale=float(np.max(np.abs(flux)))
        assert error<1e-25+1e-9*scale,(int(e),error,scale)
        max_flux_error=max(max_flux_error,float(np.max(np.abs(expected-flux))))
        storage=float(geo['weight']@offset)
        max_storage_error=max(max_storage_error,abs(storage-state['cell_storage_offset_m3'][e]))
        law=mat['retention'];se=(theta-law['residual_water_content'])/(law['saturated_water_content']-law['residual_water_content'])
        blocks.append(dict(element=e,**geo,head=np.full(len(psi),datum+h),pressure_head=psi,water_content=theta,
            effective_saturation=se,water_saturation=theta/law['saturated_water_content'],relative_conductivity=kr,flux=q,storage=storage,face_conductance_m2_s=conductance))
    if frame['points']:
        expected={b['element']:b for b in blocks};grouped={e:[] for e in active}
        for q in frame['points']:grouped[q['element_index']].append(q)
        for e,points in grouped.items():
            b=expected[e];assert len(points)==len(b['weight'])
            for native,field in [('position_m','position'),('weight_m3','weight'),('pressure_head_m','pressure_head'),('total_head_m','head'),('darcy_flux_m_s','flux')]:
                value=np.asarray([p[native] for p in points]);assert np.allclose(value,b[field],rtol=1e-10,atol=1e-14 if field!='flux' else 1e-18),(e,native)
            for native,field in [('water_content','water_content'),('water_saturation','water_saturation'),('effective_saturation','effective_saturation'),('relative_conductivity','relative_conductivity')]:
                assert np.allclose([p['retention'][native] for p in points],b[field],rtol=1e-10,atol=1e-14),(e,native)
    return blocks,dict(maximum_independent_local_flux_error_m3_s=max_flux_error,maximum_independent_storage_error_m3=max_storage_error)
