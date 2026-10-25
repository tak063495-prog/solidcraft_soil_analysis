"""Reassemble local constitutive equations and water balances from primary DOFs.

This checks the discrete solution on independent NumPy geometry. It does not
claim an independent nonlinear 3D solution; use the 1D FV/convergence suite.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import fem_reference as fem
from surface_projection_reference import upward_area, surface_area
from richards_io import geometry,fields


def audit(mesh,result):
    assert result['complete'] and result['schema']=='soil-fem-richards-result/1'
    active=np.flatnonzero(result['active_elements']);assert active.tolist()==result['cell_element_indices']
    geos={e:geometry(mesh,e) for e in active};volumes={e:geos[e]['weight'].sum() for e in active}
    materials={m['id']:m for m in result['case']['materials']};settings=result['case']['settings'];datum=result['case']['head_datum_m']
    topology={key:[e for e in owners if result['active_elements'][e]] for key,owners in fem.topology(mesh).items()}
    topology={k:v for k,v in topology.items() if v};assert set(topology)=={tuple(f['corners']) for f in result['faces']}
    face_geometry=[];internal_area_quadrature_difference=internal_projection_quadrature_difference=0.
    for f in result['faces']:
        owners=[int(active[ci]) for ci,_ in f['owners']];assert sorted(owners)==sorted(topology[tuple(f['corners'])])
        actual=fem.face(mesh,owners[0],f['corners']);coords=np.asarray(mesh['nodes'])[actual['nodes']]
        mean=coords[:3].mean(0) if len(coords)==3 else coords[3:].mean(0)
        area,area_reference_error=surface_area(coords);projected,projection_reference_error=upward_area(coords)
        assert np.allclose(mean,f['reference_position_m'],rtol=1e-13,atol=1e-13)
        assert abs(mean[2]-f['reference_elevation_m'])<1e-12
        if len(owners)==1:
            assert abs(area-f['area_m2'])<3e-12*area+area_reference_error+1e-15
            assert abs(projected-f['projected_horizontal_area_m2'])<f.get('projected_horizontal_area_error_bound_m2',0.)+projection_reference_error+1e-12*area+1e-15
        else:
            # Native order-4 vs independent order-5 surface quadrature. Internal
            # RT0 normal flux uses the Piola face integral, not these areas.
            internal_area_quadrature_difference=max(internal_area_quadrature_difference,abs(area-f['area_m2']))
            internal_projection_quadrature_difference=max(internal_projection_quadrature_difference,abs(projected-f['projected_horizontal_area_m2']))
        face_geometry.append((mean,area,projected,projection_reference_error))
    periods={p['name']:p['boundaries'] for p in result['case']['analysis'].get('periods',[])}
    accepted=[a for a in result['attempts'] if a['accepted']];assert len(accepted)==len(result['frames'])-1
    rows=[];old=None;initial_storage=None;integrated_boundary=integrated_source=integrated_rain=integrated_routed=0.
    maximum_flux_error=maximum_storage_error=0.
    for index,frame in enumerate(result['frames']):
        state=frame['state'];blocks,local=fields(mesh,result,frame,geos)
        maximum_flux_error=max(maximum_flux_error,local['maximum_independent_local_flux_error_m3_s']);maximum_storage_error=max(maximum_storage_error,local['maximum_independent_storage_error_m3'])
        assert local['maximum_independent_local_flux_error_m3_s']<1e-15
        assert local['maximum_independent_storage_error_m3']<1e-12*sum(volumes.values())
        storage={b['element']:b['storage'] for b in blocks};total_storage=sum(storage.values())
        if index==0:initial_storage=total_storage;old=storage;continue
        attempt=accepted[index-1];dt=attempt['time_step_s'];assert attempt['period']==frame['name']
        assert state['time_s']==attempt['start_time_s']+(dt or 0.)
        definitions=periods[frame['name']] if periods else result['case']['analysis']['boundaries'];conditions={}
        for b in definitions:
            for key in mesh['face_sets'][b['face_set']]:conditions[tuple(sorted(key))]=b['condition']
        maximum_cell_error=maximum_face_error=boundary=rain=routed=source=0.
        for e in active:
            q=np.asarray(state['cell_outflow_m3_s'][e]);s=volumes[e]*materials[result['material_ids'][e]].get('volume_source_per_s',0.);source+=s
            rate=(storage[e]-old[e])/dt if dt else 0.
            residual=abs(q.sum()+rate-s);tol=settings['absolute_flow_tolerance_m3_s']+settings['relative_flow_tolerance']*(np.abs(q).sum()+abs(rate)+abs(s))
            # Independent summation/quadrature roundoff bound in storage/dt.
            roundoff=2e-15*volumes[e]/dt if dt else 0.
            assert residual<=1.02*tol+roundoff,(index,e,residual,tol)
            maximum_cell_error=max(maximum_cell_error,residual)
        for fi,f in enumerate(result['faces']):
            flows=[state['cell_outflow_m3_s'][active[ci]][j] for ci,j in f['owners']];q=sum(flows)
            mean,area,projected,projection_reference_error=face_geometry[fi];condition=conditions.get(tuple(f['corners']),{'kind':'outward_flux','flux_m_s':0.})
            tol=settings['absolute_flow_tolerance_m3_s']+settings['relative_flow_tolerance']*sum(map(abs,flows))
            if len(flows)==2:assert abs(q)<=1.02*tol,(index,fi,q,tol);maximum_face_error=max(maximum_face_error,abs(q));continue
            boundary+=q;head=state['relative_face_head_m'][fi];psi=datum-mean[2]+head
            if condition['kind']=='head':
                value=condition['value'];target=value['head_m']-datum+np.dot(value.get('gradient',[0.,0.,0.]),mean-value.get('reference_point_m',[0.,0.,0.]))
                assert abs(head-target)<=settings['head_tolerance_m']
            elif condition['kind']=='pressure':assert abs(psi-condition['pressure_head_m'])<=settings['head_tolerance_m']
            elif condition['kind']=='outward_flux':
                target=condition['flux_m_s']*area;assert abs(q-target)<=1.02*(tol+settings['relative_flow_tolerance']*abs(target))+1e-18*area
            else:
                rainfall=condition['rainfall_vertical_m_s']*projected;excess=q+rainfall
                rain_geometry_error=condition['rainfall_vertical_m_s']*(f.get('projected_horizontal_area_error_bound_m2',0.)+projection_reference_error+1e-12*area)
                assert psi<=1.01*settings['head_tolerance_m'] and excess>=-1.02*(tol+settings['relative_flow_tolerance']*rainfall)-rain_geometry_error
                ci,j=f['owners'][0];t=blocks[ci]['face_conductance_m2_s'][j]
                assert abs(min(-psi,excess/t))<=1.01*settings['head_tolerance_m']+1e-14+rain_geometry_error/t,(index,fi,psi,excess/t)
                rain+=rainfall;routed+=excess
        if dt:
            integrated_boundary+=boundary*dt;integrated_source+=source*dt;integrated_rain+=rain*dt;integrated_routed+=routed*dt
        rows.append(dict(name=frame['name'],time_s=state['time_s'],maximum_independent_cell_mass_residual_m3_s=maximum_cell_error,
            maximum_independent_face_flux_jump_m3_s=maximum_face_error,storage_offset_m3=total_storage,precipitation_m3_s=rain,surface_routed_water_m3_s=routed))
        old=storage
    integrated_error=abs(total_storage-initial_storage+integrated_boundary-integrated_source) if periods else None
    if periods:assert integrated_error<1e-9*sum(volumes.values())+1e-14,integrated_error
    return dict(schema='soil-fem-richards-independent-reassembly/1',passed=True,frames=rows,
        maximum_local_constitutive_flux_error_m3_s=maximum_flux_error,maximum_storage_reassembly_error_m3=maximum_storage_error,
        integrated_mass_error_m3=integrated_error,total_precipitation_m3=integrated_rain,total_surface_routed_water_m3=integrated_routed,
        maximum_internal_face_area_order4_to_5_difference_m2=internal_area_quadrature_difference,
        maximum_internal_projected_area_order4_to_5_difference_m2=internal_projection_quadrature_difference,
        scope='Independent NumPy local RT0/geometry/retention and balance reassembly from native primary DOFs. No independent nonlinear 3D solution is inferred.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mesh',type=Path);p.add_argument('result',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    report=audit(json.loads(a.mesh.read_text(encoding='utf-8-sig')),json.loads(a.result.read_text(encoding='utf-8-sig')))
    report['sha256']={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in (a.mesh,a.result,Path(__file__))}
    with a.output.open('x',encoding='utf-8') as stream:json.dump(report,stream,indent=2)
    print('Richards discrete audit passed:',len(report['frames']),'accepted steps; mass error',report['integrated_mass_error_m3'])
