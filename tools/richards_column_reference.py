"""Independent 1D cell-centred finite volumes and full Newton, input only.

This bounded reference does not use RT0 matrices, native heads, nonlinear
iterations or fluxes. Harmonic interior conductivity, prescribed H at the
bottom and a no-ponding atmospheric top. Spatial discretization differs from
the 3D mixed FEM and must be refined before interpreting differences.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[name]='1'
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def retention(psi,law):
    x=law['alpha_per_m']*np.maximum(-psi,0);n=law['n'];m=1-1/n
    se=(1+x**n)**(-m)
    theta=law['residual_water_content']+(law['saturated_water_content']-law['residual_water_content'])*se
    kr=np.sqrt(se)*(1-(1-se**(1/m))**m)**2
    offset=theta-law['saturated_water_content']+law['saturated_specific_storage_per_m']*np.maximum(psi,0)
    return theta,kr,offset


def thomas(lower,diagonal,upper,rhs):
    diagonal=diagonal.copy();rhs=rhs.copy()
    for i in range(1,len(rhs)):
        assert abs(diagonal[i-1])>1e-30,'singular reference Jacobian'
        ratio=lower[i-1]/diagonal[i-1];diagonal[i]-=ratio*upper[i-1];rhs[i]-=ratio*rhs[i-1]
    solution=rhs.copy();solution[-1]/=diagonal[-1]
    for i in range(len(rhs)-2,-1,-1):solution[i]=(rhs[i]-upper[i]*solution[i+1])/diagonal[i]
    return solution


def solve(case,cells,time_step):
    assert 2<=cells<=4000 and time_step>0
    mats=case['materials'];law=mats[0]['retention'];k=float(mats[0]['saturated_conductivity_m_s'][0][0])
    for mat in mats:
        assert mat['retention']==law and np.array_equal(mat['saturated_conductivity_m_s'],np.eye(3)*k)
        assert mat.get('volume_source_per_s',0)==0
    assert case['initial_head']=={'head_m':0.} and case['analysis']['kind']=='transient'
    dz=1/cells;z=(np.arange(cells)+.5)*dz;head=np.zeros(cells);time=0.;frames=[];steps=0
    total_precipitation=total_routed=total_boundary=0.;maximum_residual=0.
    initial_storage=float(retention(head-z,law)[2].sum()*dz)
    for period in case['analysis']['periods']:
        boundaries={b['face_set']:b['condition'] for b in period['boundaries']}
        bottom=boundaries['bottom']['value']['head_m'];rain=boundaries['top']['rainfall_vertical_m_s']
        assert boundaries['bottom']['kind']=='head' and boundaries['top']['kind']=='atmospheric_no_ponding'
        end=period['end_time_s']
        def field(h):
            theta,kr,storage=retention(h-z,law);conductivity=k*kr
            flux=np.zeros(cells+1)
            _,kb,_=retention(np.array([bottom]),law)
            flux[0]=-2/(1/(k*kb[0])+1/conductivity[0])*(h[0]-bottom)/(dz/2)
            flux[1:-1]=-2/(1/conductivity[:-1]+1/conductivity[1:])*np.diff(h)/dz
            capacity_flux=-2/(1/k+1/conductivity[-1])*(1-h[-1])/(dz/2)
            flux[-1]=max(-rain,capacity_flux)
            return theta,storage,flux
        while time<end:
            dt=min(time_step,end-time);old_storage=field(head)[1]
            def residual(h):
                _,storage,flux=field(h)
                return (storage-old_storage)*dz/dt+np.diff(flux)
            h=head.copy();iterations=0
            for iterations in range(100):
                value=residual(h);norm=float(np.max(np.abs(value)))
                if norm<1e-15:break
                lower=np.empty(cells-1);diagonal=np.empty(cells);upper=np.empty(cells-1)
                for color in range(3):
                    indices=np.arange(color,cells,3);step=1e-6*(1+np.abs(h[indices]));p=h.copy();m=h.copy();p[indices]+=step;m[indices]-=step
                    difference=residual(p)-residual(m)
                    for index,delta in zip(indices,step):
                        diagonal[index]=difference[index]/(2*delta)
                        if index>0:upper[index-1]=difference[index-1]/(2*delta)
                        if index+1<cells:lower[index]=difference[index+1]/(2*delta)
                correction=thomas(lower,diagonal,upper,-value)
                for half in range(20):
                    candidate=h+correction*2.**(-half)
                    if np.isfinite(candidate).all() and np.max(np.abs(residual(candidate)))<norm:
                        h=candidate;break
                else:raise AssertionError(f'FV Newton line search failed at t={time} dt={dt} residual={norm}')
            else:raise AssertionError(f'FV Newton limit t={time}')
            theta,storage,flux=field(h);maximum_residual=max(maximum_residual,float(np.max(np.abs(residual(h)))))
            total_precipitation+=rain*dt;total_routed+=(rain+flux[-1])*dt;total_boundary+=(flux[-1]-flux[0])*dt
            head=h;time+=dt;steps+=1
        frames.append(dict(name=period['name'],time_s=time,z_m=z.tolist(),head_m=head.tolist(),water_content=theta.tolist(),
            bottom_flux_upward_m_s=float(flux[0]),top_outward_flux_m_s=float(flux[-1]),storage_offset_m=float(storage.sum()*dz),
            total_precipitation_m=total_precipitation,total_surface_routed_water_m=total_routed,total_boundary_outflow_m=total_boundary))
    mass_error=abs(frames[-1]['storage_offset_m']-initial_storage+total_boundary)
    assert mass_error<1e-10,mass_error
    return dict(schema='soil-fem-independent-richards-column-fv/1',cells=cells,time_step_s=time_step,steps=steps,
        maximum_cell_mass_residual_m_s=maximum_residual,integrated_mass_error_m=mass_error,frames=frames,
        scope='Independent input-only 1D finite volumes and numerical full Newton. No native RT0 fields are used. This is a distinct discretization, not an exact 3D plastic/mechanical reference.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case',type=Path);p.add_argument('output',type=Path);p.add_argument('--cells',type=int,default=200);p.add_argument('--time-step',type=float,default=25.)
    a=p.parse_args();case=json.loads(a.case.read_text(encoding='utf-8-sig'));result=solve(case,a.cells,a.time_step)
    result['inputs_sha256']={str(path.resolve()):hashlib.sha256(path.read_bytes()).hexdigest() for path in (a.case,Path(__file__))}
    with a.output.open('x',encoding='utf-8') as stream:json.dump(result,stream,indent=2)
    print('Independent FV:',a.cells,'cells,',result['steps'],'steps, integrated mass error',result['integrated_mass_error_m'],'m')
