"""Independent Cartesian FV Richards + zero-thickness contact diagnostics.

No project solver imports. See README for the explicit synthetic interface model.
Numerical completion does not confer A-04/C-04 acceptance.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

PARAMETERS = dict(ks_m_s=.02, theta_r=.08, theta_s=.4117647058823529,
                  alpha_per_m=1.2, n=1.6, ss_per_m=.02,
                  gamma_water_kn_m3=10., external_opening_stiffness_kpa_m=2000.,
                  tangential_contact_stiffness_kpa_m=200., mu=.3,
                  gap_ks_m_s=.002, one_d_elevation_m=.32, normal_load_offset_kpa=.1,
                  normal_load_amplitude_kpa=.5, boundary_head_amplitude_m=.15,
                  prescribed_head_amplitude_m=.07, shear_displacement_amplitude_m=.006,
                  schedule_peak_aperture_m=.0005, schedule_events_s=[.35,.75],
                  problem_definition='unit_cartesian_domain;left_Hbase+.15*sin(pi*t)^2*(1+.2*cos(pi*z));right_Hbase;N=.1+.5*cos(2*pi*t)^2;s=.006*sin(2*pi*t)')
EVENTS = (.35, .75)


def material(head, elevation, law='vg'):
    """w, dw/dH, K, dK/dH; same explicit VG/Mualem equations as docs/RICHARDS."""
    p = np.asarray(head)-elevation
    if law == 'linear':
        return .4117647058823529+.02*p, np.full_like(p,.02), np.full_like(p,.02), np.zeros_like(p)
    if law != 'vg':
        raise ValueError('unknown material')
    n=PARAMETERS['n']; m=1-1/n; alpha=PARAMETERS['alpha_per_m']
    a=alpha*np.maximum(-p,0); u=1+a**n; se=u**(-m)
    dse=np.where(p<0, m*n*alpha*a**(n-1)*u**(-m-1),0.)
    theta=.08+(.4117647058823529-.08)*se
    w=theta+.02*np.maximum(p,0)
    dw=(.4117647058823529-.08)*dse+.02*(p>=0)
    # Differentiate using a: B=1-a^(n-1)*(1+a^n)^(-m).
    b=1-a**(n-1)*u**(-m)
    db=np.zeros_like(a)
    mask=a>0
    db[mask]=alpha*((n-1)*a[mask]**(n-2)*u[mask]**(-m)
                          -m*n*a[mask]**(2*n-2)*u[mask]**(-m-1))
    k=.02*np.sqrt(se)*b*b
    dk=.02*(.5/np.sqrt(se)*dse*b*b+2*np.sqrt(se)*b*db)
    return w,dw,k,dk


def contact(head, z, t, plastic, dry=False):
    """Compression-positive normal force and 1D Coulomb return mapping."""
    pressure=np.zeros_like(head) if dry else 10*np.maximum(head-z,0)
    load=.1+.5*np.cos(2*np.pi*t)**2
    gap=np.maximum((pressure-load)/2000,0)
    effective=np.maximum(load-pressure,0)
    displacement=.006*np.sin(2*np.pi*t)
    trial=200*(displacement-plastic)
    shear=np.clip(trial,-.3*effective,.3*effective)
    # Open spring has zero traction; keep its archived plastic slip frozen.
    new_plastic=np.where(effective>0,displacement-shear/200,plastic)
    state=np.where(gap>0,'open',np.where(np.abs(trial)>.3*effective,'slip','stick'))
    return dict(pressure=pressure,gap=gap,effective=effective,shear=shear,
                plastic=new_plastic,state=state,load=np.full_like(head,load))


class Model:
    def __init__(self,nx=8,ny=1,law='vg',mode='ONE_WAY',policy='CONTINUOUS',dry=False,
                 switch_shift=0., schedule=True, eigen=False):
        if nx<4 or nx%2 or ny<1 or law not in ('vg','linear'):
            raise ValueError('even nx>=4, ny>=1 and known law required')
        if mode not in ('ONE_WAY','TWO_WAY') or policy not in ('CONTINUOUS','PRESCRIBED'):
            raise ValueError('unknown mode/policy; transient SOIL mapping is unconfirmed')
        if policy=='PRESCRIBED' and mode=='TWO_WAY':
            raise ValueError('contact-triggered reservoir activation is not defined')
        if law=='linear' and dry:
            pass
        self.nx=nx; self.ny=ny; self.dx=1/nx; self.dz=1/ny
        self.law=law;self.mode=mode;self.policy=policy;self.dry=dry
        self.shift=switch_shift;self.schedule=schedule;self.eigen=eigen
        self.z=np.repeat((np.arange(ny)+.5)/ny,nx) if ny>1 else np.full(nx,.32)
        self.zj=(np.arange(ny)+.5)/ny if ny>1 else np.array([.32])
        self.n=nx*ny;self.size=self.n+ny;self.volume=self.dx*self.dz
        self.base=2. if law=='linear' else .3
        self.edges=[]
        for j in range(ny):
            for i in range(nx):
                a=j*nx+i
                if i+1<nx and i!=nx//2-1:self.edges.append((a,a+1,self.dx,self.dz))
                if j+1<ny:self.edges.append((a,a+nx,self.dz,self.dx))

    def aperture(self,h,t):
        if self.mode=='ONE_WAY':
            start=round(EVENTS[0]+self.shift,12);end=round(EVENTS[1]+self.shift,12)
            on=self.schedule and start < t < end
            phase=(t-start)/(end-start)
            aperture=.0005*np.sin(np.pi*phase)**2 if on else 0.
            return np.full(self.ny,aperture),np.zeros(self.ny)
        c=contact(h,self.zj,t,np.zeros(self.ny),self.dry)
        derivative=np.where(c['gap']>0,10./2000*(h>self.zj),0.) if not self.dry else np.zeros(self.ny)
        return c['gap'],derivative

    def storage(self,x,t):
        w=material(x[:self.n],self.z,self.law)[0]*self.volume
        a=self.aperture(x[self.n:],t)[0]
        wj=material(x[self.n:],self.zj,self.law)[0]*a*self.dz
        return np.r_[w,wj]

    def residual(self,x,old,old_time,t,jacobian=True):
        dt=t-old_time;h=x[:self.n];hj=x[self.n:]
        w,dw,_,_=material(h,self.z,self.law)
        wj,dwj,kg,dkg=material(hj,self.zj,self.law)
        kg=.1*kg;dkg=.1*dkg
        a,da=self.aperture(hj,t)
        r=(self.storage(x,t)-self.storage(old,old_time))/dt
        J=lil_matrix((self.size,self.size))
        J.setdiag(np.r_[dw*self.volume/dt,(a*dwj+da*wj)*self.dz/dt])
        exterior=0.; joint_exchange=0.; qleft=[];qright=[];traces=[];joint_sources=[];joint_rates=[]
        def edge(ia,ib,distance,area):
            elevation=(self.z[ia]+self.z[ib])/2
            _,_,k,dk=material(np.array([(x[ia]+x[ib])/2]),np.array([elevation]),self.law)
            grad=(x[ia]-x[ib])*area/distance
            q=float(k[0]*grad)
            qa=float(k[0]*area/distance+.5*dk[0]*grad)
            qb=float(-k[0]*area/distance+.5*dk[0]*grad)
            r[ia]+=q;r[ib]-=q
            J[ia,ia]+=qa;J[ia,ib]+=qb;J[ib,ia]-=qa;J[ib,ib]-=qb
        for ia,ib,distance,area in self.edges:edge(ia,ib,distance,area)
        # Outer left/right specified heads, top/bottom no-flow.
        for j in range(self.ny):
            for i in (0,self.nx-1):
                c=j*self.nx+i
                boundary=self.base+(0. if self.eigen else .15*np.sin(np.pi*t)**2*(1+.2*np.cos(np.pi*self.z[c])) if i==0 else 0.)
                _,_,k,dk=material(np.array([.5*(h[c]+boundary)]),np.array([self.z[c]]),self.law)
                factor=self.dz/(self.dx/2)
                q=float(k[0]*(h[c]-boundary)*factor)
                r[c]+=q;J[c,c]+=factor*(k[0]+.5*dk[0]*(h[c]-boundary));exterior+=q
            idx=self.n+j
            flow=[];surface=[]
            for c in (j*self.nx+self.nx//2-1,j*self.nx+self.nx//2):
                # Soil half-cell + half-gap resistances in series. Kgap is an
                # explicit synthetic filled-gap Darcy coefficient, no cubic law.
                _,_,k,dk=material(np.array([h[c]]),np.array([self.zj[j]]),self.law)
                resistance=self.dx/(2*k[0])+a[j]/(2*kg[j])
                derivative_k=-self.dx*dk[0]/(2*k[0]**2)
                dh=h[c]-hj[j];q=self.dz*dh/resistance
                qc=self.dz/resistance-self.dz*dh/resistance**2*derivative_k
                qj=-self.dz/resistance-self.dz*dh/resistance**2*(da[j]/(2*kg[j])-a[j]*dkg[j]/(2*kg[j]**2))
                r[c]+=q; r[idx]-=q
                J[c,c]+=qc;J[c,idx]+=qj;J[idx,c]-=qc;J[idx,idx]-=qj
                flow.append(float(q));surface.append(float(h[c]-q/self.dz*self.dx/(2*k[0])))
            qleft.append(flow[0]);qright.append(-flow[1]);traces.append(surface)
            joint_rates.append(float((self.storage(x,t)[idx]-self.storage(old,old_time)[idx])/dt))
            source=0.
            opened=bool(a[j]>0)
            if self.policy=='PRESCRIBED' and opened:
                target=self.base+.07*np.sin(np.pi*t)**2
                # Interface reservoir supplies actual gap storage plus both soils.
                source=float(r[idx]);joint_exchange+=source
                r[idx]=hj[j]-target;J.rows[idx]=[idx];J.data[idx]=[1.]
            joint_sources.append(source)
        return r,J.tocsr(),dict(exterior_outflow=exterior,reservoir_inflow=joint_exchange,
                               q_left=qleft,q_right=qright,soil_traces=traces,aperture=a.tolist(),
                               joint_sources=joint_sources,joint_storage_rates=joint_rates)

    def solve_step(self,old,old_time,t):
        x=old.copy()
        for iteration in range(35):
            r,J,info=self.residual(x,old,old_time,t)
            residual=float(np.max(np.abs(r)))
            if residual<2e-12:return x,info,residual,iteration
            step=spsolve(J,-np.asarray(r,dtype=float))
            if not np.all(np.isfinite(step)):raise ArithmeticError('nonfinite Newton step')
            alpha=1.
            for _ in range(22):
                trial=x+alpha*step
                rr=self.residual(trial,old,old_time,t,False)[0]
                if np.max(np.abs(rr))<residual:break
                alpha*=.5
            else:raise ArithmeticError('line search failed')
            x=trial
        raise ArithmeticError('Newton iteration limit')

    def run(self,dt=.025,end=1.):
        if dt<=0 or end<=0:raise ValueError('positive physical dt/end required')
        times=sorted(set([0.,round(end,12)]+[round(min(i*dt,end),12) for i in range(1,int(np.ceil(end/dt))+1)]
                         +[round(e+self.shift,12) for e in EVENTS if 0<e+self.shift<end]))
        x=np.full(self.size,self.base,dtype=np.longdouble)
        if self.eigen:
            xx=np.tile((np.arange(self.nx)+.5)/self.nx,self.ny)
            x[:self.n]+=.1*np.sin(np.pi*xx)*(np.cos(np.pi*self.z) if self.ny>1 else 1.)
            x[self.n:]+=.1*(np.cos(np.pi*self.zj) if self.ny>1 else 1.)
        initial=self.storage(x,0).sum(); net=0.
        initial_head=x.copy();frames=[]
        plastic=np.zeros(self.ny);dissipation=np.zeros(self.ny);rows=[];maxres=0.;iterations=0
        for previous,t in zip(times,times[1:]):
            x,info,res,it=self.solve_step(x,previous,t);maxres=max(maxres,res);iterations+=it
            net+=(t-previous)*(info['reservoir_inflow']-info['exterior_outflow'])
            c=contact(x[self.n:],self.zj,t,plastic,self.dry)
            dissipation+=np.abs(c['shear']*(c['plastic']-plastic));plastic=c['plastic']
            defect=float(self.storage(x,t).sum()-initial-net)
            frames.append(dict(time_s=t,dt_s=t-previous,head=x.astype(float).tolist(),storage_m3=self.storage(x,t).astype(float).tolist(),
                               exterior_outflow_m3_s=info['exterior_outflow'],reservoir_inflow_m3_s=info['reservoir_inflow'],
                               cumulative_net_inflow_m3=float(net),water_balance_defect_m3=defect,
                               local_residual_maximum=res))
            for j in range(self.ny):
                rows.append(dict(time_s=t,joint_id='J0',point_id=j,generation=0,
                                 head_a_m=info['soil_traces'][j][0],head_b_m=info['soil_traces'][j][1],
                                 joint_head_m=float(x[self.n+j]),q_a_to_joint_m3_s=info['q_left'][j],
                                 q_joint_to_b_m3_s=info['q_right'][j],
                                 reservoir_inflow_m3_s=info['joint_sources'][j],joint_storage_rate_m3_s=info['joint_storage_rates'][j],
                                 hydraulic_open=bool(info['aperture'][j]>0),contact_state=str(c['state'][j]),
                                 gap_m=float(c['gap'][j]),hydraulic_aperture_m=float(info['aperture'][j]),
                                 physical_pressure_kpa=float(c['pressure'][j]),effective_pressure_kpa=float(c['effective'][j]),
                                 shear_kpa=float(c['shear'][j]),plastic_slip_m=float(plastic[j]),
                                 friction_dissipation_kpa_m=float(dissipation[j]),
                                 water_balance_defect_m3=defect))
        switches=[]
        for j in range(self.ny):
            prior=contact(np.array([self.base]),np.array([self.zj[j]]),0.,np.zeros(1),self.dry)['state'][0]=='open'
            last_time=0.
            for row in (r for r in rows if r['point_id']==j):
                opened=row['contact_state']=='open'
                if opened!=prior:switches.append(dict(point_id=j,bracket_s=[last_time,row['time_s']],new_open=opened))
                prior=opened;last_time=row['time_s']
        return dict(configuration=dict(nx=self.nx,ny=self.ny,law=self.law,mode=self.mode,policy=self.policy,
                                       dry=self.dry,dt_s=dt,end_s=end,switch_shift_s=self.shift,schedule=self.schedule,eigen=self.eigen),
                    complete=True,accepted=None,A04='未達',C04='未達',rows=rows,contact_switch_brackets=switches,
                    initial_head=initial_head.astype(float).tolist(),hydraulic_frames=frames,final_head=x.astype(float).tolist(),initial_storage_m3=float(initial),final_storage_m3=float(self.storage(x,times[-1]).sum()),
                    cumulative_net_inflow_m3=float(net),maximum_water_balance_defect_m3=max(abs(r['water_balance_defect_m3']) for r in rows),
                    maximum_residual_m3_s=maxres,newton_iterations=iterations)


def write(path,value):
    Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True)
    p.add_argument('--nx',type=int,default=8);p.add_argument('--ny',type=int,default=1)
    p.add_argument('--law',choices=['vg','linear'],default='vg');p.add_argument('--dt',type=float,default=.025)
    p.add_argument('--mode',choices=['ONE_WAY','TWO_WAY'],default='ONE_WAY')
    p.add_argument('--policy',choices=['CONTINUOUS','PRESCRIBED'],default='CONTINUOUS');p.add_argument('--dry',action='store_true')
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    try:
        result=Model(a.nx,a.ny,a.law,a.mode,a.policy,a.dry).run(a.dt)
        write(a.output/'result.json',result)
    except Exception as error:
        write(a.output/'failure.json',dict(complete=False,accepted=None,error=repr(error)));raise
    write(a.output/'execution.json',dict(elapsed_seconds=time.monotonic()-started,complete=True))

if __name__=='__main__':main()
