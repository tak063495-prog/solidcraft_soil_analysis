"""Input-only 3D reference for prescribed changing hydrostatic water.

The C-03 independent bulk/Coulomb solve is reused without modifying its source.
This observation extension recomputes cut-volume, exterior and interface water
forces at every specified level. No native displacement, stress, partition or
accepted schedule is used to generate the independent solution. Stage clock is
an explicitly supplied reporting clock, not a new transient storage equation.
"""
import copy
import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'tools'))
C03=HERE.parent/'2026-11-20-wet-joint-plastic-bulk-integration'
sys.path.insert(0,str(C03))
spec=importlib.util.spec_from_file_location('c04_c03_reference',C03/'reference.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
from verify_high_order import physical
from hydrostatic_cut_volume_reference import force as cut_force
from verify_joint_coupling import contact
from wet_joint_cycle_reference import contact_work
integrate=old.integrate


class Model(old.Model):
    def __init__(self,mesh,ground,case,user_pc=False):
        if any(s['kind'] not in ('boundary','hydrostatic_water') for s in case['stages']):
            raise ValueError('only prescribed hydrostatic/boundary stages; Richards JOINT unsupported')
        if ground['initialization'].get('water_level_m') is None:
            raise ValueError('explicit initial hydrostatic level required')
        fixed=copy.deepcopy(case)
        fixed['stages']=[s for s in fixed['stages'] if s['kind']=='boundary']
        super().__init__(mesh,ground,fixed,user_pc)
        self.case=case
        self.dry_base=physical(mesh,{**ground,'initialization':{**ground['initialization'],'water_level_m':None}},self.active,ground['initialization']['surface_surcharge_kpa'],case.get('pressure_on_original_surface_only',False),joint_system=self.system)[0]
        self.set_level(self.level)

    def set_level(self,level):
        self.level=float(level)
        self.pore=self.water.gamma*np.maximum(self.level-np.asarray(self.positions)[:,2],0.)
        self.pore_force=np.zeros(self.nd)
        for ei,(q,ids) in enumerate(self.ops):
            self.pore_force[ids]+=cut_force(self.mesh,ei,self.level,self.water.gamma)['force'].ravel()
        self.interface=self.water.hydrostatic(self.active,self.level)
        frame_ground={**self.ground,'initialization':{**self.ground['initialization'],'water_level_m':self.level}}
        loads,_=physical(self.mesh,frame_ground,self.active,self.ground['initialization']['surface_surcharge_kpa'],self.case.get('pressure_on_original_surface_only',False),joint_system=self.system)
        self.exterior_water=loads-self.dry_base
        self.base=(loads+self.interface['load']).ravel()

    def fields(self,*args,**kwargs):
        f=super().fields(*args,**kwargs)
        f['total_stress']=f['stress'].copy();f['total_stress'][:,:3]+=self.pore[:,None]
        f['bulk_water_force']=self.pore_force.reshape(-1,3).copy()
        f['exterior_water_load']=self.exterior_water.copy()
        f['applied_load']=np.asarray(args[9]).reshape(-1,3).copy()
        normal=self.system.normal
        f['effective_contact_pressure']=-np.einsum('qa,qa->q',f['traction'],normal)
        f['total_normal_pressure']=f['effective_contact_pressure']+f['physical_pressure']
        f['friction_limit']=np.maximum(f['effective_contact_pressure'],0)*self.system.mu
        f['nodal_prescribed_level']=np.full(len(self.xyz),self.level)
        f['joint_head_a']=np.full(len(normal),self.level)
        f['joint_head_b']=f['joint_head_a'].copy()
        return f


def run(model,steps=16,material_maximum=1e-3,contact_maximum=5e-4,sink=None,tolerance=1e-10):
    assert isinstance(steps,int) and steps>=2 and 0<tolerance<=1e-10
    assert 1e-7<=material_maximum<=1e-3 and contact_maximum>=5e-6
    m=model;s=m.system;u=np.zeros(m.nd);states=copy.deepcopy(m.seeds);jump=np.zeros((len(s.points),3));traction=jump.copy()
    plastic=jump.copy();accumulated=np.zeros(len(jump));pd=accumulated.copy();nd=accumulated.copy();work=accumulated.copy()
    rows=[];origin=0.
    level_origin=m.level
    phases=[dict(name='geostatic',kind='initial',target_factor=0.,initial_increment=1/steps)]+m.case['stages']
    for phase,stage in enumerate(phases):
        name=stage['name'];end=stage.get('target_factor',origin)
        count=round(1/stage['initial_increment'])
        for index in range(1,count+1):
            progress=index/count;factor=origin+(end-origin)*progress
            level=level_origin+(stage.get('target_level_m',level_origin)-level_origin)*progress
            m.set_level(level)
            load=m.base+(1-progress)*(m.seed_force-m.base) if phase==0 else m.base+factor*m.boundary
            fixed=m.offset(factor);free=m.t.T@(u-fixed);scale=max(np.linalg.norm(load),1.)
            counts=np.ones(len(states),dtype=int);jc=np.ones(len(jump),dtype=int)
            def required(candidate,target):
                de=-m.strains(candidate-u)
                return np.maximum(1,np.ceil(np.abs(de).sum(axis=1)/material_maximum).astype(int)),np.maximum(1,np.ceil(np.linalg.norm(target-jump,axis=1)/contact_maximum).astype(int))
            def projected(target):
                q=np.zeros((len(jump),1,3) if target.ndim==2 else target.shape);dp=q.copy();da=np.zeros(q.shape[:2])
                for c in np.unique(jc):
                    mask=jc==c;sub=SimpleNamespace(**{k:getattr(s,k)[mask] for k in ('normal','kn','kt','mu')})
                    q[mask],dp[mask],da[mask]=contact(sub,jump[mask],traction[mask],target[mask],int(c))
                return q,dp,da
            def evaluate(x,tangent=False):
                candidate=fixed+m.t@x;target=np.einsum('qag,g->qa',s.c,candidate);de=-m.strains(candidate-u)
                updated=[integrate(mat,old,e,c) for mat,old,e,c in zip(m.mats,states,de,counts)]
                q,dp,da=projected(target);res=m.t.T@(m.internal(updated,q[:,0])-load)
                if not tangent:return res,candidate,target,updated,q[:,0],dp[:,0],da[:,0]
                jac=np.zeros((m.t.shape[1],)*2);h=1e-8
                for b,ids,w,mat,old,e,c in zip(m.b,m.ids,m.weights,m.mats,states,de,counts):
                    d=np.empty((6,6))
                    for a in range(6):
                        perturb=np.zeros(6);perturb[a]=h
                        sp=integrate(mat,old,e+perturb,c)['effective_stress_compression_positive_kpa'];sm=integrate(mat,old,e-perturb,c)['effective_stress_compression_positive_kpa']
                        d[:,a]=(np.asarray(sp)-sm)/(2*h)
                    br=b@m.t[ids];jac+=br.T@d@br*w
                h=2e-9;targets=np.repeat(target[:,None,:],7,axis=1)
                for a in range(3):targets[:,1+2*a,a]+=h;targets[:,2+2*a,a]-=h
                pq,_,_=projected(targets);d=np.stack([(pq[:,1+2*a]-pq[:,2+2*a])/(2*h) for a in range(3)],axis=2)
                jac+=np.einsum('qai,qab,qbj,q->ij',m.cr,d,m.cr,s.weight,optimize=True)
                return res,jac
            for iteration in range(100):
                res,candidate,target,updated,q,dp,da=evaluate(free);relative=float(np.linalg.norm(res)/scale)
                if relative<tolerance:
                    needed,jneeded=required(candidate,target)
                    if np.any(needed>counts) or np.any(jneeded>jc):counts=np.maximum(counts,needed);jc=np.maximum(jc,jneeded);continue
                    status=np.empty(len(jump),dtype=int)
                    for c in np.unique(jc):
                        mask=jc==c;sub=SimpleNamespace(**{k:getattr(s,k)[mask] for k in ('normal','kn','kt','mu')})
                        a,b,cw=contact_work(sub,jump[mask],traction[mask],target[mask],int(c));pd[mask]+=a;nd[mask]+=b;work[mask]+=cw
                        before=jump[mask];before_q=traction[mask]
                        if c>1:
                            before=jump[mask]+(target[mask]-jump[mask])*(c-1)/c
                            before_q=contact(sub,jump[mask],traction[mask],before,int(c-1))[0][:,0]
                        last_flow=contact(sub,before,before_q,target[mask],1)[1][:,0]
                        gap=(target[mask]*sub.normal).sum(axis=1)
                        status[mask]=np.where(gap>0,0,np.where(np.linalg.norm(last_flow,axis=1)>0,2,1))
                    u=candidate;states=updated;jump=target;traction=q;plastic+=dp;accumulated+=da
                    row=dict(stage=name,progress=progress,factor=factor,level_m=m.level,time_s=phase+progress,iterations=iteration,relative_force_residual=relative,
                             material_partitions=counts.tolist(),contact_partitions=jc.tolist())
                    rows.append(row)
                    if sink:sink(row,m.fields(u,states,jump,traction,plastic,accumulated,pd,nd,work,load,status))
                    break
                res,jac=evaluate(free,True);correction=np.linalg.solve(jac,-res)
                for refinement in range(5):
                    linear=float(np.linalg.norm(jac@correction+res)/max(np.linalg.norm(res),1e-30))
                    if linear<=1e-10:break
                    correction+=np.linalg.solve(jac,-res-jac@correction)
                else:raise ArithmeticError(('true linear residual',linear))
                for cut in range(18):
                    trial=free+correction*2.**(-cut)
                    try:trial_res=evaluate(trial)[0]
                    except ArithmeticError:continue
                    if np.linalg.norm(trial_res)<np.linalg.norm(res):free=trial;break
                else:raise ArithmeticError(('independent line search',name,progress,iteration,relative))
            else:raise ArithmeticError(('independent Newton limit',name,progress,relative))
        origin=end;level_origin=m.level
    return dict(complete=True,reference_uses_native_results=False,rows=rows,steps=steps,material_maximum=material_maximum,
                contact_maximum=contact_maximum,volume_rule='input Gauss-Jacobi; independent geometry/assembly')
