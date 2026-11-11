"""Observation-only input-built 3D MC plastic bulk + wet Coulomb interface.

Reuses independent geometry, pressure quadrature and Coulomb disk projection.
MC uses projected Gauss-Seidel complementarity (production enumerates active
sets). No production executable, tangent, results or accepted path is an input.
Only fixed hydrostatic CONTINUOUS, all-active affine Tet4/Tet10 boundary cycles.
The finite volume/contact quadrature is prescribed by the input, not a continuum
claim. Partition growth follows this reference's converged-candidate policy;
production's predictor-dependent doubling is deliberately not read from results.
"""
import copy
from types import SimpleNamespace
import numpy as np
import fem_reference as fem
import mc_plastic as mc
from joint_water_reference import Water
from verify_high_order import physical
from verify_joint_coupling import contact
from wet_joint_cycle_reference import contact_work
from hydrostatic_cut_volume_reference import force as cut_force

STAGES = ('geostatic', 'shear', 'open', 'recontact', 'reverse', 'reload', 'return')
RETURNS = ('elastic','shear_face','compression_edge','extension_edge','shear_apex','tension_face','tension_edge','tension_apex','mixed')
CONTACT = ('open','stick','slip')


def integrate(material, old, strain, count):
    state = old
    for _ in range(int(count)):
        state = mc.step(material, state, np.asarray(strain)/count)
    return state


class Model:
    def __init__(self, mesh, ground, case, user_pc=False):
        cfg = ground['initialization']; geo = cfg['geostatic']
        if (cfg['method'] != 'equilibrate' or geo['seed']['kind'] != 'body_states'
                or geo['initial_boundary_factor'] != 0 or case['settings']['kinematics'] != 'explicit'
                or tuple(s['name'] for s in case['stages']) != STAGES[1:]
                or any(s['kind'] != 'boundary' for s in case['stages'])
                or any(case.get(k) for k in ('beam_system','anchor_system','initially_inactive_bodies'))
                or case['joint_system'].get('structure_bodies')
                or any(j.get('water',{}).get('kind') != 'continuous' for j in case['joint_system']['joints'])):
            raise ValueError('outside fixed-hydrostatic MC boundary-cycle scope')
        self.mesh=mesh; self.ground=ground; self.case=case; self.xyz=np.asarray(mesh['nodes']); self.nd=self.xyz.size
        if self.nd > (15000 if user_pc else 300):
            raise ValueError('cloud dense-reference budget exceeded')
        self.water=Water(mesh,ground,case); self.system=self.water.contact
        self.level=cfg['water_level_m']; assert self.level is not None
        self.active=np.ones(len(mesh['tetrahedra']),dtype=bool)
        self.interface=self.water.hydrostatic(self.active,self.level)
        self.ops=[]; seeds=[]; self.b=[]; self.ids=[]; self.weights=[]; self.pore=[]; self.mats=[]; self.positions=[]
        regions={r['id']:r['body'] for r in mesh['regions']}; materials={m['id']:m for m in ground['materials']}
        for ei,mid in enumerate(mesh['material_ids']):
            material=materials[mid]
            if 'mcc' in material or 'invariant_double_hardening' in material:raise ValueError('MC reference only')
            mat=material['mohr_coulomb']; seed=geo['seed']['states'][regions[mesh['region_ids'][ei]]]
            assert seed['equivalent_plastic_strain']==seed['plastic_work_kpa']==0 and not np.any(seed['plastic_strain_compression_positive'])
            q=fem.volume(mesh,ei); ids=(3*q['nodes'][:,None]+np.arange(3)).ravel(); self.ops.append((q,ids))
            for b,w,x in zip(q['b'],q['weight'],q['position']):
                self.b.append(b);self.ids.append(ids);self.weights.append(w);self.pore.append(self.water.gamma*max(self.level-x[2],0));self.positions.append(x)
                self.mats.append(mat);seeds.append(copy.deepcopy(seed))
        self.weights=np.asarray(self.weights);self.pore=np.asarray(self.pore);self.seeds=seeds
        self.pore_force=np.zeros(self.nd)
        for ei,(q,ids) in enumerate(self.ops):self.pore_force[ids]+=cut_force(mesh,ei,self.level,self.water.gamma)["force"].ravel()
        self.seed_force=self.internal(seeds,np.zeros((len(self.system.points),3)))
        base,_=physical(mesh,ground,self.active,cfg['surface_surcharge_kpa'],case.get('pressure_on_original_surface_only',False),joint_system=self.system)
        self.base=(base+self.interface['load']).ravel();self.boundary=fem.conditions_load(mesh,case['boundary'],self.active).ravel()
        basis,error=fem.staged_constraint_space(mesh,case,dict(active_elements=self.active,boundary_factor=0.,displacement_m=np.zeros_like(self.xyz)))
        assert error<1e-12
        self.t=np.zeros((self.nd,sum(b.shape[1] for b in basis.values())));start=0
        for node,b in basis.items():self.t[3*node:3*node+3,start:start+b.shape[1]]=b;start+=b.shape[1]
        self.cr=self.system.c@self.t

    def internal(self,states,traction):
        force=self.pore_force.copy()
        for b,ids,w,pore,s in zip(self.b,self.ids,self.weights,self.pore,states):
            total=np.asarray(s['effective_stress_compression_positive_kpa'])
            force[ids]-=b.T@total*w
        force+=np.einsum('qag,qa,q->g',self.system.c,traction,self.system.weight)
        return force

    def strains(self,u):return np.asarray([b@u[ids] for b,ids in zip(self.b,self.ids)])

    def offset(self,factor):
        equations={};fixed=np.zeros(self.nd)
        for cond in self.case['boundary']['constraints']:
            d=np.asarray(cond['direction']);gradient=np.asarray(cond.get('gradient',[0.,0.,0.]))
            for node in self.mesh['node_sets'][cond['node_set']]:
                a,b=equations.setdefault(node,([],[]));a.append(d);b.append(factor*(cond['value_m']+gradient@(self.xyz[node]-cond.get('reference_point_m',[0.,0.,0.]))))
        for n,(a,b) in equations.items():
            fixed[3*n:3*n+3]=np.linalg.lstsq(a,b,rcond=None)[0]
            np.testing.assert_allclose(np.asarray(a)@fixed[3*n:3*n+3],b,rtol=0,atol=1e-12)
        return fixed

    def fields(self,u,states,jump,traction,plastic,accumulated,pd,nd,work,load,status=None):
        s=self.system;qn=np.einsum('qa,qa->q',traction,s.normal);tau=traction-qn[:,None]*s.normal
        return dict(material_return=np.asarray([RETURNS.index(p['last_return']) for p in states]),contact_status=np.asarray(status if status is not None else np.zeros(len(jump),dtype=int)),u=u.reshape(-1,3),stress=np.asarray([p['effective_stress_compression_positive_kpa'] for p in states]),
                    strain=self.strains(u),plastic_strain=np.asarray([p['plastic_strain_compression_positive'] for p in states]),
                    equivalent_plastic_strain=np.asarray([p['equivalent_plastic_strain'] for p in states]),
                    plastic_work=np.asarray([p['plastic_work_kpa'] for p in states]),void=np.asarray([p['void_ratio'] for p in states]),pore=self.pore.copy(),
                    jump=jump.copy(),traction=traction.copy(),plastic_slip=plastic.copy(),accumulated=accumulated.copy(),
                    gap=np.einsum('qa,qa->q',jump,s.normal),open=np.einsum('qa,qa->q',jump,s.normal)>0,
                    elastic_energy=qn**2/(2*s.kn)+(tau*tau).sum(axis=1)/(2*s.kt),plastic_dissipation=pd.copy(),
                    numerical_dissipation=nd.copy(),endpoint_work=work.copy(),physical_pressure=self.interface['pressure'].copy(),
                    water_load=self.interface['load'].copy(),reaction=(self.internal(states,traction)-load).reshape(-1,3))


def run(model,steps=16,material_maximum=1e-3,contact_maximum=5e-4,sink=None,tolerance=1e-10):
    assert isinstance(steps,int) and steps>=2 and 0<tolerance<=1e-10
    assert 1e-7<=material_maximum<=1e-3 and contact_maximum>=5e-6
    m=model;s=m.system;u=np.zeros(m.nd);states=copy.deepcopy(m.seeds);jump=np.zeros((len(s.points),3));traction=jump.copy()
    plastic=jump.copy();accumulated=np.zeros(len(jump));pd=accumulated.copy();nd=accumulated.copy();work=accumulated.copy()
    rows=[];origin=0.
    for phase,(name,end) in enumerate([('geostatic',0.)]+[(x['name'],x['target_factor']) for x in m.case['stages']]):
        for index in range(1,steps+1):
            progress=index/steps;factor=origin+(end-origin)*progress;load=m.base+(1-progress)*(m.seed_force-m.base) if phase==0 else m.base+factor*m.boundary
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
                    row=dict(stage=name,progress=progress,factor=factor,iterations=iteration,relative_force_residual=relative,
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
        origin=end
    return dict(complete=True,reference_uses_native_results=False,rows=rows,steps=steps,material_maximum=material_maximum,
                contact_maximum=contact_maximum,volume_rule='input Gauss-Jacobi; independent geometry/assembly')
