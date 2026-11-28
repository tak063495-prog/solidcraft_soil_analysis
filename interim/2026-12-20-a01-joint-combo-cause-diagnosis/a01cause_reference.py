"""Opt-in A-01 diagnostic snapshot of frozen 2026-11-23 reference.
Adds input-only contact partition controls and read-only observations.
No production import, native data input, new physics or relaxed tolerance.
"""
"""Input-only small 3D elasticity/Coulomb/structure solve for A-01.

Independent dense SVD constraints, GL125 volume assembly, frame flexibility,
plate quadrature and embedded tendon virtual work. No production solve, saved
state, native accepted fractions, tangent or partitions enters this solver.
The fixed grid is an input. This is a finite-discretization comparison, never
a contact-history convergence or general plastic-bulk certificate.
"""
import copy
from types import SimpleNamespace
import numpy as np
from a01cause_binding import fem
from joint_reference import System
from joint_ext_reference import constraint_space
from beam_reference import System as Frame
from anchor_reference import System as Tendon
from verify_joint_coupling import contact
from wet_joint_cycle_reference import contact_work
from verify_high_order import physical
from mc_reference import normalized_yield

STATES = ('open', 'stick', 'slip')


class Model:
    def __init__(self, mesh, ground, case, user_pc=False):
        self.mesh, self.ground, self.case = mesh, ground, case
        self.xyz = np.asarray(mesh['nodes']); self.nd = self.xyz.size
        if self.nd > (3000 if user_pc else 350):
            raise ValueError('dense independent-reference budget')
        if ground['initialization']['water_level_m'] is not None:
            raise ValueError('wet declared-tip BIRTH remains a production refusal')
        if case.get('initially_inactive_bodies') or any(d.get('crack_tip') is not None or d.get('water') is not None for d in case['joint_system']['joints']):
            raise ValueError('this fresh solve covers dry non-tip all-active contact only')
        self.system = System(mesh, case['joint_system'])
        self.frame = Frame(mesh, case['beam_system']) if case.get('beam_system') else None
        self.tendon = Tendon(mesh, case['anchor_system']) if case.get('anchor_system') else None
        self.active = np.ones(len(mesh['tetrahedra']), dtype=bool)
        self.k = np.zeros((self.nd, self.nd)); self.seed = np.zeros(self.nd); self.ops = []
        materials = {m['id']: m['mohr_coulomb'] for m in ground['materials']}
        regions = {r['id']: r['body'] for r in mesh['regions']}
        geo = ground['initialization']['geostatic']
        assert geo['initial_boundary_factor']==0 and geo['seed']['kind']=='body_states'
        for e, mid in enumerate(mesh['material_ids']):
            mat = materials[mid]; d = fem.elastic_tensor(mat['young_kpa'], mat['poisson'])
            state = geo['seed']['states'][regions[mesh['region_ids'][e]]]
            assert state['equivalent_plastic_strain']==state['plastic_work_kpa']==0 and not np.any(state['plastic_strain_compression_positive'])
            s0 = np.asarray(state['effective_stress_compression_positive_kpa'])
            op = fem.volume(mesh,e,independent=True,independent_order=5)
            ids = (3*op['nodes'][:,None]+np.arange(3)).ravel()
            self.k[np.ix_(ids,ids)] += np.einsum('qai,ab,qbj,q->ij',op['b'],d,op['b'],op['weight'],optimize=True)
            self.seed[ids] -= np.einsum('qai,a,q->i',op['b'],s0,op['weight'])
            self.ops.append((fem.volume(mesh,e),ids,d,s0,state['void_ratio'],mat))
        load,_ = physical(mesh,ground,self.active,ground['initialization']['surface_surcharge_kpa'],case['pressure_on_original_surface_only'],joint_system=self.system)
        self.base = load.ravel()
        self.boundary = fem.conditions_load(mesh,case['boundary'],self.active).ravel()
        self.spaces = {tied: constraint_space(mesh,case,self.active,self.system,tied) for tied in (False,True)}
        self.index = np.array([next(i for i,d in enumerate(case['joint_system']['joints']) if d['id']==p['joint_id']) for p in self.system.points])

    def fields(self, u, joint, born, generations, archived, load, internal, frame_result, tendon_result, holding, frozen):
        s=self.system;stress=[];strain=[];void=[]
        for op,ids,d,s0,e0,mat in self.ops:
            eps=np.einsum('qai,i->qa',op['b'],u[ids]);sig=s0-eps@d.T
            if np.max(normalized_yield(mat,sig))>1e-8:
                raise ValueError('fresh elastic reference reached bulk yield; plastic generalization is outside A-01 scope')
            stress.extend(sig);strain.extend(eps);void.extend((1+e0)*np.exp(eps[:,:3].sum(axis=1))-1)
        q=joint['traction'];qn=np.einsum('qa,qa->q',q,s.normal);tau=q-qn[:,None]*s.normal
        gap=np.einsum('qa,qa->q',joint['jump']-joint['rest'],s.normal)
        result=dict(u=u.reshape(-1,3), stress=stress, strain=strain, void=void,
                    plastic_strain=np.zeros_like(stress), equivalent_plastic_strain=np.zeros(len(stress)),
                    plastic_work=np.zeros(len(stress)), material_return=np.zeros(len(stress),dtype=int),
                    jump=joint['jump'], traction=q, plastic_slip=joint['plastic'], accumulated=joint['accumulated'],
                    plastic_dissipation=joint['pd'], numerical_dissipation=joint['nd'], endpoint_work=joint['work'],
                    elastic_energy=qn**2/(2*s.kn)+(tau*tau).sum(axis=1)/(2*s.kt), gap=gap,
                    opening=np.maximum(gap,0), contact_open=gap>0, contact_status=joint['status'],
                    initial_traction=joint['initial'], reference_jump=joint['reference'], elastic_offset=joint['rest'],
                    initial_energy=joint['initial_energy'], generation=np.asarray(generations),
                    joint_active=np.full(len(generations),born), joint_tied=np.full(len(generations),not born),
                    archived_generation_count=np.asarray(archived),
                    internal=internal.reshape(-1,3), applied_load=load.reshape(-1,3),
                    reaction=(internal-load).reshape(-1,3), holding_load=holding.reshape(-1,3),
                    normal_pressure=-qn, shear_traction=tau)
        for key,value in frozen.items():result['frozen_'+key]=value
        if frame_result is not None:
            result['structure_u']=frame_result['state']['node_displacement_m']
            result['structure_rotation']=frame_result['state']['node_rotation_rad']
            result['structure_support_force']=frame_result['support_force']
            result['structure_support_moment']=frame_result['support_moment']
            result['structure_external_force']=frame_result['external_force']
            result['structure_external_moment']=frame_result['external_moment']
            for key, member in frame_result['members'].items():
                for name in ('local_end_actions','local_end_relative_rotation_rad','section_forces_at_start_mid_end'):
                    result['beam_'+key+'_'+name]=member[name]
            for key,p in frame_result.get('plates',{}).items():
                for name in ('membrane_energy_kn_m','bending_energy_kn_m','drilling_penalty_energy_kn_m'):
                    result['plate_'+key+'_'+name]=p[name]
                for name in p['samples'][0]:
                    result['plate_'+key+'_'+name]=[v[name] for v in p['samples']]
                # This specimen uses the existing three-node DKT/CST thin
                # plate. Kirchhoff transverse shear is identically zero;
                # it is not a missing Mindlin shear recovery or a fitted value.
                result['plate_'+key+'_transverse_shear_strain']=np.zeros((len(p['samples']),2))
                result['plate_'+key+'_transverse_shear_force_kn_m']=np.zeros((len(p['samples']),2))
        if tendon_result is not None:
            state,_,_,reports=tendon_result
            for i,r in enumerate(reports):
                for name in ('point_displacement_m','point_soil_force_kn','segment_forces_kn','segment_extension_m','elastic_energy_kn_m'):
                    result['anchor_'+str(i)+'_'+name]=r[name]
                result['anchor_'+str(i)+'_taut']=state['anchors'][i]['segment_taut']
        return {k:np.asarray(v) for k,v in result.items()}


def run(m,steps,sink,contact_maximum=None,*,fixed_parts=None,boost_nontrivial=False,observer=None):
    if fixed_parts is not None and (type(fixed_parts) is not int or not 1<=fixed_parts<=64):
        raise ValueError('a01cause fixed parts outside 1..64')
    assert steps>=2
    def observe(kind,**data):
        if observer is not None:observer(dict(kind=kind,**data))
    s=m.system;nq=len(s.points);nd=m.nd
    maximum=contact_maximum or m.case['joint_system']['maximum_increment_m']
    u=np.zeros(nd);born=False;generations=[0]*len(m.case['joint_system']['joints']);archived=[0]*len(generations)
    joint={k:np.zeros((nq,3)) for k in ('jump','traction','plastic','initial','reference','rest')}
    joint.update({k:np.zeros(nq) for k in ('accumulated','pd','nd','work','initial_energy')});joint['status']=np.ones(nq,dtype=int)
    frozen={k:np.zeros_like(v) for k,v in joint.items()}
    bs=m.frame.initial_state() if m.frame else None;ats=m.tendon.initial_state() if m.tendon else None
    rows=[];factor=0.;last_load=None;max_residual=0.

    def structure_at(candidate,factor,bs):
        br=m.frame.solve_interface(bs,candidate.reshape(-1,3),factor) if m.frame else None
        ar=m.tendon.evaluate(ats,candidate.reshape(-1,3)) if m.tendon else None
        force=np.zeros(nd);attached=np.zeros(nd)
        if br:
            for node,v in br['interface_residual'].items(): force[3*node:3*node+3]+=v
            for node,v in br['attached_load'].items(): attached[3*node:3*node+3]+=v
        if ar:force+=ar[1].ravel()
        return force,attached,br,ar

    # Structure interface solve above is used for physical nodal conventions
    # and birth forces. Newton below retains all structure free rotations.
    for phase,stage in enumerate([dict(kind='geostatic',name='geostatic')]+m.case['stages']):
        kind=stage['kind'];origin=factor;end=stage.get('target_factor',factor)
        if kind=='joint_birth':
            assert not born
            ref=np.einsum('qag,g->qa',s.c,u);initial=np.zeros_like(ref)
            if all(d['birth_policy']=='transfer_force' for d in m.case['joint_system']['joints']):
                free,_=m.spaces[False]
                sf,attached,_,_=structure_at(u,factor,bs)
                rhs=last_load-(m.k@u+m.seed+sf+attached)
                pairs=sorted({(a,b) for p in s.points for a,b in zip(p['nodes_a'],p['nodes_b']) if a!=b})
                pm=np.zeros((nd,3*len(pairs)))
                for i,(a,b) in enumerate(pairs):
                    pm[3*a:3*a+3,3*i:3*i+3]=-np.eye(3);pm[3*b:3*b+3,3*i:3*i+3]=np.eye(3)
                multipliers=np.linalg.lstsq(free.T@pm,free.T@rhs,rcond=1e-12)[0];target=pm@multipliers
                if np.max(np.abs(free.T@(target-rhs)),initial=0)>1e-10:raise ArithmeticError('independent nodal multiplier residual')
                a=np.einsum('qag,q->gqa',s.c,np.sqrt(s.weight)).reshape(nd,-1)
                projection=np.linalg.lstsq(a,target,rcond=1e-12)[0]
                if np.max(np.abs(a@projection-target),initial=0)>1e-10:raise ArithmeticError('independent GP transfer residual')
                initial=projection.reshape(-1,3)/np.sqrt(s.weight)[:,None]
                pressure=-np.einsum('qa,qa->q',s.normal,initial)
                if pressure.min()<-1e-8 or np.any(np.linalg.norm(initial+pressure[:,None]*s.normal,axis=1)>s.mu*pressure+1e-8):
                    raise ArithmeticError('independent initial transfer outside Coulomb cone')
            for i,g in enumerate(generations):
                if g and archived[i]<g:archived[i]+=1
                generations[i]+=1
            qn=np.einsum('qa,qa->q',initial,s.normal);tau=initial-qn[:,None]*s.normal
            joint={k:np.zeros_like(v) for k,v in joint.items()}
            joint.update(jump=ref.copy(),traction=initial.copy(),reference=ref.copy(),initial=initial.copy(),
                         rest=ref-qn[:,None]*s.normal/s.kn[:,None]-tau/s.kt[:,None],
                         initial_energy=qn**2/(2*s.kn)+(tau*tau).sum(axis=1)/(2*s.kt))
            joint['status']=np.ones(nq,dtype=int)
            born=True
        elif kind=='joint_death':
            assert born;born=False
            frozen=copy.deepcopy(joint)
            for i in range(len(archived)):archived[i]+=1
        t,offset=m.spaces[not born]
        holding=np.zeros(nd)
        sf,attached,_,_=structure_at(u,factor,bs)
        initial_internal=m.k@u+m.seed+sf+attached
        if born:initial_internal+=np.einsum('qag,qa,q->g',s.c,joint['traction'],s.weight)
        base=m.base+factor*m.boundary+attached
        if phase==0:holding=initial_internal-base
        elif kind in ('joint_birth','joint_death'):
            holding=initial_internal-base
            if kind=='joint_birth' and all(d['birth_policy']=='transfer_force' for d in m.case['joint_system']['joints']):holding=np.zeros(nd)
        for index in range(1,steps+1):
            progress=index/steps;factor=origin+(end-origin)*progress
            fixed=factor*offset
            # Build a combined soil/frame kinematic transformation. The frame
            # free DOFs remain explicit; production uses a condensed interface.
            if m.frame:
                kb,rest,fb,eff,extf,extm,info=m.frame.assemble(bs,factor)
                bo,binding,qb,normals=m.frame.space(bs,fixed.reshape(-1,3),soil_basis=t)
                for i,node in enumerate(m.frame.attach):
                    if node is not None:bo[6*i:6*i+3]+=fixed[3*node:3*node+3]
                transform=np.column_stack((binding,qb));soil_t=np.column_stack((t,np.zeros((nd,qb.shape[1]))))
                matrix=soil_t.T@m.k@soil_t+transform.T@kb@transform
                constant=soil_t.T@(m.k@fixed+m.seed)+transform.T@(kb@bo+rest-fb)
                free=np.r_[t.T@(u-fixed),np.zeros(qb.shape[1])]
            else:
                soil_t=t;matrix=t.T@m.k@t;constant=t.T@(m.k@fixed+m.seed);free=t.T@(u-fixed)
            physical_load=m.base+factor*m.boundary
            external=physical_load+(1-progress)*holding
            # Attached frame loads cancel from the reduced equations because
            # frame residual includes -fb. Holding is a physical nodal vector.
            constant-=soil_t.T@external
            counts=np.full(nq,fixed_parts if fixed_parts is not None else 1,dtype=int)
            def projected(target):
                q=np.zeros((nq,1,3) if target.ndim==2 else target.shape);dp=q.copy();da=np.zeros(q.shape[:2])
                for count in np.unique(counts):
                    mask=counts==count;sub=SimpleNamespace(**{key:getattr(s,key)[mask] for key in ('normal','kn','kt','mu')})
                    shifted=target[mask]-joint['rest'][mask] if target.ndim==2 else target[mask]-joint['rest'][mask,None,:]
                    q[mask],dp[mask],da[mask]=contact(sub,joint['jump'][mask]-joint['rest'][mask],joint['traction'][mask],shifted,int(count))
                return q[:,0],dp[:,0],da[:,0]
            cr=s.c@soil_t
            def evaluate(x,tangent=False):
                candidate=fixed+soil_t@x;target=np.einsum('qag,g->qa',s.c,candidate)
                q,dp,da=projected(target) if born else (np.zeros((nq,3)),None,None)
                residual=matrix@x+constant
                if born:residual+=np.einsum('qai,qa,q->i',cr,q,s.weight)
                ar=m.tendon.evaluate(ats,candidate.reshape(-1,3)) if m.tendon else None
                if ar:residual+=soil_t.T@ar[1].ravel()
                if not tangent:return residual,candidate,target,q,dp,da,ar
                jac=matrix.copy()
                if born:
                    h=2e-9;d=np.zeros((nq,3,3))
                    for a in range(3):
                        e=np.zeros(3);e[a]=h
                        d[:,:,a]=(projected(target+e)[0]-projected(target-e)[0])/(2*h)
                    jac+=np.einsum('qai,qab,qbj,q->ij',cr,d,cr,s.weight,optimize=True)
                if ar:jac+=soil_t.T@ar[2]@soil_t
                return residual,jac
            scale=max(np.linalg.norm(external),1.)
            for iteration in range(100):
                residual,candidate,target,q,dp,da,ar=evaluate(free);relative=np.linalg.norm(residual)/scale
                if relative<1e-10:
                    needed=np.maximum(1,np.ceil(np.linalg.norm(target-joint['jump'],axis=1)/maximum).astype(int))
                    if boost_nontrivial:needed=np.where(needed>1,2*needed,needed)
                    if born and fixed_parts is None and np.any(needed>counts):counts=np.maximum(counts,needed);continue
                    break
                residual,jac=evaluate(free,True)
                observe('linear',stage=stage['name'],progress=progress,iteration=iteration,
                        matrix=jac.copy(),rhs=-residual.copy(),free=free.copy(),soil_t=soil_t.copy(),
                        target=target.copy(),gap=np.einsum('qa,qa->q',target-joint['rest'],s.normal),
                        traction=q.copy(),contact_partitions=counts.copy(),relative=float(relative))
                correction=np.linalg.solve(jac,-residual)
                for cut in range(18):
                    trial=free+correction*2.**(-cut)
                    trial_norm=float(np.linalg.norm(evaluate(trial)[0]))
                    observe('line_search',stage=stage['name'],progress=progress,iteration=iteration,cut=cut,residual_norm=trial_norm)
                    if trial_norm<np.linalg.norm(residual):free=trial;break
                else:raise ArithmeticError(('independent line search',stage['name'],progress,float(relative)))
            else:raise ArithmeticError(('independent Newton limit',stage['name'],progress,float(relative)))
            if born:
                status=np.zeros(nq,dtype=int)
                for count in np.unique(counts):
                    mask=counts==count;sub=SimpleNamespace(**{key:getattr(s,key)[mask] for key in ('normal','kn','kt','mu')})
                    before=joint['jump'][mask]-joint['rest'][mask];after=target[mask]-joint['rest'][mask]
                    pd,num,work=contact_work(sub,before,joint['traction'][mask],after,int(count))
                    joint['pd'][mask]+=pd;joint['nd'][mask]+=num;joint['work'][mask]+=work
                    if count>1:
                        penultimate=before+(after-before)*(count-1)/count
                        priorq=contact(sub,before,joint['traction'][mask],penultimate,int(count-1))[0][:,0]
                    else:penultimate=before;priorq=joint['traction'][mask]
                    lastflow=contact(sub,penultimate,priorq,after,1)[1][:,0]
                    gap=np.einsum('qa,qa->q',after,sub.normal)
                    status[mask]=np.where(gap>0,0,np.where(np.linalg.norm(lastflow,axis=1)>0,2,1))
                joint.update(jump=target,traction=q,plastic=joint['plastic']+dp,accumulated=joint['accumulated']+da,status=status)
            u=candidate
            br=None;attached=np.zeros(nd);sf=np.zeros(nd)
            if m.frame:
                br=m.frame.results(bs,bo+transform@free,kb,rest,fb,eff,extf,extm,info,normals);bs=br['state']
                for node,v in br['interface_residual'].items():sf[3*node:3*node+3]+=v
                for node,v in br['attached_load'].items():attached[3*node:3*node+3]+=v
            if ar:ats=ar[0];sf+=ar[1].ravel()
            internal=m.k@u+m.seed+sf+attached
            if born:internal+=np.einsum('qag,qa,q->g',s.c,joint['traction'],s.weight)
            last_load=external+attached
            row=dict(stage=stage['name'],stage_index=phase,progress=progress,factor=factor,relative_force_residual=float(relative),contact_partitions=counts.tolist())
            max_residual=max(max_residual,float(relative));rows.append(row)
            sink(row,m.fields(u,joint,born,generations,archived,last_load,internal,br,ar,(1-progress)*holding,frozen))
    return dict(complete=True,reference_uses_native_results=False,accepted_states=len(rows),rows=rows,maximum_force_residual=max_residual)
