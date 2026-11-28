"""Fresh branch tangents, localization, and independent descent decomposition."""
import argparse, copy, itertools, json
from pathlib import Path
from types import SimpleNamespace
import subprocess
import numpy as np
from a01follow_binding import fixtures, native_events, matrix, assess, rigid_fit, norm, dot
from a01follow_reference import Model
from a01follow_contact import evaluate as analytic, assemble
from a01follow_run import write
from verify_joint_coupling import contact

def events(path):return [json.loads(line) for line in path.read_text(encoding='utf-8').splitlines()]

def joint_from_native(points):
    result={name:np.asarray([p[key] for p in points]) for name,key in
        dict(jump='jump_m',traction='traction_kpa').items()}
    result['rest']=np.asarray([p.get('birth',{}).get('elastic_offset_m',[0.,0.,0.]) for p in points])
    return result

def projected(s,joint,target,counts):
    q=np.zeros_like(target)
    for count in np.unique(counts):
        mask=counts==count;sub=SimpleNamespace(**{k:getattr(s,k)[mask] for k in ('normal','kn','kt','mu')})
        q[mask]=contact(sub,joint['jump'][mask]-joint['rest'][mask],joint['traction'][mask],target[mask]-joint['rest'][mask],int(count))[0][:,0]
    return q

def fd(s,joint,target,counts,h):
    d=np.zeros((len(target),3,3))
    for a in range(3):
        step=np.eye(3)[a]*h
        d[:,:,a]=(projected(s,joint,target+step,counts)-projected(s,joint,target-step,counts))/(2*h)
    return d

def vector_record(mesh,t,v,a,parts):
    physical=(t@v).reshape(-1,3);xyz=np.asarray(mesh['nodes'])
    central=np.flatnonzero((xyz[:,2]>=1/3-1e-14)&(xyz[:,2]<=2/3+1e-14))
    # The split mesh has coincident nodes from adjacent blocks. Choose the
    # actual middle face-connected component, not a coordinate slab.
    cells=[cell for cell in mesh['tetrahedra'] if 1/3<xyz[cell,2].mean()<2/3]
    central=np.asarray(sorted(set(itertools.chain.from_iterable(cells))))
    x=xyz[central]-xyz[central].mean(axis=0)
    basis=np.zeros((len(xyz)*3,6))
    for node,r in zip(central,x):basis[3*node:3*node+3]=np.c_[np.eye(3),np.cross(np.eye(3),r).T]
    reduced=t.T@basis;q,_=np.linalg.qr(reduced)
    return dict(vector=v.tolist(),rigid_fits=rigid_fit(mesh,t,[v])[0],
        middle_motion_fraction=norm(physical[central])**2/max(norm(physical)**2,1e-300),
        middle_rigid_subspace_overlap_squared=norm(q.T@v)**2/max(norm(v)**2,1e-300),
        total_action_norm=norm(a@v),
        contributions={key:dict(action_norm=norm(k@v),quadratic_form=dot(v,k@v),
            action=(k@v).tolist()) for key,k in parts.items()},
        action_sum_error=norm(sum((k@v for k in parts.values()),np.zeros(len(v)))-a@v))

def modes(folder,out,probe):
    m=Model(*(fixtures.read(folder/(key+'.json')) for key in ('mesh','ground','case')))
    trace=events(folder/'native/progress.jsonl');e=native_events(folder/'native/progress.jsonl')[-1]
    original=e['original_tangent'];a0=np.asarray(original['matrix_columns']).T
    it=[v for v in trace if v['kind']=='iterate'][-1]
    t=np.zeros((m.nd,len(a0)))
    for i,row in enumerate(original['dof_rows']):
        for j,c in row:t[i,j]=c
    old=it['committed']['joints']['points'];joint=joint_from_native(old)
    target=np.einsum('qag,g->qa',m.system.c,np.asarray(it['displacement_m']).ravel())
    counts=np.asarray(it['partitions'][len(it['committed']['material_states']):],dtype=int)
    h=m.case['joint_system']['tangent_step_m']
    kernel_input=dict(mesh=m.mesh,configuration=m.case['joint_system'],old_state=it['committed']['joints'],
        displacement=it['displacement_m'],counts=counts.tolist())
    write(out/'native-kernel-input.json',kernel_input)
    process=subprocess.run([str(probe),str(out/'native-kernel-input.json'),str(out/'native-kernel-output.json')],check=True)
    kernel=fixtures.read(out/'native-kernel-output.json');raw=kernel['records']
    df=np.asarray([v['tangent_columns'] for v in raw]).transpose(0,2,1)
    target=np.asarray([v['target'] for v in kernel['local_inputs']])
    native_c=np.zeros_like(m.system.c)
    for i,coefficients in enumerate(kernel['jump_coefficients']):
        for node,c in coefficients:
            for axis in range(3):native_c[i,axis,3*node+axis]=c
    cr=native_c@t
    m.system.normal=np.asarray([v['normal'] for v in kernel['layout']])
    m.system.weight=np.asarray([v['weight_m2'] for v in kernel['layout']])
    q,da,branches=analytic(m.system,joint,target,counts)
    native_contact=assemble(cr,df,m.system.weight)
    # Keep the actual stopped global matrix; only replace the re-evaluated
    # contact derivative. This retains its original bulk and assembly roundoff.
    skeleton=a0-native_contact
    independent_bulk=t.T@m.k@t
    mats={'analytic':da,'native_central_fd':df,'independent_central_fd':fd(m.system,joint,target,counts,h),
        'central_fd_half_h':fd(m.system,joint,target,counts,h/2),
        'analytic_contact_and_elastic_bulk':da,'native_fd_with_analytic_elastic_bulk':df}
    _,sv,vh=np.linalg.svd(a0);records={}
    branch_indices={key:[i for i,b in enumerate(branches) if b['branch']==key] for key in ('open','stick','slip')}
    xyz=np.asarray(m.mesh['nodes']);cells=[c for c in m.mesh['tetrahedra'] if 1/3<xyz[c,2].mean()<2/3]
    middle=sorted(set(itertools.chain.from_iterable(cells)));center=xyz[middle].mean(axis=0)
    rb=np.zeros((m.nd,6))
    for node in middle:rb[3*node:3*node+3]=np.c_[np.eye(3),np.cross(np.eye(3),xyz[node]-center).T]
    rigid,_=np.linalg.qr(t.T@rb)
    for name,d in mats.items():
        selected_bulk=independent_bulk if name in ('analytic_contact_and_elastic_bulk','native_fd_with_analytic_elastic_bulk') else skeleton
        parts={'bulk':selected_bulk}
        for key,indices in branch_indices.items():
            mask=np.zeros(len(target),dtype=bool);mask[indices]=True
            parts[key]=assemble(cr,d,m.system.weight,mask)
        total=sum(parts.values(),np.zeros_like(skeleton));u,s,right=np.linalg.svd(total)
        contact_rigid=(total-selected_bulk)@rigid;csv=np.linalg.svd(contact_rigid,compute_uv=False)
        records[name]=dict(smallest_singular_values=s[-8:].tolist(),
            numerical_rank=int(np.sum(s>len(s)*np.finfo(float).eps*s[0])),
            source_pre_predictor_matrix_difference=norm(total-a0),
            middle_rigid_contact_singular_values=csv.tolist(),
            middle_rigid_contact_rank=int(np.sum(csv>max(contact_rigid.shape)*np.finfo(float).eps*csv[0])),
            per_joint_contact_action={jid:norm(assemble(cr,d,m.system.weight,np.asarray([p['joint_id']==jid for p in m.system.points]))@rigid)
                for jid in sorted({p['joint_id'] for p in m.system.points})},
            source_two_directions=[vector_record(m.mesh,t,v,total,parts) for v in vh[-2:]],
            own_two_directions=[vector_record(m.mesh,t,v,total,parts) for v in right[-2:]],
            matrices={key:value.tolist() for key,value in parts.items()})
    qnative=np.asarray([v['state']['traction_kpa'] for v in raw])
    write(out/'mode-separation.json',dict(passed=False,diagnostic_only=True,
        context=e['context'],source_pre_predictor_smallest_singular_values=sv[-8:].tolist(),
        source_directions=vh[-2:].tolist(),branches=branches,branch_gp_indices=branch_indices,
        captured_bulk_plus_assembly_roundoff_vs_independent_bulk=norm(skeleton-independent_bulk),
        native_jump_coefficients_vs_independent_maximum=float(np.max(abs(native_c-m.system.c))),
        counts=counts.tolist(),target=target.tolist(),analytic_traction_vs_native_maximum=float(np.max(abs(q-qnative))),
        gp_layout=kernel['layout'],
        native_vs_independent_fd_maximum=float(np.max(abs(df-mats['independent_central_fd']))),
        local_tangent_analytic_vs_fd_maximum=float(np.max(abs(da-df))),
        linearizations=records,existing_predictor_change=(matrix(e)-a0).tolist(),
        existing_predictor_change_action_norms=[norm((matrix(e)-a0)@v) for v in vh[-2:]],
        new_support_or_regularization=False,physical_instability_proved=False))

def direction_metrics(residual,jac,correct_jac,evaluate,x):
    result=dict(rank=int(np.linalg.matrix_rank(jac)),condition=float(np.linalg.cond(jac)),
        residual_norm=norm(residual))
    try:
        p=np.linalg.solve(jac,-residual)
        candidates=[norm(evaluate(x+p*2.**(-cut))) for cut in range(18)]
        result.update(direction=p.tolist(),true_linear_relative=norm(jac@p+residual)/max(norm(residual),1e-300),
            actual_analytic_merit_slope=dot(residual,correct_jac@p),
            candidate_norms=candidates,accepted_cuts=[i for i,v in enumerate(candidates) if v<norm(residual)],
            best_norm_change=min(candidates)-norm(residual))
    except np.linalg.LinAlgError as error:result['error']=str(error)
    return result

def descent(folder,out):
    m=Model(*(fixtures.read(folder/(key+'.json')) for key in ('mesh','ground','case')));s=m.system
    trace=events(folder/'trace.jsonl');birth=[v for v in trace if v['kind']=='linear' and v['stage']=='birth']
    records=[]
    for label,e in (('first',birth[0]),('last',birth[-1])):
        joint={k:np.asarray(v) for k,v in e['joint'].items()};counts=np.asarray(e['contact_partitions']);x=np.asarray(e['free'])
        t=np.asarray(e['soil_t']);cr=s.c@t;fixed=np.asarray(e['fixed']);bulk=np.asarray(e['bulk_matrix'])
        constant=np.asarray(e['constant']);target=np.asarray(e['target']);r=-np.asarray(e['rhs']);a=np.asarray(e['matrix'])
        q,da,branches=analytic(s,joint,target,counts);a_ad=bulk+assemble(cr,da,s.weight)
        residual=lambda z:bulk@z+constant+np.einsum('qai,qa,q->i',cr,projected(s,joint,np.einsum('qag,g->qa',s.c,fixed+t@z),counts),s.weight)
        p=np.linalg.solve(a,-r);defect=dot(r,(a_ad-a)@p)
        no_reference=copy.deepcopy(joint)
        no_reference['jump']-=joint['reference'];no_reference['rest']-=joint['reference']
        q_without_reference=projected(s,no_reference,target,counts)
        contact_without_reference=np.einsum('qai,qa,q->i',cr,q_without_reference,s.weight)
        reference_term=np.einsum('qai,qa,q->i',cr,q-q_without_reference,s.weight)
        residual_terms=dict(bulk_and_seed=bulk@x+t.T@(m.k@fixed+m.seed),
            physical_external=-t.T@m.base,
            birth_holding_release=-(1-e['progress'])*(t.T@np.asarray(e['holding'])),
            contact_without_birth_reference=contact_without_reference,birth_reference_shift=reference_term,
            initial_force_projection=np.zeros_like(r))
        components={'bulk':bulk}
        for branch in ('open','stick','slip'):
            mask=np.asarray([b['branch']==branch for b in branches]);components[branch]=assemble(cr,da,s.weight,mask)
        variants=[]
        for use_ad,remove_projection,remove_generation,half_step in itertools.product((False,True),repeat=4):
            j=copy.deepcopy(joint);c=constant.copy()
            if remove_projection:
                # B has exactly zero projected initial traction. No fitted force.
                initial=j['initial'].copy();j['traction']-=initial
                qn=np.einsum('qa,qa->q',initial,s.normal)
                j['rest']+=qn[:,None]*s.normal/s.kn[:,None]+(initial-qn[:,None]*s.normal)/s.kt[:,None]
            if remove_generation:
                j['jump']-=j['reference'];j['rest']-=j['reference']
                c+=t.T@((1-e['progress'])*np.asarray(e['holding']))
            elif half_step:c-=t.T@(e['progress']/2*np.asarray(e['holding']))
            def evaluate(z):
                ztarget=np.einsum('qag,g->qa',s.c,fixed+t@z)
                return bulk@z+c+np.einsum('qai,qa,q->i',cr,projected(s,j,ztarget,counts),s.weight)
            rv=evaluate(x);_,dv,_=analytic(s,j,target,counts)
            correct=bulk+assemble(cr,dv,s.weight)
            tangent=correct if use_ad else bulk+assemble(cr,fd(s,j,target,counts,2e-9),s.weight)
            variants.append(dict(use_analytic=use_ad,remove_initial_projection=remove_projection,
                remove_birth_reference_and_holding=remove_generation,halve_release_step=half_step,
                residual_change=norm(rv-r),**direction_metrics(rv,tangent,correct,evaluate,x)))
        # These are read-only evaluations below the existing 18-candidate grid.
        # None enters line search or a committed state.
        extra=[dict(cut=cut,alpha=2.**(-cut),residual_norm=norm(residual(x+p*2.**(-cut))),
            decrease=norm(residual(x+p*2.**(-cut)))-norm(r),adopted=False) for cut in (18,20,24,28,32)]
        records.append(dict(label=label,iteration=e['iteration'],progress=e['progress'],
            source_jacobian=assess(a,-r),analytic_jacobian=assess(a_ad,-r),
            residual_reconstruction_error=norm(residual(x)-r),analytic_traction_maximum_error=float(np.max(abs(q-np.asarray(e['traction'])))),
            source_direction=p.tolist(),source_actual_merit_slope=dot(r,a_ad@p),
            source_direction_norm=norm(p),outside_candidate_grid_observations=extra,
            analytic_branch_boundary_flags=branches,
            residual_component_sum_error=norm(sum(residual_terms.values(),np.zeros_like(r))-r),
            residual_components={key:dict(norm=norm(v),merit_slope_contribution=dot(v,a_ad@p)) for key,v in residual_terms.items()},
            ideal_newton_merit_slope=-dot(r,r),tangent_defect_merit_slope=defect,
            slope_identity_error=abs(dot(r,a_ad@p)-(-dot(r,r)+defect+dot(r,a@p+r))),
            branch_slope_components={key:dot(r,k@p) for key,k in components.items()},
            initial_projected_force_norm=norm(joint['initial']),birth_reference_norm=norm(joint['reference']),
            holding_force_norm=norm(e['holding']),factorial=variants,
            counterfactuals_are_residual_evaluations_only=True,physical_instability_proved=False))
    write(out/'descent-decomposition.json',dict(passed=False,diagnostic_only=True,records=records,
        strict_decrease_and_18_candidates_unchanged=True,new_equilibrium_from_counterfactuals=False))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--formal',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--contact-probe',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False)
    modes(a.formal/'stacked__native__Pnative__Rnative__g8',a.output,a.contact_probe.resolve())
    descent(a.formal/'stacked__independent__Pindependent__Rindependent__g8',a.output)

if __name__=='__main__':main()
