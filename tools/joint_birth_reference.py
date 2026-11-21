"""Independent contact-BIRTH evidence. No existing verify allowlist is changed.

Point integration uses a scalar normal/two-vector Coulomb return and a fresh
birth reference. The small 3D reference assembles its own elastic stiffness,
exact ties, initial multipliers, pore forces and physical surface loads from
inputs. It refuses branches outside its declared scope.
"""
import copy
import math
import numpy as np
import fem_reference as fem
from joint_reference import System
from joint_water_reference import Water
from verify_high_order import physical


def point_path(data):
    m = data['material']; kn = m['normal_stiffness_kpa_m']; kt = m['tangential_stiffness_kpa_m']; mu = m['friction_coefficient']
    n = np.asarray(data['normal'], dtype=float); n /= np.linalg.norm(n)
    birth = np.asarray(data['reference'], dtype=float); q = np.asarray(data['traction'], dtype=float)
    normal = n@q; shear = q-normal*n
    assert normal <= 0 and np.linalg.norm(shear) <= -mu*normal
    rest = birth-normal*n/kn-shear/kt
    def energy(t):
        p=n@t; return p*p/(2*kn)+np.linalg.norm(t-p*n)**2/(2*kt)
    initial_energy=energy(q); jump=birth.copy(); slip=np.zeros(3); accumulated=plastic=numerical=work=0.
    def snapshot():
        return dict(jump_m=jump.copy(),traction_kpa=q.copy(),plastic_slip_m=slip.copy(),accumulated_plastic_slip_m=accumulated,
                    plastic_dissipation_kpa_m=plastic,numerical_dissipation_kpa_m=numerical,endpoint_work_kpa_m=work,initial_energy=initial_energy)
    results=[snapshot()]
    for increment in data['increments']:
        target=birth+increment; start=jump.copy(); count=max(1,math.ceil(np.linalg.norm(target-start)/data['maximum_increment_m']))
        for k in range(1,count+1):
            end=target if k==count else start+(target-start)*k/count
            g0=n@(jump-rest); g1=n@(end-rest)
            path=[end]
            if (g0<0<g1) or (g1<0<g0):
                cross=jump+(end-jump)*g0/(g0-g1); cross-=n*(n@(cross-rest)); path=[cross,end]
            for segment_index,next_jump in enumerate(path):
                g=n@(next_jump-rest); delta=next_jump-jump
                if g>0 or (g0>0 and segment_index==0 and (len(path)==2 or g1==0)):
                    q=np.zeros(3);jump=next_jump.copy();continue
                old_q=q.copy(); old_energy=energy(q)
                trial=q-(n@q)*n+kt*(delta-(n@delta)*n)
                cap=-mu*kn*min(g,0.); length=np.linalg.norm(trial)
                gamma=max(length-cap,0.)/kt
                direction=trial/length if length else np.zeros(3)
                plastic_increment=gamma*direction
                q=kn*min(g,0.)*n+(min(length,cap)*direction)
                loss=(q-(n@q)*n)@plastic_increment
                endpoint=q@delta
                num=endpoint-(energy(q)-old_energy)-loss
                assert num>=-1e-10*(1+abs(endpoint)+old_energy)
                slip+=plastic_increment;accumulated+=gamma;plastic+=loss;numerical+=max(num,0);work+=endpoint
                jump=next_jump.copy()
        results.append(snapshot())
    return results


def basis(mesh,case,tied,contact):
    """Global constraint equations and SVD nullspace, independent of Rust T."""
    xyz=np.asarray(mesh['nodes']); rows=[]; values=[]; ndof=xyz.size
    def add(node,d,value):
        row=np.zeros(ndof);row[3*node:3*node+3]=d;rows.append(row);values.append(value)
    if case['settings']['kinematics']=='side_rollers':
        for node,p in enumerate(xyz):
            for a in range(3):
                fixed=abs(p[2]-xyz[:,2].min())<1e-8 if a==2 else min(abs(p[a]-xyz[:,a].min()),abs(p[a]-xyz[:,a].max()))<1e-8
                if fixed:add(node,np.eye(3)[a],0.)
    else:assert case['settings']['kinematics']=='explicit'
    for c in case.get('boundary',{}).get('constraints',[]):
        for node in mesh['node_sets'][c['node_set']]:
            add(node,c['direction'],c['value_m']+np.asarray(c.get('gradient',[0,0,0]))@(xyz[node]-c.get('reference_point_m',[0,0,0])))
    if tied:
        pairs={(int(a),int(b)) for p in contact.points for a,b in zip(p['nodes_a'],p['nodes_b'])}
        for a,b in sorted(pairs):
            for d in range(3):
                row=np.zeros(ndof);row[3*a+d]=-1;row[3*b+d]=1;rows.append(row);values.append(0.)
    a=np.asarray(rows).reshape(-1,ndof);b=np.asarray(values)
    _,sv,vh=np.linalg.svd(a,full_matrices=True);rank=int((sv>1e-10*max(sv.max(initial=0),1.)).sum())
    offset=np.linalg.lstsq(a,b,rcond=1e-12)[0] if rows else np.zeros(ndof)
    assert np.max(np.abs(a@offset-b),initial=0)<1e-10
    return vh[rank:].T,offset


def elastic_reference(mesh,ground,case,order=5):
    assert not case.get('beam_system') and not case.get('anchor_system') and not case.get('initially_inactive_bodies')
    assert all(s['kind'] in {'pressure','boundary','joint_birth','hydrostatic_water','steady_seepage'} for s in case['stages'])
    assert len(case['joint_system']['joints'])==1
    definition=case['joint_system']['joints'][0]; assert definition['birth_policy'] in ('transfer_force','stress_free') and not definition.get('crack_tip')
    wet=definition.get('water') is not None
    water=Water(mesh,ground,case) if wet else None
    contact=water.contact if wet else System(mesh,case['joint_system'])
    xyz=np.asarray(mesh['nodes']);ndof=xyz.size;assert ndof<2000
    active=np.ones(len(mesh['tetrahedra']),dtype=bool);regions={r['id']:r['body'] for r in mesh['regions']}
    materials={m['id']:m['mohr_coulomb'] for m in ground['materials']}
    k=np.zeros((ndof,ndof));operators=[]
    seed=ground['initialization']['geostatic']['seed'];assert seed['kind']=='body_states'
    for ei,mid in enumerate(mesh['material_ids']):
        m=materials[mid];d=fem.elastic_tensor(m['young_kpa'],m['poisson'])
        q=fem.volume(mesh,ei,independent=True,independent_order=order);dofs=(3*q['nodes'][:,None]+np.arange(3)).ravel()
        k[np.ix_(dofs,dofs)]+=np.einsum('qai,ab,qbj,q->ij',q['b'],d,q['b'],q['weight'],optimize=True)
        s=seed['states'][regions[mesh['region_ids'][ei]]];assert s['equivalent_plastic_strain']==0
        operators.append((ei,q,dofs,d,np.asarray(s['effective_stress_compression_positive_kpa'])))
    constitutive=np.asarray([kt*np.eye(3)+(kn-kt)*np.outer(n,n) for kt,kn,n in zip(contact.kt,contact.kn,contact.normal)])
    cj=np.einsum('qai,qab,qbj,q->ij',contact.c,constitutive,contact.c,contact.weight,optimize=True)
    free,offset=basis(mesh,case,False,contact);tied,tied_offset=basis(mesh,case,True,contact)
    born=False;birth_jump=None;initial_q=None;reference_force=np.zeros(ndof)
    factor=0.;pressure=ground['initialization']['surface_surcharge_kpa'];level=ground['initialization']['water_level_m'];field=None;frames=[]
    u=np.zeros(ndof)
    for stage in [None]+case['stages']:
        if stage:
            kind=stage['kind']
            if kind=='pressure':pressure=stage['target_pressure_kpa']
            elif kind=='boundary':factor=stage['target_factor']
            elif kind=='hydrostatic_water':level=stage['target_level_m'];field=None
            elif kind=='steady_seepage':
                assert wet;level=stage['exterior_water_level_m'];field=water.steady(stage['flow'],active,mesh['material_ids'],level)['field']
            elif kind=='joint_birth':
                assert not born
                birth_jump=np.einsum('qag,g->qa',contact.c,u)
                if definition['birth_policy']=='stress_free':initial_q=np.zeros_like(birth_jump)
                else:
                    # Work-conjugate nodal multiplier recovery; square-root
                    # surface measure gives the minimum L2 GP field.
                    a=np.einsum('qag,q->gqa',contact.c,np.sqrt(contact.weight)).reshape(ndof,-1)
                    rhs=last_load-last_seed-k@u
                    solution=np.linalg.lstsq(free.T@a,free.T@rhs,rcond=1e-12)[0]
                    assert np.max(np.abs(free.T@(a@solution-rhs)),initial=0)<1e-10
                    initial_q=solution.reshape(-1,3)/np.sqrt(contact.weight)[:,None]
                reference_force=np.einsum('qag,qa,q->g',contact.c,initial_q,contact.weight)-np.einsum('qag,qab,qb,q->g',contact.c,constitutive,birth_jump,contact.weight)
                born=True
        seed_force=np.zeros(ndof)
        for ei,q,dofs,d,s0 in operators:
            if not wet or not water.soil[ei]:pore=np.zeros(len(q['weight']))
            elif field is None:
                assert level is None or level>=xyz[:,2].max(),'partial volume wetting outside elastic reference'
                pore=np.zeros(len(q['weight'])) if level is None else water.gamma*np.maximum(level-q['position'][:,2],0.)
            else:
                h=field['relative_head'][q['nodes']];pore=water.gamma*(field['datum'][q['nodes'][0]]-q['position'][:,2]+h[0]+q['shape']@(h-h[0]));assert pore.min()>=0
            total=np.broadcast_to(s0,(len(pore),6)).copy();total[:,:3]+=pore[:,None]
            seed_force[dofs]-=np.einsum('qai,qa,q->i',q['b'],total,q['weight'])
        fg=copy.deepcopy(ground);fg['initialization']['water_level_m']=level
        load,_=physical(mesh,fg,active,pressure,case.get('pressure_on_original_surface_only',False),joint_system=contact)
        if case.get('boundary'):load+=factor*fem.conditions_load(mesh,case['boundary'],active)
        interface=water.hydrostatic(active,level) if wet and field is None else water.from_heads(active,field,level) if wet else None
        if wet:load+=interface['load']
        load=load.ravel();t=free if born else tied;o=factor*(offset if born else tied_offset);matrix=k+cj if born else k
        rhs=t.T@(load-seed_force-reference_force-matrix@o)
        reduced=t.T@matrix@t;solution=np.linalg.solve(reduced,rhs);u=o+t@solution
        linear=np.linalg.norm(reduced@solution-rhs)/max(np.linalg.norm(rhs),1.);assert linear<1e-10
        traction=np.einsum('qab,qb->qa',constitutive,np.einsum('qag,g->qa',contact.c,u)-birth_jump)+initial_q if born else np.zeros((len(contact.points),3))
        if born:
            pn=-np.einsum('qa,qa->q',contact.normal,traction);shear=traction+pn[:,None]*contact.normal
            assert pn.min()>1e-6 and np.all(np.linalg.norm(shear,axis=1)<contact.mu*pn-1e-6),'nonlinear contact branch needs separate integration evidence'
        frames.append(dict(displacement=u.reshape(-1,3),traction=traction,load=load,seed_force=seed_force,bulk=k@u,linear=linear,interface=interface,field=field,born=born))
        last_load=load;last_seed=seed_force
    return frames


def audit_elastic(mesh,ground,run):
    assert run['complete'],'incomplete, cancelled or timed-out run cannot pass'
    expected=elastic_reference(mesh,ground,run['case']);assert len(expected)==len(run['frames'])
    head_errors={}
    for audit in run.get('flow_audits',[]):
        reference=expected[audit['stage_index']+1]['field'];native=audit['joint_flow']['original_field']
        defined=reference['defined'];error=float(np.max(np.abs(np.asarray(native['head_m'])[defined]-reference['head'][defined]),initial=0.))
        assert error<1e-8
        head_errors[audit['stage_index']+1]=error
    rows=[]
    for frame_index,(reference,frame) in enumerate(zip(expected,run['frames'])):
        native=np.asarray(frame['state']['displacement_m']);delta=np.max(np.abs(reference['displacement']-native));relative=np.linalg.norm(reference['displacement']-native)/max(np.linalg.norm(reference['displacement']),1e-30)
        assert relative<1e-7,(frame['name'],relative)
        load_error=np.max(np.abs(reference['load']-frame['state']['applied_nodal_load_kn']));assert load_error<1e-10
        assert frame['relative_force_residual']<1e-8 and frame['maximum_constraint_error_m']<1e-10
        q=np.asarray([p['traction_kpa'] for p in frame['state']['joints']['points']]);q_error=np.max(np.abs(q-reference['traction']));assert q_error<1e-7
        if reference['born']:
            for p in frame['state']['joints']['points']:assert p['plastic_slip_m']==[0.,0.,0.] and p['plastic_dissipation_kpa_m']==0.
        if reference['interface'] is not None:
            water=reference['interface'];native_water=frame['state']['joint_water']
            assert np.max(np.abs(water['load'].ravel()-native_water['nodal_load_kn']))<1e-10
            assert np.max(np.abs(water['pressure']-native_water['physical_pressure_kpa']))<1e-8
        rows.append(dict(name=frame['name'],displacement_absolute_m=float(delta),displacement_relative=float(relative),nodal_load_error_kn=float(load_error),gp_traction_error_kpa=float(q_error),independent_linear_residual=float(reference['linear']),water_head_error_m=head_errors.get(frame_index,0.)))
    return dict(passed=True,frames=rows)


def audit_hydrostatic_reassembly(mesh,ground,run):
    """Independent force reassembly, including nonlinear/open native histories.

    Material stresses are audit data here. This is distinct from the independent
    elastic solve above and never labels a plastic 3D solve as independent.
    """
    from verify_high_order import prepare,internal
    assert run['complete'] and run['schema']=='soil-fem-stage-result/15'
    case=run['case']
    assert not case.get('beam_system') and not case.get('anchor_system') and not case.get('initially_inactive_bodies')
    assert len(case['joint_system']['joints'])==1 and not case['joint_system']['joints'][0].get('crack_tip')
    assert all(s['kind'] in {'pressure','boundary','hydrostatic_water','joint_birth','joint_death'} for s in case['stages']),'domain changes/flow need a separate reassembly scope'
    assert len(run['frames'])==len(case['stages'])+1
    assert all(all(f['state']['active_elements']) and f['state'].get('pore_pressure_kpa') is None for f in run['frames'])
    water=Water(mesh,ground,case);q,offset=prepare(mesh,run['integration_points']);rows=[]
    for frame in run['frames']:
        state=frame['state'];active=np.asarray(state['active_elements']);level=state.get('hydrostatic_water_level_m',ground['initialization']['water_level_m'])
        fg=copy.deepcopy(ground);fg['initialization']['water_level_m']=level
        pore=np.concatenate([water.gamma*np.maximum(level-v['position'][:,2],0.) if water.soil[e] else np.zeros(len(v['weight'])) for e,v in enumerate(q)])
        force,_=internal(mesh,fg,state['material_states'],active,q,offset,state.get('element_material_ids'),pore)
        interface=water.hydrostatic(active,level)
        load,_=physical(mesh,fg,active,state['surface_pressure_kpa'],case['pressure_on_original_surface_only'],state.get('element_material_ids'),water.contact)
        if case.get('boundary'):load+=state['boundary_factor']*fem.conditions_load(mesh,case['boundary'],active)
        load+=interface['load'];load=load.ravel();force=force.ravel()
        life=state['joints']['lifecycle'][0];traction=np.asarray([p['traction_kpa'] for p in state['joints']['points']])
        if life['active']:force+=np.einsum('qag,qa,q->g',water.contact.c,traction,water.contact.weight)
        t,o=basis(mesh,case,life['tied'],water.contact)
        residual=np.linalg.norm(t.T@(load-force))/max(np.linalg.norm(load),1.);assert residual<1e-8
        load_error=np.max(np.abs(load-state['applied_nodal_load_kn']));assert load_error<1e-10
        u=np.asarray(state['displacement_m']).ravel()-state['boundary_factor']*o
        constraint=np.max(np.abs(u-t@(t.T@u)));assert constraint<1e-10
        field=state['joint_water'];field_error=np.max(np.abs(interface['load'].ravel()-field['nodal_load_kn']));assert field_error<1e-10
        rows.append(dict(name=frame['name'],relative_force_residual=float(residual),nodal_load_error_kn=float(load_error),constraint_error_m=float(constraint),water_face_error_kn=float(field_error)))
    return dict(passed=True,frames=rows)
