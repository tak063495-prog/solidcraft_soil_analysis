"""Bounded independent references for extension (a).

Does not change any historical verifier's allowlist. The tip reference uses
input-only GL125 volume integration, Tri6 surface integration and a global
constraint nullspace. Saved material histories are used only in the separately
labelled nodal reassembly, never as an independent elastic solution.
"""
import copy
import json
import numpy as np
import fem_reference as fem
from joint_reference import System
from joint_water_reference import Water
from verify_high_order import physical, prepare, internal


def audit_clean_path(full,clean):
    """Compare the complete accepted Run, including excavation/geostatic audits.

    Input injection controls, the initial step used to generate the same
    accepted sequence, and their input fingerprint necessarily differ. Remove
    only those declared controls, rejected attempts and elapsed measurements.
    Every physical input and every other output field remains in the equality.
    """
    selected=full['case']['verification_reject_once_stages']
    assert clean['case']['verification_reject_once_stages']==[]
    assert full['complete'] and clean['complete']
    def without_clock(v):
        if isinstance(v,dict):return {k:without_clock(x) for k,x in v.items() if k!='elapsed_seconds'}
        if isinstance(v,list):return [without_clock(x) for x in v]
        return v
    def accepted(r):
        r=copy.deepcopy(r)
        r['case']['verification_reject_once_stages']=[]
        for s in r['case']['stages']:
            if s['name'] in selected:s.pop('initial_increment',None)
        r['provenance'].pop('case_sha256')
        r['attempts']=[a for a in r['attempts'] if a['accepted']]
        return json.dumps(without_clock(r),sort_keys=True,separators=(',',':'),allow_nan=False)
    assert accepted(full)==accepted(clean),'complete accepted Run differs'
    return dict(passed=True,complete_accepted_run_exact=True,excluded_only=['elapsed_seconds','rejected attempts','verification_reject_once_stages','selected stages initial_increment','case input fingerprint'])


def constraint_space(mesh, case, active, contact=None, tied=False):
    xyz = np.asarray(mesh['nodes']); nd = xyz.size
    used = set(int(n) for e in np.flatnonzero(active) for n in fem.cells(mesh)[e])
    rows, values = [], []
    def add(row, value=0.):
        rows.append(row); values.append(value)
    assert case['settings']['kinematics'] in ('explicit', 'side_rollers')
    for n in range(len(xyz)):
        for d in range(3):
            fixed = n not in used
            if case['settings']['kinematics'] == 'side_rollers':
                fixed |= abs(xyz[n,2]-xyz[:,2].min())<1e-8 if d==2 else min(abs(xyz[n,d]-xyz[:,d].min()),abs(xyz[n,d]-xyz[:,d].max()))<1e-8
            if fixed:
                row = np.zeros(nd); row[3*n+d] = 1.; add(row)
    for c in case.get('boundary', {}).get('constraints', []):
        for n in mesh['node_sets'][c['node_set']]:
            if n in used:
                row = np.zeros(nd); row[3*n:3*n+3] = c['direction']
                add(row,c['value_m']+np.asarray(c.get('gradient',[0.,0.,0.]))@(xyz[n]-c.get('reference_point_m',[0.,0.,0.])))
    if tied:
        assert contact is not None
        pairs = {(a,b) for p in contact.points if active[p['owner_a']] and active[p['owner_b']] for a,b in zip(p['nodes_a'],p['nodes_b']) if a!=b}
        for a,b in sorted(pairs):
            for d in range(3):
                row=np.zeros(nd); row[3*a+d]-=1.; row[3*b+d]+=1.; add(row)
    matrix=np.asarray(rows).reshape(-1,nd); values=np.asarray(values)
    _,sv,vh=np.linalg.svd(matrix,full_matrices=True)
    rank=int((sv>1e-10*max(sv.max(initial=0.),1.)).sum())
    offset=np.linalg.lstsq(matrix,values,rcond=1e-12)[0]
    assert np.max(np.abs(matrix@offset-values),initial=0.)<1e-10
    return vh[rank:].T, offset


def elastic_operators(mesh, ground):
    materials={m['id']:m['mohr_coulomb'] for m in ground['materials']}
    result=[]
    for e,mid in enumerate(mesh['material_ids']):
        q=fem.volume(mesh,e,independent=True,independent_order=5)
        d=fem.elastic_tensor(materials[mid]['young_kpa'],materials[mid]['poisson'])
        ids=(3*q['nodes'][:,None]+np.arange(3)).ravel()
        k=np.einsum('qai,ab,qbj,q->ij',q['b'],d,q['b'],q['weight'],optimize=True)
        result.append((q,d,ids,k))
    return result


def tip_elastic(mesh, ground, case):
    assert ground['initialization']['water_level_m'] is None
    assert not case.get('initially_inactive_bodies') and not case.get('beam_system') and not case.get('anchor_system')
    assert len(case['joint_system']['joints'])==1
    definition=case['joint_system']['joints'][0]
    assert definition.get('crack_tip') and definition['birth_policy'] in ('transfer_force','stress_free') and not definition.get('water')
    assert all(s['kind'] in {'joint_birth','joint_death','pressure','boundary'} for s in case['stages'])
    contact=System(mesh,case['joint_system']); nd=3*len(mesh['nodes']); assert nd<=2000
    operators=elastic_operators(mesh,ground); k=np.zeros((nd,nd)); seed_force=np.zeros(nd)
    bodies={r['id']:r['body'] for r in mesh['regions']}
    seed=ground['initialization']['geostatic']['seed']; assert seed['kind']=='body_states'
    for e,(q,d,ids,ke) in enumerate(operators):
        k[np.ix_(ids,ids)]+=ke
        s=seed['states'][bodies[mesh['region_ids'][e]]]
        assert s['equivalent_plastic_strain']==s['plastic_work_kpa']==0.
        seed_force[ids]-=np.einsum('qai,a,q->i',q['b'],s['effective_stress_compression_positive_kpa'],q['weight'])
    tensor=np.asarray([kt*np.eye(3)+(kn-kt)*np.outer(n,n) for kt,kn,n in zip(contact.kt,contact.kn,contact.normal)])
    cj=np.einsum('qag,qab,qbh,q->gh',contact.c,tensor,contact.c,contact.weight,optimize=True)
    active=np.ones(len(operators),bool)
    free,offset=constraint_space(mesh,case,active)
    tied,tied_offset=constraint_space(mesh,case,active,contact,True)
    born=False; initial_q=None; birth_jump=None; reference_force=np.zeros(nd)
    factor=ground['initialization']['geostatic'].get('initial_boundary_factor',0.)
    pressure=ground['initialization']['surface_surcharge_kpa']; frames=[]; u=np.zeros(nd)
    for stage in [None]+case['stages']:
        if stage:
            kind=stage['kind']
            if kind=='boundary': factor=stage['target_factor']
            elif kind=='pressure': pressure=stage['target_pressure_kpa']
            elif kind=='joint_death':
                assert born; born=False; reference_force=np.zeros(nd)
            elif kind=='joint_birth':
                assert not born
                birth_jump=np.einsum('qag,g->qa',contact.c,u)
                initial_q=np.zeros_like(birth_jump)
                if definition['birth_policy']=='transfer_force':
                    # Recover physical tie forces first. Shared tip nodes never
                    # enter this unknown set, and no GP value is clipped.
                    pairs=sorted({(a,b) for p in contact.points for a,b in zip(p['nodes_a'],p['nodes_b']) if a!=b})
                    pmat=np.zeros((nd,3*len(pairs)))
                    for j,(a,b) in enumerate(pairs):
                        for axis in range(3):pmat[3*a+axis,3*j+axis]-=1.;pmat[3*b+axis,3*j+axis]+=1.
                    multiplier=np.linalg.lstsq(free.T@pmat,free.T@(last_load-seed_force-k@u),rcond=1e-12)[0]
                    target=pmat@multiplier
                    assert np.max(np.abs(free.T@(target-last_load+seed_force+k@u)),initial=0.)<1e-10
                    a=np.einsum('qag,q->gqa',contact.c,np.sqrt(contact.weight)).reshape(nd,-1)
                    projection=np.linalg.lstsq(a,target,rcond=1e-12)[0]
                    assert np.max(np.abs(a@projection-target),initial=0.)<1e-10
                    initial_q=projection.reshape(-1,3)/np.sqrt(contact.weight)[:,None]
                    pn=-np.einsum('qa,qa->q',contact.normal,initial_q)
                    assert pn.min()>=0. and np.all(np.linalg.norm(initial_q+pn[:,None]*contact.normal,axis=1)<=contact.mu*pn)
                reference_force=np.einsum('qag,qa,q->g',contact.c,initial_q,contact.weight)-cj@u
                born=True
        load,_=physical(mesh,ground,active,pressure,case['pressure_on_original_surface_only'],joint_system=contact)
        load+=factor*fem.conditions_load(mesh,case.get('boundary',{}),active)
        load=load.ravel(); t=free if born else tied; o=factor*(offset if born else tied_offset)
        matrix=k+cj if born else k
        rhs=t.T@(load-seed_force-reference_force-matrix@o); reduced=t.T@matrix@t
        v=np.linalg.solve(reduced,rhs); u=o+t@v
        linear=np.linalg.norm(reduced@v-rhs)/max(np.linalg.norm(rhs),1.)
        assert linear<1e-10
        traction=np.einsum('qab,qb->qa',tensor,np.einsum('qag,g->qa',contact.c,u)-birth_jump)+initial_q if born else np.zeros((len(contact.points),3))
        if born:
            pn=-np.einsum('qa,qa->q',contact.normal,traction)
            assert pn.min()>0. and np.all(np.linalg.norm(traction+pn[:,None]*contact.normal,axis=1)<contact.mu*pn),'nonlinear branch outside independent elastic scope'
        frames.append(dict(u=u.copy(),traction=traction,load=load,linear=float(linear),born=born))
        last_load=load
    return frames


def audit_tip_elastic(mesh,ground,run):
    assert run['complete'] and len(run['frames'])==len(run['case']['stages'])+1
    expected=tip_elastic(mesh,ground,run['case']); rows=[]
    operators=elastic_operators(mesh,ground)
    bodies={r['id']:r['body'] for r in mesh['regions']}
    seed=ground['initialization']['geostatic']['seed']['states']
    for frame,r in zip(run['frames'],expected):
        u=np.asarray(frame['state']['displacement_m']).ravel()
        error=float(np.linalg.norm(u-r['u'])/max(np.linalg.norm(r['u']),1e-30)); assert error<1e-7,(frame['name'],error)
        native_q=np.asarray([p['traction_kpa'] for p in frame['state']['joints']['points']])
        if r['born']:
            qe=float(np.linalg.norm(native_q-r['traction'])/max(np.linalg.norm(r['traction']),1.)); assert qe<1e-7
        else: qe=0. # Dormant native history is frozen, not a contact force.
        stresses=[]
        for e,(q,d,ids,_) in enumerate(operators):
            q_native=fem.volume(mesh,e)
            s0=seed[bodies[mesh['region_ids'][e]]]['effective_stress_compression_positive_kpa']
            stresses.extend(np.asarray(s0)-np.einsum('ab,qbi,i->qa',d,q_native['b'],r['u'][ids]))
        native_s=np.asarray([s['effective_stress_compression_positive_kpa'] for s in frame['state']['material_states']])
        se=float(np.linalg.norm(native_s-stresses)/max(np.linalg.norm(stresses),1.)); assert se<1e-7
        assert all(s['equivalent_plastic_strain']==s['plastic_work_kpa']==0. for s in frame['state']['material_states'])
        le=float(np.max(np.abs(r['load']-frame['state']['applied_nodal_load_kn']))); assert le<1e-10
        rows.append(dict(name=frame['name'],displacement_relative=error,volume_gp_stress_relative=se,contact_gp_traction_relative=qe,nodal_load_error_kn=le,independent_linear_residual=r['linear']))
    return dict(passed=True,reference='input-only GL125 elasticity, Tri6, two independent SVD projections; all volume GP stresses',frames=rows)


def audit_reassembly(mesh,ground,run):
    assert run['complete'];case=run['case']
    assert not case.get('beam_system') and not case.get('anchor_system')
    assert len(case['joint_system']['joints'])==1
    assert all(s['kind'] in {'joint_birth','joint_death','pressure','boundary','birth','excavate'} for s in case['stages'])
    definition=case['joint_system']['joints'][0];wet=definition.get('water') is not None
    if definition.get('crack_tip'):
        assert not wet and not case.get('initially_inactive_bodies')
        assert all(s['kind'] not in {'birth','excavate'} for s in case['stages'])
        assert all(all(f['state']['active_elements']) for f in run['frames'])
    elif any(s['kind']=='birth' for s in case['stages']):
        assert wet and definition['birth_policy']=='stress_free' and definition['water']['kind']=='continuous'
        assert not case['joint_system'].get('structure_bodies') and len(case['initially_inactive_bodies'])==1
    water=Water(mesh,ground,case) if wet else None
    contact=water.contact if wet else System(mesh,case['joint_system'])
    q,offset=prepare(mesh,run['integration_points']);rows=[]
    for frame in run['frames']:
        s=frame['state']; active=np.asarray(s['active_elements']); assert s.get('pore_pressure_kpa') is None
        level=s.get('hydrostatic_water_level_m',ground['initialization']['water_level_m'])
        fg=copy.deepcopy(ground);fg['initialization']['water_level_m']=level
        pore=np.concatenate([ground['initialization']['water_unit_weight_kn_m3']*np.maximum(level-op['position'][:,2],0.) if wet and water.soil[e] else np.zeros(len(op['weight'])) for e,op in enumerate(q)])
        force,_=internal(mesh,fg,s['material_states'],active,q,offset,s.get('element_material_ids'),pore)
        life=s['joints']['lifecycle'][0]
        if life['active']:
            traction=np.asarray([p['traction_kpa'] for p in s['joints']['points']])
            live=np.asarray([active[p['owner_a']] and active[p['owner_b']] for p in contact.points])
            force=force.ravel()+np.einsum('qag,qa,q->g',contact.c,traction,contact.weight*live)
        else: force=force.ravel()
        load,_=physical(mesh,fg,active,s['surface_pressure_kpa'],case['pressure_on_original_surface_only'],s.get('element_material_ids'),contact)
        load+=s['boundary_factor']*fem.conditions_load(mesh,case.get('boundary',{}),active)
        water_error=0.
        if wet:
            field=water.hydrostatic(active,level);load+=field['load']
            water_error=float(np.max(np.abs(field['load'].ravel()-s['joint_water']['nodal_load_kn'])));assert water_error<1e-10
            assert np.max(np.abs(field['pressure']-s['joint_water']['physical_pressure_kpa']))<1e-8
        load=load.ravel();t,o=constraint_space(mesh,case,active,contact,life['tied'])
        residual=float(np.linalg.norm(t.T@(force-load))/max(np.linalg.norm(load),1.));assert residual<1e-10,(frame['name'],residual)
        le=float(np.max(np.abs(load-s['applied_nodal_load_kn'])));assert le<1e-10
        u=np.asarray(s['displacement_m']).ravel()-s['boundary_factor']*o
        used=sorted({int(n) for e in np.flatnonzero(active) for n in fem.cells(mesh)[e]})
        dofs=(3*np.asarray(used)[:,None]+np.arange(3)).ravel()
        ce=float(np.max(np.abs(u-t@(t.T@u))[dofs]));assert ce<1e-10
        assert frame['relative_force_residual']<1e-8 and frame['maximum_constraint_error_m']<1e-10
        for p in s['joints']['points']:
            if p.get('birth'):
                n=np.asarray(p['normal']);traction=np.asarray(p['traction_kpa']);normal=n@traction
                elastic=normal**2/(2*definition['material']['normal_stiffness_kpa_m'])+np.linalg.norm(traction-normal*n)**2/(2*definition['material']['tangential_stiffness_kpa_m'])
                balance=p['birth']['initial_energy_kpa_m']+p['endpoint_work_kpa_m']-elastic-p['plastic_dissipation_kpa_m']-p['numerical_dissipation_kpa_m']
                assert abs(balance)<1e-10
        rows.append(dict(name=frame['name'],reassembled_relative_residual=residual,nodal_load_error_kn=le,constraint_error_m=ce,water_load_error_kn=water_error))
    linear=[v['true_relative_residual'] for a in run['geostatic']['attempts']+run['attempts'] for v in a.get('linear_solves',[])]
    assert linear and max(linear)<1e-10
    return dict(passed=True,reference='independent all-node reassembly of saved stress histories; separate from input-only elastic solution',maximum_native_linear_residual=max(linear),frames=rows)


def audit_wet_birth_elastic(mesh,ground,run):
    """Input-only placed-soil elasticity, including the new material origin.

    Fully submerged hydrostatic MC soil; one inactive body, one mode-B joint,
    and birth/contact/pressure only. No saved stress or displacement generates
    the reference solution or its birth placement.
    """
    assert run['complete'];case=run['case'];stages=case['stages']
    assert [s['kind'] for s in stages]==['birth','joint_birth','pressure']
    assert not case.get('beam_system') and not case.get('anchor_system')
    assert len(case['initially_inactive_bodies'])==1 and len(case['joint_system']['joints'])==1
    definition=case['joint_system']['joints'][0]
    assert definition['birth_policy']=='stress_free' and not definition.get('crack_tip')
    assert definition['water']['kind']=='continuous' and not case['joint_system'].get('structure_bodies')
    birth=stages[0];assert birth['bodies']==case['initially_inactive_bodies'] and birth['placement']=='cad_new_nodes'
    level=ground['initialization']['water_level_m'];assert birth['hydraulic_initialization']==dict(kind='hydrostatic',level_m=level)
    assert level>=np.asarray(mesh['nodes'])[:,2].max()
    water=Water(mesh,ground,case);assert water.soil.all()
    contact=water.contact;operators=elastic_operators(mesh,ground);nd=3*len(mesh['nodes'])
    bodies={r['id']:r['body'] for r in mesh['regions']}
    names=[bodies[r] for r in mesh['region_ids']]
    active=np.asarray([b not in birth['bodies'] for b in names])
    seed=ground['initialization']['geostatic']['seed']['states']
    origins=np.zeros((len(operators),nd));initial_stress=[np.asarray(seed[b]['effective_stress_compression_positive_kpa']) for b in names]
    factor=ground['initialization']['geostatic'].get('initial_boundary_factor',0.)
    assert factor==0.
    tensor=np.asarray([kt*np.eye(3)+(kn-kt)*np.outer(n,n) for kt,kn,n in zip(contact.kt,contact.kn,contact.normal)])
    cj=np.einsum('qag,qab,qbh,q->gh',contact.c,tensor,contact.c,contact.weight,optimize=True)
    u=np.zeros(nd);born=False;ref_force=np.zeros(nd);pressure=ground['initialization']['surface_surcharge_kpa'];rows=[]
    for index,stage in enumerate([None]+stages):
        if stage:
            if stage['kind']=='birth':
                existing={int(n) for e in np.flatnonzero(active) for n in fem.cells(mesh)[e]}
                placed=u.reshape(-1,3).copy()
                for n in set(range(len(mesh['nodes'])))-existing:placed[n]=0.
                for p in contact.points:
                    for a,b in zip(p['nodes_a'],p['nodes_b']):
                        if a in existing and b not in existing:placed[b]=placed[a]
                        elif b in existing and a not in existing:placed[a]=placed[b]
                for e,b in enumerate(names):
                    if not active[e]:
                        origins[e]=placed.ravel();state=birth['states'][b]
                        assert state['equivalent_plastic_strain']==state['plastic_work_kpa']==0. and state['plastic_strain_compression_positive']==[0.]*6
                        initial_stress[e]=np.asarray(state['effective_stress_compression_positive_kpa'])
                active[:]=True
                native_prepared=run['transition_audits'][0]['prepared_state']
                assert np.max(np.abs(placed-np.asarray(native_prepared['displacement_m'])))<1e-10
                for e,(q,d,ids,ke) in enumerate(operators):
                    if names[e] in birth['bodies']:
                        native=[s for p,s in zip(run['integration_points'],native_prepared['material_states']) if p['element_index']==e]
                        assert all(np.array_equal(s['effective_stress_compression_positive_kpa'],initial_stress[e]) for s in native)
            elif stage['kind']=='joint_birth':
                born=True;ref_force=-cj@u;birth_jump=np.einsum('qag,g->qa',contact.c,u)
            else:pressure=stage['target_pressure_kpa']
        k=np.zeros((nd,nd));seed_force=np.zeros(nd)
        for e,(q,d,ids,ke) in enumerate(operators):
            if active[e]:
                k[np.ix_(ids,ids)]+=ke
                total=np.broadcast_to(initial_stress[e],(len(q['weight']),6)).copy()
                total[:,:3]+=water.gamma*(level-q['position'][:,2,None])
                seed_force[ids]-=np.einsum('qai,qa,q->i',q['b'],total,q['weight'])+ke@origins[e,ids]
        t,o=constraint_space(mesh,case,active,contact,not born)
        load,_=physical(mesh,ground,active,pressure,case['pressure_on_original_surface_only'],joint_system=contact)
        field=water.hydrostatic(active,level);load=(load+field['load']).ravel()
        matrix=k+cj if born else k;rhs=t.T@(load-seed_force-ref_force);reduced=t.T@matrix@t
        v=np.linalg.solve(reduced,rhs);u=t@v
        linear=float(np.linalg.norm(reduced@v-rhs)/max(np.linalg.norm(rhs),1.));assert linear<1e-10
        frame=run['frames'][index];s=frame['state'];used=sorted({int(n) for e in np.flatnonzero(active) for n in fem.cells(mesh)[e]});dofs=(3*np.asarray(used)[:,None]+np.arange(3)).ravel()
        native_u=np.asarray(s['displacement_m']).ravel()
        ue=float(np.linalg.norm((u-native_u)[dofs])/max(np.linalg.norm(u[dofs]),1e-30));assert ue<1e-7,(frame['name'],ue)
        stress_error=0.;offset=0
        for e,(q,d,ids,ke) in enumerate(operators):
            op=fem.volume(mesh,e);count=len(op['weight'])
            if active[e]:
                expected=initial_stress[e]-np.einsum('ab,qbi,i->qa',d,op['b'],u[ids]-origins[e,ids])
                actual=np.asarray([s['effective_stress_compression_positive_kpa'] for s in s['material_states'][offset:offset+count]])
                se=float(np.linalg.norm(expected-actual)/max(np.linalg.norm(expected),1.));assert se<1e-7
                stress_error=max(stress_error,se)
            offset+=count
        le=float(np.max(np.abs(load-s['applied_nodal_load_kn'])));assert le<1e-10
        if born:
            traction=np.einsum('qab,qb->qa',tensor,np.einsum('qag,g->qa',contact.c,u)-birth_jump)
            pn=-np.einsum('qa,qa->q',contact.normal,traction)
            assert pn.min()>0. and np.all(np.linalg.norm(traction+pn[:,None]*contact.normal,axis=1)<contact.mu*pn)
            qe=float(np.linalg.norm(traction-np.asarray([p['traction_kpa'] for p in s['joints']['points']]))/max(np.linalg.norm(traction),1.));assert qe<1e-7
        else:qe=0.
        rows.append(dict(name=frame['name'],displacement_relative=ue,volume_gp_stress_relative=stress_error,contact_gp_traction_relative=qe,nodal_load_error_kn=le,independent_linear_residual=linear))
    return dict(passed=True,reference='input-only GL125 placed-soil elasticity with independent initial pore/seed forces and placement origin',frames=rows)


def audit_wet_birth_holding(mesh,ground,run):
    """Birth at lambda=0: independently assemble the NEW physical water RHS.

    prepared_state is the constitutive replacement snapshot, whose stored
    applied load is still the preceding stage's load. Rebuild the new domain
    loads rather than using that snapshot's old cached surface-water field.
    The holding vector must balance every physical component, including support
    components, before the first mechanical increment.
    """
    assert run['complete'];case=run['case'];water=Water(mesh,ground,case)
    assert water.soil.all() and case['joint_system']['joints'][0]['birth_policy']=='stress_free'
    q,offset=prepare(mesh,run['integration_points']);rows=[]
    for audit in run['transition_audits']:
        stage=next(s for s in case['stages'] if s['name']==audit['stage'])
        if stage['kind']!='birth':continue
        state=audit['prepared_state'];active=np.asarray(state['active_elements'])
        level=stage['hydraulic_initialization']['level_m']
        assert stage['hydraulic_initialization']['kind']=='hydrostatic' and state.get('pore_pressure_kpa') is None
        assert state.get('hydrostatic_water_level_m',ground['initialization']['water_level_m'])==level
        pore=np.concatenate([water.gamma*(level-op['position'][:,2]) for op in q]);assert pore.min()>=0.
        force,_=internal(mesh,ground,state['material_states'],active,q,offset,state.get('element_material_ids'),pore)
        assert not state['joints']['lifecycle'][0]['active']
        load,_=physical(mesh,ground,active,state['surface_pressure_kpa'],case['pressure_on_original_surface_only'],state.get('element_material_ids'),water.contact)
        load+=state['boundary_factor']*fem.conditions_load(mesh,case.get('boundary',{}),active)
        field=water.hydrostatic(active,level);load+=field['load']
        expected=(force-load).ravel();holding=np.asarray(audit['initial_holding_load_kn'])
        error=float(np.max(np.abs(expected-holding)));assert error<1e-10,(audit['stage'],error)
        t,_=constraint_space(mesh,case,active,water.contact,True)
        residual=float(np.linalg.norm(t.T@(load.ravel()+holding-force.ravel()))/max(np.linalg.norm(load),1.));assert residual<1e-10
        assert audit['zero_release_relative_force_residual']<1e-10
        selected=np.concatenate([pore[offset[e]:offset[e+1]] for e in audit['selected_elements']])
        rows.append(dict(stage=audit['stage'],all_nodal_holding_load_error_kn=error,zero_release_relative_residual=residual,initial_pore_pressure_min_kpa=float(selected.min()),initial_pore_pressure_max_kpa=float(selected.max()),fresh_joint_water_load_norm_kn=float(np.linalg.norm(field['load']))))
    assert rows
    return dict(passed=True,frames=rows)
