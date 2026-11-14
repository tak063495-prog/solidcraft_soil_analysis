"""All-component comparisons with the unchanged A-02 difference convention."""
import itertools
import json
import numpy as np
CATEGORIES={'material_return','contact_status','contact_open','joint_active','joint_tied','generation',
            'archived_generation_count','active_elements','exterior_faces','wet_exterior_faces'}


def row(stage,progress):
    return dict(stage=stage,progress=progress,hydraulic_accounting=dict(
        model='prescribed_hydrostatic',flux=None,storage=None,birth_water_transfer=None,balance_acceptance=None))


def native_fields(m,state,d,run,prepared=False):
    from joint_ext_reference import constraint_space
    materials=state['material_states'];points=state['joints']['points'];s=m.system
    u=np.asarray(state['displacement_m']).ravel();active=np.asarray(state['active_elements'])
    origins=np.zeros((len(m.ops),m.nd))
    birth=m.case['stages'][0]
    if active.all():
        audit=next(a for a in run['transition_audits'] if a['stage']==birth['name'])
        placed=np.asarray(audit['prepared_state']['displacement_m']).ravel()
        for e,b in enumerate(m.names):
            if b=='UpperClay':origins[e]=placed
    strain=np.concatenate([np.einsum('qai,i->qa',op['b'],u[ids]-origins[e,ids]) if active[e] else np.zeros((len(op['weight']),6))
                           for e,((_,_,ids,_),op) in enumerate(zip(m.ops,m.sample))])
    v=lambda k:np.asarray([p[k] for p in points]);b=lambda k:np.asarray([p[k] for p in materials])
    traction=v('traction_kpa');qn=np.einsum('qa,qa->q',traction,s.normal);tau=traction-qn[:,None]*s.normal
    reference=np.asarray([p.get('birth',{}).get('jump_m',[0.,0.,0.]) for p in points])
    rest=np.asarray([p.get('birth',{}).get('elastic_offset_m',[0.,0.,0.]) for p in points])
    gap=np.einsum('qa,qa->q',v('jump_m')-rest,s.normal)
    holding=state.get('geostatic_release_load_kn') if state['stage_index']==0 else state.get('transition_release_load_kn')
    if prepared:holding=next(a['initial_holding_load_kn'] for a in run['transition_audits'] if a['stage']==birth['name'])
    if not holding:holding=np.zeros(m.nd)
    load=np.asarray(state['applied_nodal_load_kn'])
    if prepared:
        load=sum(np.asarray(d[k]) for k in ('gravity_load','base_load','surcharge_load','exterior_water_load','joint_water_load'))+holding
    force=np.asarray(d['internal']);life=state['joints']['lifecycle']
    assert d['face_ids']==[list(f) for f,_ in m.faces]
    assert len(materials)==sum(len(q['weight']) for q in m.sample) and len(points)==len(s.points)
    pore=np.asarray(d['pore_pressure'])
    result=dict(u=u.reshape(-1,3),stress=b('effective_stress_compression_positive_kpa'),strain=strain,void=b('void_ratio'),
        total_stress=b('effective_stress_compression_positive_kpa')+pore[:,None]*np.array([1,1,1,0,0,0]),pore_pressure=pore,head=d['head'],
        plastic_strain=b('plastic_strain_compression_positive'),equivalent_plastic_strain=b('equivalent_plastic_strain'),plastic_work=b('plastic_work_kpa'),
        material_return=[0 if p['last_return']=='elastic' else 1 for p in materials],active_elements=active,
        jump=v('jump_m'),traction=traction,reference_jump=reference,initial_traction=[p.get('birth',{}).get('initial_traction_kpa',[0.,0.,0.]) for p in points],
        plastic_slip=v('plastic_slip_m'),accumulated=v('accumulated_plastic_slip_m'),plastic_dissipation=v('plastic_dissipation_kpa_m'),
        numerical_dissipation=v('numerical_dissipation_kpa_m'),endpoint_work=v('endpoint_work_kpa_m'),elastic_offset=rest,
        initial_energy=[p.get('birth',{}).get('initial_energy_kpa_m',0.) for p in points],
        contact_status=[('open','stick','slip').index(p['status']) for p in points],contact_open=gap>0,opening=np.maximum(gap,0),gap=gap,
        normal_pressure=-qn,shear_traction=tau,elastic_energy=qn**2/(2*s.kn)+(tau*tau).sum(axis=1)/(2*s.kt),
        joint_active=[p['active'] for p in life],joint_tied=[p['tied'] for p in life],generation=[p['generation'] for p in life],
        archived_generation_count=[len(p['archived']) for p in life],
        internal=force.reshape(-1,3),applied_load=load.reshape(-1,3),reaction=(force-load).reshape(-1,3),holding_load=np.asarray(holding).reshape(-1,3))
    for key in ('effective_internal','pore_internal','contact_internal','gravity_load','base_load','surcharge_load','exterior_water_load','joint_water_load'):
        result[key]=np.asarray(d[key]).reshape(-1,3)
    for key in ('joint_water_pressure','joint_face_water_force','exterior_faces','wet_exterior_faces'):result[key]=np.asarray(d[key])
    t,offset=constraint_space(m.mesh,m.case,active,s,life[0]['tied'])
    residual=float(np.linalg.norm(t.T@(force-load))/max(np.linalg.norm(load),1.))
    assert residual<1e-10,(state['stage_index'],state['stage_progress'],residual)
    used=sorted({int(n) for e in np.flatnonzero(active) for n in m.mesh.get('tetrahedra10',m.mesh['tetrahedra'])[e]})
    dofs=(3*np.asarray(used)[:,None]+np.arange(3)).ravel()
    constraint=float(np.max(np.abs(u-t@(t.T@u))[dofs],initial=0.));assert constraint<1e-10
    identity=np.asarray(result['initial_energy'])+result['endpoint_work']-result['elastic_energy']-result['plastic_dissipation']-result['numerical_dissipation']
    work_error=float(np.max(np.abs(identity),initial=0.));assert work_error<1e-10
    # Check native load decomposition, including the full holding vector once.
    physical=sum(result[k] for k in ('gravity_load','base_load','surcharge_load','exterior_water_load','joint_water_load'))
    # Committed release vectors are already scaled by the remaining fraction.
    decomposition=float(np.max(np.abs(physical+result['holding_load']-result['applied_load'])))
    assert decomposition<1e-10,decomposition
    return {k:np.asarray(v) for k,v in result.items()},dict(relative_force_residual=residual,load_decomposition_error_kn=decomposition,constraint_error_m=constraint,local_work_identity_error_kpa_m=work_error)


class Sink:
    def __init__(self,stream,work=True):
        self.stream=stream;self.previous=None;self.work=0.;self.support_work=0.;self.include_work=work

    def placement(self,u):
        # A change of the new nodes' CAD/material reference is not mechanical work.
        self.previous=np.asarray(u).copy()

    def __call__(self,row,fields):
        if self.include_work:
            previous=np.zeros_like(fields['u']) if self.previous is None else self.previous
            du=fields['u']-previous
            fields['nodal_load_work_increment']=fields['applied_load']*du
            fields['nodal_reaction_work_increment']=fields['reaction']*du
            self.work+=float(fields['nodal_load_work_increment'].sum())
            self.support_work+=float(fields['nodal_reaction_work_increment'].sum())
            fields['load_work_cumulative']=np.asarray(self.work)
            fields['reaction_work_cumulative']=np.asarray(self.support_work)
            self.previous=fields['u'].copy()
        self.stream.write(json.dumps(dict(row=row,fields={k:np.asarray(v).tolist() for k,v in fields.items()}),allow_nan=False,separators=(',',':'))+'\n')
        self.stream.flush()


def load_fields(path):
    with path.open(encoding='utf-8') as stream:
        for line in stream:
            v=json.loads(line);yield v['row'],{k:np.asarray(x) for k,x in v['fields'].items()}


def compare(a,b):
    stats={};stages={};count=0
    def add(store,key,x,y):
        assert x.shape==y.shape and np.isfinite(x).all() and np.isfinite(y).all(),key
        if key in CATEGORIES:
            v=store.setdefault(key,dict(mismatches=0,components=0));v['mismatches']+=int(np.sum(x!=y));v['components']+=x.size
        else:
            v=store.setdefault(key,dict(difference_squared=0.,reference_squared=0.,maximum_absolute=0.,components=0))
            d=x-y;v['difference_squared']+=float(np.square(d).sum());v['reference_squared']+=float(np.square(y).sum());v['components']+=x.size
            v['maximum_absolute']=max(v['maximum_absolute'],float(np.max(np.abs(d),initial=0.)))
    for left,right in itertools.zip_longest(load_fields(a),load_fields(b)):
        assert left is not None and right is not None,'missing accepted state'
        ra,fa=left;rb,fb=right
        assert ra==rb,(ra,rb)
        assert set(fa)==set(fb),(set(fa)-set(fb),set(fb)-set(fa))
        for key in fa:
            add(stats,key,fa[key],fb[key]);add(stages.setdefault(ra['stage'],{}),key,fa[key],fb[key])
        count+=1
    for store in [stats]+list(stages.values()):
        for v in store.values():
            if 'difference_squared' in v:
                v['absolute_l2']=v['difference_squared']**.5
                v['relative_l2']=v['absolute_l2']/max(v['reference_squared']**.5,1e-12)
    return dict(accepted_pairs=count,field_count=len(stats),all_fields=stats,stages=stages,
                numerical_acceptance=False,reason='all-field numerical thresholds are not invented from these observations')
