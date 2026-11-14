"""Native observations and exhaustive full-field differences, no fitted gates."""
import itertools
import json
import numpy as np
from reference import STATES

CATEGORIES={'material_return','contact_status','contact_open','joint_active','joint_tied',
            'generation','archived_generation_count','anchor_0_taut','frozen_status'}


def native_fields(model,state,diagnostic,reports):
    """Read actual saved material/contact/structure values, not a fresh solve.

    Strain, gap and energy are named derived observations. Internal and target
    are the actual converged production diagnostics, not independent reassembly.
    Structure actions are production local reports at its accepted solution.
    """
    s=model.system;j=state['joints']['points'];u=np.asarray(state['displacement_m']).ravel()
    v=lambda k:np.asarray([p[k] for p in j])
    materials=state['material_states'];b=lambda k:np.asarray([p[k] for p in materials])
    traction=v('traction_kpa');qn=np.einsum('qa,qa->q',traction,s.normal);tau=traction-qn[:,None]*s.normal
    reference=np.array([p.get('birth',{}).get('jump_m',[0.,0.,0.]) for p in j])
    rest=np.array([p.get('birth',{}).get('elastic_offset_m',[0.,0.,0.]) for p in j])
    initial=np.array([p.get('birth',{}).get('initial_traction_kpa',[0.,0.,0.]) for p in j])
    gap=np.einsum('qa,qa->q',v('jump_m')-rest,s.normal)
    force=np.asarray(diagnostic['internal']).reshape(-1,3)
    load=np.asarray(state['applied_nodal_load_kn']).reshape(-1,3)
    np.testing.assert_array_equal(np.asarray(diagnostic['target']).ravel(),load.ravel())
    strain=np.concatenate([np.einsum('qai,i->qa',op['b'],u[ids]) for op,ids,*_ in model.ops])
    life=state['joints']['lifecycle']
    holding=state.get('geostatic_release_load_kn') if state['stage_index']==0 else state.get('transition_release_load_kn')
    if not holding:holding=np.zeros(u.size)
    result=dict(u=u.reshape(-1,3),stress=b('effective_stress_compression_positive_kpa'),strain=strain,void=b('void_ratio'),
        plastic_strain=b('plastic_strain_compression_positive'),equivalent_plastic_strain=b('equivalent_plastic_strain'),
        plastic_work=b('plastic_work_kpa'),material_return=np.array([0 if p['last_return']=='elastic' else 1 for p in materials]),
        jump=v('jump_m'),traction=traction,plastic_slip=v('plastic_slip_m'),accumulated=v('accumulated_plastic_slip_m'),
        plastic_dissipation=v('plastic_dissipation_kpa_m'),numerical_dissipation=v('numerical_dissipation_kpa_m'),endpoint_work=v('endpoint_work_kpa_m'),
        elastic_energy=qn**2/(2*s.kn)+(tau*tau).sum(axis=1)/(2*s.kt),gap=gap,opening=np.maximum(gap,0),
        contact_open=gap>0,contact_status=np.array([STATES.index(p['status']) for p in j]),
        initial_traction=initial,reference_jump=reference,elastic_offset=rest,
        initial_energy=np.array([p.get('birth',{}).get('initial_energy_kpa_m',0.) for p in j]),
        generation=np.array([p['generation'] for p in life]),joint_active=np.array([p['active'] for p in life]),
        joint_tied=np.array([p['tied'] for p in life]),archived_generation_count=np.array([len(p['archived']) for p in life]),
        internal=force,applied_load=load,reaction=force-load,holding_load=np.asarray(holding).reshape(-1,3),
        normal_pressure=-qn,shear_traction=tau)
    frozen={k:np.zeros((len(j),3)) for k in ('jump','traction','plastic','initial','reference','rest')}
    frozen.update({k:np.zeros(len(j)) for k in ('accumulated','pd','nd','work','initial_energy','status')})
    for lifecycle in life:
        assert len(lifecycle['archived'])<=1,'this bounded specimen has two live generations and one archive'
        for archive in lifecycle['archived']:
            for pi,p in archive['points']:
                for key,raw in dict(jump='jump_m',traction='traction_kpa',plastic='plastic_slip_m',accumulated='accumulated_plastic_slip_m',
                                    pd='plastic_dissipation_kpa_m',nd='numerical_dissipation_kpa_m',work='endpoint_work_kpa_m').items():frozen[key][pi]=p[raw]
                for key,raw in dict(initial='initial_traction_kpa',reference='jump_m',rest='elastic_offset_m',initial_energy='initial_energy_kpa_m').items():frozen[key][pi]=p['birth'][raw]
                frozen['status'][pi]=STATES.index(p['status'])
    for key,value in frozen.items():result['frozen_'+key]=value
    if state.get('beams'):
        raw=state['beams'];report=reports['beams']
        result.update(structure_u=raw['node_displacement_m'],structure_rotation=raw['node_rotation_rad'],
                      structure_support_force=report['support_force_kn'],structure_support_moment=report['support_moment_origin_kn_m'],
                      structure_external_force=report['external_force_kn'],structure_external_moment=report['external_moment_origin_kn_m'])
        for member in report['members']:
            for name in ('local_end_actions','local_end_relative_rotation_rad','section_forces_at_start_mid_end'):
                result['beam_'+member['id']+'_'+name]=member[name]
        for p in report.get('plates',[]):
            for name in ('membrane_energy_kn_m','bending_energy_kn_m','drilling_penalty_energy_kn_m'):
                result['plate_'+p['id']+'_'+name]=p[name]
            for name in p['samples'][0]:result['plate_'+p['id']+'_'+name]=[v[name] for v in p['samples']]
    if state.get('anchors'):
        for i,r in enumerate(reports['anchors']['anchors']):
            for name in ('point_displacement_m','point_soil_force_kn','segment_forces_kn','segment_extension_m','elastic_energy_kn_m'):
                result['anchor_'+str(i)+'_'+name]=r[name]
            result['anchor_'+str(i)+'_taut']=state['anchors']['anchors'][i]['segment_taut']
    return {k:np.asarray(v) for k,v in result.items()}


class Sink:
    def __init__(self,stream):
        self.stream=stream;self.previous=None;self.work=0.;self.support_work=0.

    def __call__(self,row,fields):
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
        assert (ra['stage'],ra['progress'])==(rb['stage'],rb['progress']),(ra,rb)
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
