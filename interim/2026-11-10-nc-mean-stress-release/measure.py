"""Remeasure every saved field/system; incomplete evidence never passes."""
import argparse
import json
from pathlib import Path
import numpy as np
import sys
import gzip

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
from mean_stress_release_control import fractions, partition, load as force_path
from make_free_contact import make
from linear_evidence_codec import read_rows as binary_rows
from bounded_mean_release_control import partition as bounded_partition
import joint_reference
import probe_low_preload_nc as low


def rows(p):
    binary=Path(str(p)+'.f64.gz')
    if binary.exists():
        yield from binary_rows(binary)
        return
    opener=p.open if p.exists() else lambda:gzip.open(str(p)+'.gz','rt',encoding='utf-8')
    with opener() as f:
        for line in f:yield json.loads(line)


def read(p):return json.loads(p.read_text())


def audit(folder,allow_platform_float_differences=False):
    c=read(folder/'controls.json')
    if (folder/'report.json').exists():
        r=read(folder/'report.json')
    else:
        watchdog=read(folder/'external-watchdog.json')
        r=dict(passed=False,static_limit_confirmed=False,production_adopted=False,native_committed_increments=0,
               full_removal_reached=False,stop_reason=watchdog['reason'],origin='external watchdog, solver report absent',watchdog=watchdog)
    if r['passed'] or r.get('static_limit_confirmed') or r.get('production_adopted') or r['native_committed_increments']:
        raise ValueError('diagnostic falsely certifies production or stability')
    B=np.array(c['B']);w=np.array(c['weights']);R=np.array(c['transform']);pi=np.array(c['active_gp'])
    initial=read(folder/'initial.json')['payload'];native=c['binding']['map']=='native'
    old=initial['states'];active=np.array(c['active_elements']);contact=initial['contacts']
    final=np.array(c['final_load_kn']);h=np.array(c['horizontal_holding_kn']);v=np.array(c['vertical_holding_kn'])
    mesh,ground,case=make(mcc=True,pc=100.,preload_kpa=5.,force=(1.,.5,0.))
    joints=joint_reference.System(mesh,case['joint_system'])
    if any(active[p['owner_a']] and active[p['owner_b']] for p in joints.points):
        raise ValueError('this diagnostic requires all contact owners cut')
    if native:
        joints.audit(dict(joints=contact,displacement_m=initial['displacement_m'],active_elements=active))
    else:
        q=np.array(contact['traction_kpa']);jump=np.array(contact['jump_m']);normal=joints.normal
        qn=np.einsum('qa,qa->q',q,normal);gn=np.einsum('qa,qa->q',jump,normal);tau=q-qn[:,None]*normal
        energy=qn**2/(2*joints.kn)+np.einsum('qa,qa->q',tau,tau)/(2*joints.kt)
        if not np.allclose(qn,joints.kn*np.minimum(gn,0),atol=1e-7,rtol=1e-9) or np.any(np.linalg.norm(tau,axis=1)>joints.mu*np.maximum(-qn,0)+1e-7):
            raise ValueError('independent frozen contact law differs')
        if not np.allclose(energy,contact['elastic_energy_kpa_m'],atol=1e-9,rtol=1e-8):
            raise ValueError('independent frozen elastic contact energy differs')
    seed_path = (folder.parent/'prepared/origin/prior/fresh-head-midpoint/job/recovery-checkpoint-0.json' if native
                 else folder.parent/'prepared'/('reference800' if folder.name.endswith('-pre800') else 'reference400')/'unload.json')
    seed=read(seed_path)
    if native:
        seed=next(f['state'] for f in seed['run']['frames'] if f['name']=='unload')
        source_stress=np.array([seed['material_states'][i]['effective_stress_compression_positive_kpa'] for i in pi])
    else:
        source_stress=np.array(seed['material_state'])[active].reshape(-1,14)[:,:6]
    if not np.array_equal(source_stress,np.array(c['split_source_stress_kpa'])):
        raise ValueError('split did not originate in own pre-cut history')
    bounded=c['binding'].get('force_decomposition')=='sign_preserving_bounded_mean'
    split=bounded_partition if bounded else partition
    first,second=split(mesh,active,source_stress,np.array(c['holding_kn']),c['binding']['order'])
    partition_difference=float(max(np.max(np.abs(first-h)),np.max(np.abs(second-v))))
    if partition_difference and not allow_platform_float_differences:
        raise ValueError('saved mean-stress/complement partition differs from own history')
    if not c['final_load_exact'] or not np.array_equal(force_path(final,c['holding_kn'],h,v,1.),final):
        raise ValueError('final endpoint differs')
    maximum_force=0.;samples=0;frozen=True;plastic_changed=set();last=None
    materials={x['id']:x['mcc'] for x in ground['materials']}
    m2=np.concatenate([np.full(27,materials[mesh['material_ids'][ei]]['critical_state_m']**2) for ei in np.flatnonzero(active)])
    minimum_mean_pressure=float('inf');minimum_pc=float('inf');minimum_void=float('inf');maximum_positive_yield=0.
    maximum_mean_pressure=0.;maximum_pc=0.;maximum_extra_component=0.;maximum_opposite_component=0.;extra_components=set();original_holding=np.array(c['holding_kn'])
    maximum_displacement=0.;maximum_compatible_strain=0.;reduce_u=np.linalg.pinv(R)
    def fields(state):
        if native:
            return np.array([state[i]['effective_stress_compression_positive_kpa'] for i in pi]),np.array([state[i]['plastic_strain_compression_positive'] for i in pi])
        a=np.array(state)[active].reshape(-1,14);return a[:,:6],a[:,8:]
    _,ep0=fields(old)
    for x in rows(folder/'samples.jsonl'):
        contacts=x['contacts']
        contact_count=len(contacts['points']) if native else len(contacts['jump_m'])
        if x['fractions']!=list(fractions(x['coordinate'])) or np.array(x['displacement_m']).shape!=(54,3) or contact_count!=32:
            raise ValueError('coordinate/full 3D field layout differs')
        if (native and len(x['states'])!=324) or (not native and np.array(x['states']).shape!=(12,27,14)):
            raise ValueError('full material field missing')
        stress,ep=fields(x['states'])
        if native:
            pc=np.array([x['states'][i]['preconsolidation_pressure_kpa'] for i in pi]);void=np.array([x['states'][i]['void_ratio'] for i in pi])
        else:
            state=np.array(x['states'])[active].reshape(-1,14);pc=np.exp(state[:,6]);void=np.exp(state[:,7])-1
        mean=stress[:,:3].mean(axis=1);dev=stress.copy();dev[:,:3]-=mean[:,None]
        q2=1.5*np.square(dev[:,:3]).sum(axis=1)+3*np.square(dev[:,3:]).sum(axis=1)
        yf=(q2+m2*mean*(mean-pc))/(m2*pc**2)
        if not np.isfinite(stress).all() or np.any(mean<=0) or np.any(pc<=0) or np.any(void<=0):raise ValueError('saved active MCC state left its existing domain')
        if not native and np.max(np.maximum(yf,0.))>2e-6:raise ValueError('saved independent state violates existing yield drift gate')
        maximum_mean_pressure=max(maximum_mean_pressure,float(mean.max()));maximum_pc=max(maximum_pc,float(pc.max()))
        minimum_mean_pressure=min(minimum_mean_pressure,float(mean.min()));minimum_pc=min(minimum_pc,float(pc.min()));minimum_void=min(minimum_void,float(void.min()));maximum_positive_yield=max(maximum_positive_yield,float(np.maximum(yf,0.).max()))
        a,b=x['fractions'];load=force_path(final,c['holding_kn'],h,v,x['coordinate'])
        remaining=original_holding if x['coordinate']==0. else (1-a)*h+(1-b)*v
        extra=np.maximum(np.abs(remaining)-np.abs(original_holding),0.)
        opposite=np.maximum(-remaining*np.sign(original_holding),0.)
        maximum_extra_component=max(maximum_extra_component,float(extra.max()));maximum_opposite_component=max(maximum_opposite_component,float(opposite.max()))
        extra_components.update(np.flatnonzero((extra>0.)|(opposite>0.)).tolist())
        if bounded and (np.any(opposite>0.) or np.any(np.abs(remaining)>np.nextafter(np.abs(original_holding),np.inf))):raise ValueError('bounded force path adds original holding (beyond float64 rounding)')
        u=np.array(x['displacement_m']);maximum_displacement=max(maximum_displacement,float(np.linalg.norm(u,axis=1).max()))
        compatible=-np.einsum('qai,i->qa',B,reduce_u@u.ravel())
        maximum_compatible_strain=max(maximum_compatible_strain,float(np.abs(compatible).sum(axis=1).max()))
        f=-np.einsum('qai,qa,q->i',B,stress,w)
        residual=float(np.linalg.norm(R.T@load-f)/max(np.linalg.norm(load),1.))
        if residual>=1e-8:raise ValueError('saved field violates unchanged force gate')
        maximum_force=max(maximum_force,residual)
        if native:
            frozen=all(x['states'][i]==old[i] for i in range(len(old)) if i not in pi)
        else:frozen=bool(np.array_equal(np.array(x['states'])[~active],np.array(old)[~active]))
        if not frozen or x['contacts']!=contact:raise ValueError('removed history changed')
        plastic_changed.update(np.flatnonzero(np.any(ep!=ep0,axis=1)).tolist())
        samples+=1;last=x
    maximum_linear=0.;systems=0
    for x in rows(folder/'linear-systems.jsonl'):
        A=np.array(x['matrix']);b=np.array(x['rhs']);solution=np.array(x['solution'])
        residual=float(np.linalg.norm(A@solution-b)/max(np.linalg.norm(b),1e-30))
        # Rejected systems are evidence, not accepted systems. Keep their true
        # residual even if above the unchanged gate; never reclassify the stop.
        maximum_linear=max(maximum_linear,residual);systems+=1
    attempts=list(rows(folder/'attempts.jsonl'))
    for attempt in attempts:
        is_rate=attempt.get('events') and attempt['events'][0].get('method','').startswith('global_load_rk4')
        if attempt['accepted'] and ((not is_rate and len(attempt['events'])<2) or any(e['true_linear_residual']>1e-10 for e in attempt['events'])):
            raise ValueError('accepted Newton path violates unchanged linear gate')
    if any(not x.get('rejected_committed_state_unchanged',True) for x in attempts):raise ValueError('rejected trial changed state')
    if r['passed'] or r.get('static_limit_confirmed') or r.get('production_adopted') or r['native_committed_increments']:
        raise ValueError('diagnostic falsely certifies production or stability')
    trial_count=0;trial_maximum_force=0.
    if (folder/'rate-trials.jsonl').exists() or Path(str(folder/'rate-trials.jsonl')+'.gz').exists():
        for trial in rows(folder/'rate-trials.jsonl'):
            if not trial['outer_trial_not_yet_committed'] or trial['contacts']!=contact:
                raise ValueError('rate-trial status or frozen contact differs')
            state=np.array(trial['states'])
            if not np.array_equal(state[~active],np.array(old)[~active]):raise ValueError('rate trial changed removed material')
            stress=state[active].reshape(-1,14)[:,:6]
            f=-np.einsum('qai,qa,q->i',B,stress,w)
            load=force_path(final,c['holding_kn'],h,v,trial['coordinate'])
            residual=float(np.linalg.norm(R.T@load-f)/max(np.linalg.norm(load),1.))
            trial_maximum_force=max(trial_maximum_force,residual);trial_count+=1
    release=last['fractions'];remaining=(1-release[0])*h+(1-release[1])*v
    full=release==[1.,1.] and r['full_removal_reached']
    return dict(report=r,grouping='mean-stress nodal force/complement; legacy horizontal/vertical keys mean first/second',partition_reassembly_max_absolute_difference_kn=partition_difference,decomposition_roundoff_norm_kn=float(np.linalg.norm(h+v-np.array(c['holding_kn']))),saved_fields=samples,linear_systems=systems,maximum_reassembled_force_residual=maximum_force,
        maximum_remeasured_linear_residual_including_rejections=maximum_linear,internal_rate_trial_fields=trial_count,maximum_internal_rate_trial_force_residual_including_rejections=trial_maximum_force,removed_material_and_contacts_frozen=frozen,
        active_plastic_points_changed=len(plastic_changed),rejected_attempts=sum(not a['accepted'] for a in attempts),
        initial_coordinate=initial['coordinate'],last_coordinate=last['coordinate'],last_fractions=release,
        remaining_holding_norm_kn=float(np.linalg.norm(remaining)),full_release=full,full_removal_lambda=1. if full else None,nonproportional_partial_coordinate_is_not_lambda=True,force_decomposition_bounded=bounded,maximum_active_mean_pressure_kpa=maximum_mean_pressure,maximum_active_pc_kpa=maximum_pc,maximum_extra_holding_component_kn=maximum_extra_component,maximum_opposite_holding_component_kn=maximum_opposite_component,original_holding_component_bound_exceedances=len(extra_components),minimum_active_mean_pressure_kpa=minimum_mean_pressure,minimum_active_pc_kpa=minimum_pc,minimum_active_void_ratio=minimum_void,maximum_active_normalized_positive_yield=maximum_positive_yield,
        maximum_nodal_displacement_from_CAD_m=maximum_displacement,
        maximum_compatible_strain_l1_from_zero_displacement=maximum_compatible_strain,
        lambda_zero_newton=next((x for x in attempts if x.get('label')=='lambda-zero'),None),
        full_release_freeze=read(folder/'full-release.json') if (folder/'full-release.json').exists() else None)


def vector(x,field,native):
    if field=='displacement':return np.array(x['displacement_m'])
    if field.startswith('contact_'):
        key={'contact_jump':'jump_m','contact_traction':'traction_kpa','contact_slip':'plastic_slip_m','contact_elastic_energy':'elastic_energy_kpa_m'}[field]
        if not native:return np.array(x['contacts'][key])
        points=x['contacts']['points']
        if field!='contact_elastic_energy':return np.array([p[key] for p in points])
        _,_,case=make(mcc=True,pc=100.,preload_kpa=5.,force=(1.,.5,0.))
        material=case['joint_system']['joints'][0]['material'];kn=material['normal_stiffness_kpa_m'];kt=material['tangential_stiffness_kpa_m']
        q=np.array([p['traction_kpa'] for p in points]);n=np.array([p['normal'] for p in points]);qn=np.einsum('qa,qa->q',q,n);tau=q-qn[:,None]*n
        return qn**2/(2*kn)+np.einsum('qa,qa->q',tau,tau)/(2*kt)
    if native:
        key={'stress':'effective_stress_compression_positive_kpa','plastic':'plastic_strain_compression_positive','pc':'preconsolidation_pressure_kpa','void':'void_ratio'}[field]
        return np.array([s[key] for s in x['states']])
    a=np.array(x['states']).reshape(-1,14)
    if field=='stress':return a[:,:6]
    if field=='plastic':return a[:,8:]
    if field=='pc':return np.exp(a[:,6])
    return np.exp(a[:,7])-1


def compare(a,b):
    aa_paths=a if isinstance(a,list) else [a];bb_paths=b if isinstance(b,list) else [b]
    ca=read(aa_paths[0]/'controls.json');cb=read(bb_paths[0]/'controls.json')
    for folders,control in [(aa_paths,ca),(bb_paths,cb)]:
        for p in folders[1:]:
            other=read(p/'controls.json')
            if other['binding']!=control['binding'] or other['holding_kn']!=control['holding_kn'] or other['final_load_kn']!=control['final_load_kn']:
                raise ValueError('continued history model/control differs')
    native=ca['binding']['map']=='native'
    native_b=cb['binding']['map']=='native'
    aa={x['coordinate']:x for p in aa_paths for x in rows(p/'samples.jsonl')};count=0;maximum={};absolute={};end=None;load_difference=0.;reduced_load_difference=0.
    bb={x['coordinate']:x for p in bb_paths for x in rows(p/'samples.jsonl')}
    Ra=np.array(ca['transform']);Rb=np.array(cb['transform'])
    def load(c,x):
        h,v=x['fractions']
        return force_path(c['final_load_kn'],c['holding_kn'],c['horizontal_holding_kn'],c['vertical_holding_kn'],x['coordinate'])
    for y in bb.values():
        if y['coordinate'] not in aa:continue
        x=aa[y['coordinate']];count+=1;end=y['coordinate']
        la,lb=load(ca,x),load(cb,y)
        load_difference=max(load_difference,float(np.linalg.norm(la-lb)/max(np.linalg.norm(lb),1.)))
        # Full nodal load includes differing support/reaction conventions.
        # Compare generalized force in each model's compatible free basis too.
        reduced_load_difference=max(reduced_load_difference,float(np.linalg.norm(Ra.T@la-Rb.T@lb)/max(np.linalg.norm(Rb.T@lb),1.)))
        for f in ['displacement','stress','plastic','pc','void','contact_jump','contact_traction','contact_slip','contact_elastic_energy']:
            av=vector(x,f,native);bv=vector(y,f,native_b)
            d=float(np.linalg.norm(av-bv)/max(np.linalg.norm(bv),1e-30))
            maximum[f]=max(maximum.get(f,0.),d)
            absolute[f]=max(absolute.get(f,0.),float(np.max(np.abs(av-bv))))
    return dict(common_coordinates=count,last_common_coordinate=end,maximum_relative_L2=maximum,maximum_absolute=absolute,
                maximum_load_relative_L2=load_difference,
                maximum_reduced_load_relative_L2=reduced_load_difference,
                holding_difference_norm_kn=float(np.linalg.norm(np.array(ca['holding_kn'])-np.array(cb['holding_kn']))),
                reduced_holding_difference_norm_kn=float(np.linalg.norm(Ra.T@np.array(ca['holding_kn'])-Rb.T@np.array(cb['holding_kn']))),
                complete_history_converged=False,new_acceptance_threshold=False)


def restart_compare(evidence, direct_name, resumed_name):
    direct=evidence/direct_name;resumed=evidence/resumed_name
    aa={x['coordinate']:x for x in rows(direct/'samples.jsonl')};count=0
    for x in rows(resumed/'samples.jsonl'):
        y=aa[x['coordinate']]
        xx=dict(x);yy=dict(y)
        if count==0:
            xx.pop('reduced_increment_m',None);yy.pop('reduced_increment_m',None)
        if xx!=yy:raise ValueError('restart accepted full field differs')
        count+=1;last=x
    attempts=list(rows(resumed/'attempts.jsonl'));lo=min(x['step'] for x in attempts);hi=max(x['step'] for x in attempts)
    direct_attempts=[x for x in rows(direct/'attempts.jsonl') if lo<=x.get('step',-1)<=hi]
    rate_offsets=[]
    for group in [direct_attempts,attempts]:
        for attempt in group:
            for event in attempt['events']:
                if 'diagnostics' in event:
                    event['diagnostics']=dict(event['diagnostics'])
    for a,b in zip(direct_attempts,attempts):
        for ea,eb in zip(a['events'],b['events']):
            if 'diagnostics' in ea:
                da=ea['diagnostics'];db=eb['diagnostics']
                rate_offsets.append(da.pop('rate_assemblies')-db.pop('rate_assemblies'))
    if attempts!=direct_attempts:raise ValueError('restart actual attempt path differs')
    if rate_offsets and len(set(rate_offsets))!=1:raise ValueError('restart informational rate counter offset varies')
    def systems(p):
        for x in rows(p/'linear-systems.jsonl'):
            if x['label'].startswith('step-') and lo<=int(x['label'][5:])<=hi:yield x
    from itertools import zip_longest
    systems_count=0
    for a,b in zip_longest(systems(direct),systems(resumed)):
        if a!=b:raise ValueError('restart actual linear systems differ')
        systems_count+=1
    return dict(exact_common_full_fields=count,entire_accepted_rows_exact=True,
                attempt_path_exact_except_recorded_informational_rate_counter=bool(rate_offsets),
                entire_attempt_path_exact=not rate_offsets,informational_rate_assemblies_constant_offset=rate_offsets[0] if rate_offsets else None,
                actual_linear_systems_exact=True,compared_actual_linear_systems=systems_count,
                completion_verified=last['coordinate']==1.)


def post_release(evidence,name):
    folder=evidence/(name+'-post')
    if not (folder/'post-full-release.json').exists():return None
    import hashlib
    result=read(folder/'post-full-release.json');last=read(evidence/name/'restart.json')['payload'];c=read(evidence/name/'controls.json')
    if result['source_checkpoint_sha256']!=hashlib.sha256((evidence/name/'restart.json').read_bytes()).hexdigest():raise ValueError('post-release source differs')
    B=np.array(c['B']);w=np.array(c['weights']);R=np.array(c['transform']);pi=np.array(c['active_gp']);native=c['binding']['map']=='native';active=np.array(c['active_elements'])
    maxima=[]
    for row in result['full_release_zero_increment_checks']:
        if row['states']!=last['states'] or row['contacts']!=last['contacts'] or row['displacement_m']!=last['displacement_m'] or np.any(row['holding_load_kn']):raise ValueError('post-release frozen full field differs')
        stress=(np.array([row['states'][i]['effective_stress_compression_positive_kpa'] for i in pi]) if native else np.array(row['states'])[active].reshape(-1,14)[:,:6])
        f=-np.einsum('qai,qa,q->i',B,stress,w);load=np.array(c['final_load_kn'])
        res=float(np.linalg.norm(R.T@load-f)/max(np.linalg.norm(load),1.))
        if res>=1e-8:raise ValueError('post-release force differs')
        maxima.append(res)
    return dict(actual_zero_increment_checks=len(maxima),all_material_contact_displacement_fields_exact=True,maximum_force_residual=max(maxima),holding_load_exactly_zero=True)


def prehistory(evidence):
    import gzip,hashlib
    from low_preload_nc_reference import ClosedStickSystem,expected_inputs
    import mcc_rate_fem as rate
    import fem_reference as fem
    from verify_high_order import physical
    mesh,ground,case=expected_inputs();result={}
    for count in [400,800]:
        folder=evidence/'prepared'/f'reference{count}'
        history=folder.parent/(folder.name+'-history.jsonl.gz')
        if not history.exists():continue
        full=np.ones(len(mesh['tetrahedra']),dtype=bool)
        cut=np.array(mesh['region_ids'])!=next(r['id'] for r in mesh['regions'] if r['body']=='UpperClay')
        systems={True:ClosedStickSystem(mesh,ground,full,case),False:ClosedStickSystem(mesh,ground,cut,case)}
        seed=rate.seed(mesh,ground);u0=np.zeros(3*len(mesh['nodes']));initial_before=systems[True].total_force(seed.reshape(-1,14),u0)
        base={key:physical(mesh,ground,sys.active,0.,True,joint_system=sys.joint)[0].ravel() for key,sys in systems.items()}
        boundary=fem.conditions_load(mesh,case['boundary'],full,skip_removed_surfaces=True).ravel()
        unload=read(folder/'unload.json');unload_state=np.array(unload['material_state']);unload_u=np.array(unload['displacement_m']).ravel()
        cut_before=systems[False].total_force(unload_state[cut].reshape(-1,14),unload_u)
        counts={};maximum=0.;frozen=True
        with gzip.open(history,'rt') as f:
            for line in f:
                row=json.loads(line);stage=row['stage'];t=row['progress'];intact=stage!='excavate_upper';sys=systems[intact]
                state=np.array(row['material_state']);u=np.array(row['displacement_m']).ravel()
                if row['active_elements']!=sys.active.tolist():raise ValueError('independent prehistory active set differs')
                if stage=='initial_equilibrated':external=initial_before+t*(base[True]-initial_before)
                elif stage=='load':external=base[True]+t*boundary
                elif stage=='unload':external=base[True]+(1-t)*boundary
                else:external=cut_before+t*(base[False]-cut_before)
                force=sys.total_force(state[sys.active].reshape(-1,14),u)
                residual=float(np.linalg.norm(sys.reduce(force-external))/max(np.linalg.norm(external),1.))
                if residual>=1e-8:raise ValueError('independent prehistory force violates existing gate')
                maximum=max(maximum,residual);counts[stage]=counts.get(stage,0)+1
                evaluated=sys.contacts(u)
                if intact:
                    for key in evaluated:
                        if not np.allclose(evaluated[key],row['contacts'][key],atol=1e-12,rtol=1e-12):raise ValueError('independent prehistory contact field differs')
                else:
                    frozen=frozen and np.array_equal(state[~cut],unload_state[~cut]) and row['contacts']==unload['contacts']
                    if not frozen:raise ValueError('unused subsequent proportional history is not frozen')
        result[str(count)]=dict(input_only=True,all_accepted_material_and_contact_fields_saved=True,accepted_rows_by_stage=counts,
            maximum_reassembled_force_residual=maximum,initial_load_unload_complete=all((folder/(n+'.json')).exists() for n in ['initial_equilibrated','load','unload']),
            subsequent_proportional_excavation_not_used_for_new_path=True,subsequent_proportional_report=read(folder/'reference.json')['failure_reason'],
            history_gzip_sha256=hashlib.sha256(history.read_bytes()).hexdigest())
    return result


def measure(evidence,allow_platform_float_differences=False):
    paths=sorted(p for p in evidence.iterdir() if p.is_dir() and (p/'controls.json').exists() and ((p/'report.json').exists() or (p/'external-watchdog.json').exists()))
    result=dict(passed=False,full_removal_accepted=False,history_converged=False,static_limit_confirmed=False,production_adopted=False,
                native_committed_increments=0,paths={p.name:audit(p,allow_platform_float_differences) for p in paths},step_sensitivity={},restart={},cross_map={},material_sensitivity={},prehistory_sensitivity={},post_full_release={})
    for order in ['mean','complement']:
        for kind,suffix in [('native',''),('independent',''),('independent','-regrowth'),('independent','-rate')]:
            a=evidence/f'{kind}-{order}128{suffix}';b=evidence/f'{kind}-{order}256{suffix}'
            if a in paths and b in paths:result['step_sensitivity'][kind+'-'+order+suffix]=compare(a,b)
        for step in ['128','256']:
            a=evidence/f'native-{order}{step}';b=evidence/f'independent-{order}{step}-rate'
            if a in paths and b in paths:result['cross_map'][order+step+'-global-rate']=compare(a,b)
    a=evidence/'independent-complement128-rate'
    for suffix,key in [('-m8','material_sensitivity'),('-pre800','prehistory_sensitivity')]:
        b=evidence/('independent-complement128-rate'+suffix)
        if a in paths and b in paths:result[key][suffix[1:]]=compare(a,b)
    for direct,resumed in [('native-complement128','native-complement-resumed'),('native-complement256','native-complement256-resumed'),('independent-complement128-rate','independent-complement-rate-resumed'),('bounded-native-complement128','bounded-native-complement-resumed'),('bounded-independent-complement128','bounded-independent-complement-resumed')]:
        if evidence/direct in paths and evidence/resumed in paths:result['restart'][direct]=restart_compare(evidence,direct,resumed)
    result['bounded_withdrawal']={}
    for order in ['mean','complement']:
        for kind in ['native','independent']:
            a=evidence/f'bounded-{kind}-{order}128';b=evidence/f'bounded-{kind}-{order}256'
            if a in paths and b in paths:result['step_sensitivity']['bounded-'+kind+'-'+order]=compare(a,b)
        for step in ['128','256']:
            a=evidence/f'bounded-native-{order}{step}';b=evidence/f'bounded-independent-{order}{step}'
            if a in paths and b in paths:result['cross_map']['bounded-'+order+step]=compare(a,b)
    for suffix,key in [('-m8','material_sensitivity'),('-pre800','prehistory_sensitivity')]:
        a=evidence/'bounded-independent-complement128';b=evidence/('bounded-independent-complement128'+suffix)
        if a in paths and b in paths:result[key]['bounded'+suffix]=compare(a,b)
    bounded_paths=[p for p in paths if p.name.startswith('bounded-') and p.name[-3:] in ['128','256']]
    result['bounded_withdrawal']['all_eight_paths_full_release']=len(bounded_paths)==8 and all(result['paths'][p.name]['full_release'] for p in bounded_paths)
    result['bounded_withdrawal']['complete_orders']=[order for order in ['mean','complement'] if all((evidence/f'bounded-{kind}-{order}{step}') in paths and result['paths'][f'bounded-{kind}-{order}{step}']['full_release'] for kind in ['native','independent'] for step in [128,256])]
    result['bounded_withdrawal']['unbounded_auxiliary_force_paths_not_transferred_to_goal']=True
    tail=evidence/'independent-complement-tail'
    coarse=evidence/'independent-complement128';fine=evidence/'independent-complement256'
    if tail in paths:
        result['continued_history_sensitivity']=dict(independent_complement=compare([coarse,tail],fine))
        result['cross_map']['complement128-local-strain-continued']=compare(evidence/'native-complement128',[coarse,tail])
    for step in ['128','256']:
        a=evidence/f'independent-complement{step}-regrowth';b=evidence/f'independent-complement{step}-rate'
        if a in paths and b in paths:result['cross_map']['local-strain-to-global-rate-'+step]=compare(a,b)
    for folder in paths:
        post=post_release(evidence,folder.name)
        if post:result['post_full_release'][folder.name]=post
    result['independent_prehistories']=prehistory(evidence)
    prefix=evidence/'fresh-native-prefix'
    if (prefix/'verification.json').exists():
        record=read(prefix/'verification.json');run=read(prefix/'result.json')
        fresh=next(f['state'] for f in run['frames'] if f['name']=='unload')
        before=read(evidence/'prepared/origin/prior/fresh-head-midpoint/job/recovery-checkpoint-0.json')['run']
        if fresh!=next(f['state'] for f in before['frames'] if f['name']=='unload'):raise ValueError('fresh native unload differs')
        result['fresh_native_prefix']=record
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--evidence',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.write_text(json.dumps(measure(a.evidence),indent=2,allow_nan=False)+'\n')
