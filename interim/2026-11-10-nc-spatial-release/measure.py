"""Remeasure every saved field/system; incomplete evidence never passes."""
import argparse
import json
from pathlib import Path
import numpy as np
import sys
import gzip

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
from spatial_release_control import fractions, partition
from make_free_contact import make
import joint_reference
import probe_low_preload_nc as low


def rows(p):
    opener=p.open if p.exists() else lambda:gzip.open(str(p)+'.gz','rt',encoding='utf-8')
    with opener() as f:
        for line in f:yield json.loads(line)


def read(p):return json.loads(p.read_text())


def audit(folder):
    c=read(folder/'controls.json');r=read(folder/'report.json')
    B=np.array(c['B']);w=np.array(c['weights']);R=np.array(c['transform']);pi=np.array(c['active_gp'])
    initial=read(folder/'initial.json')['payload'];native=c['binding']['map']=='native'
    old=initial['states'];active=np.array(c['active_elements']);contact=initial['contacts']
    final=np.array(c['final_load_kn']);h=np.array(c['horizontal_holding_kn']);v=np.array(c['vertical_holding_kn'])
    mesh,_,case=make(mcc=True,pc=100.,preload_kpa=5.,force=(1.,.5,0.))
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
    first,second,nodes=partition(mesh['nodes'],np.array(c['holding_kn']),c['binding']['order'])
    if not np.array_equal(first,h) or not np.array_equal(second,v) or nodes!=c['first_group_nodes']:
        raise ValueError('saved spatial partition differs from existing mesh coordinates')
    if not np.array_equal(h+v,np.array(c['holding_kn'])) or not c['final_load_exact']:
        raise ValueError('load decomposition/final endpoint differs')
    maximum_force=0.;samples=0;frozen=True;plastic_changed=set();last=None
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
        stress,ep=fields(x['states']);a,b=x['fractions'];load=final+(1-a)*h+(1-b)*v
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
        if attempt['accepted'] and (len(attempt['events'])<2 or any(e['true_linear_residual']>1e-10 for e in attempt['events'])):
            raise ValueError('accepted Newton path violates unchanged linear gate')
    if any(not x.get('rejected_committed_state_unchanged',True) for x in attempts):raise ValueError('rejected trial changed state')
    if r['passed'] or r.get('static_limit_confirmed') or r.get('production_adopted') or r['native_committed_increments']:
        raise ValueError('diagnostic falsely certifies production or stability')
    release=last['fractions'];remaining=(1-release[0])*h+(1-release[1])*v
    full=release==[1.,1.] and r['full_removal_reached']
    return dict(report=r,grouping='spatial nodes, all XYZ; legacy horizontal/vertical JSON keys mean first/second group',first_group_nodes=nodes,saved_fields=samples,linear_systems=systems,maximum_reassembled_force_residual=maximum_force,
        maximum_remeasured_linear_residual_including_rejections=maximum_linear,removed_material_and_contacts_frozen=frozen,
        active_plastic_points_changed=len(plastic_changed),rejected_attempts=sum(not a['accepted'] for a in attempts),
        initial_coordinate=initial['coordinate'],last_coordinate=last['coordinate'],last_fractions=release,
        remaining_holding_norm_kn=float(np.linalg.norm(remaining)),full_release=full,
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
        return np.array(c['final_load_kn'])+(1-h)*np.array(c['horizontal_holding_kn'])+(1-v)*np.array(c['vertical_holding_kn'])
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


def measure(evidence):
    paths=sorted(p for p in evidence.iterdir() if p.is_dir() and (p/'controls.json').exists() and (p/'report.json').exists())
    result=dict(passed=False,full_removal_accepted=False,history_converged=False,static_limit_confirmed=False,
                production_adopted=False,native_committed_increments=0,paths={p.name:audit(p) for p in paths},step_sensitivity={},restart={},cross_map={})
    for kind in ['native','independent']:
        for order in ['left','center']:
            a=evidence/(kind+'-'+order+'128');b=evidence/(kind+'-'+order+'256')
            if a in paths and b in paths:result['step_sensitivity'][kind+'-'+order]=compare(a,b)
            resumed=evidence/(kind+'-'+order+'-resumed')
            if a in paths and resumed in paths:
                aa={x['coordinate']:x for x in rows(a/'samples.jsonl')};n=0
                for x in rows(resumed/'samples.jsonl'):
                    # Initial restart row has no increment, compare full fields.
                    if x['coordinate'] not in aa:continue
                    y=aa[x['coordinate']]
                    for key in ['states','displacement_m','contacts','fractions','reassembled_relative_force_residual']:
                        if x[key]!=y[key]:raise ValueError('restart full field differs: '+key)
                    xx=dict(x);yy=dict(y)
                    if n==0:
                        xx.pop('reduced_increment_m',None);yy.pop('reduced_increment_m',None)
                    if xx!=yy:raise ValueError('restart accepted row differs')
                    n+=1
                attempts=list(rows(resumed/'attempts.jsonl'));lo=min(x['step'] for x in attempts);hi=max(x['step'] for x in attempts)
                direct=[x for x in rows(a/'attempts.jsonl') if lo<=x.get('step',-1)<=hi]
                if attempts!=direct:raise ValueError('restart full Newton attempt path differs')
                def systems(p):
                    return [x for x in rows(p/'linear-systems.jsonl') if x['label'].startswith('step-') and lo<=int(x['label'][5:])<=hi]
                linear=systems(resumed)
                if linear!=systems(a):raise ValueError('restart actual linear systems differ')
                result['restart'][kind+'-'+order]=dict(exact_common_full_fields=n,entire_accepted_rows_exact=True,entire_attempt_path_exact=True,actual_linear_systems_exact=True,compared_actual_linear_systems=len(linear),completion_verified=False)
    for order in ['left','center']:
        for step in ['128','256']:
            a=evidence/('native-'+order+step);b=evidence/('independent-'+order+step)
            if a in paths and b in paths:result['cross_map'][order+step]=compare(a,b)
    result['material_sensitivity']={}
    a=evidence/'independent-left128';b=evidence/'independent-left128-m8'
    if a in paths and b in paths:result['material_sensitivity']['left4-to-8']=compare(a,b)
    b=evidence/'independent-event-left128'
    if a in paths and b in paths:result['event_floor_sensitivity']=compare(a,b)
    result['continued_history_sensitivity']={}
    for kind,order in [('native','left'),('independent','center')]:
        a=evidence/f'{kind}-{order}128';tail=evidence/f'{kind}-{order}-tail';b=evidence/f'{kind}-{order}256'
        if a in paths and tail in paths and b in paths:
            result['continued_history_sensitivity'][kind+'-'+order]=compare([a,tail],b)
    a=evidence/'native-center128';b=evidence/'independent-center128';tail=evidence/'independent-center-tail'
    if a in paths and b in paths and tail in paths:
        result['cross_map']['center128-continued']=compare(a,[b,tail])
    phase=evidence/'phase-restarts'
    if (phase/'verification.json').exists():
        verification=read(phase/'verification.json')
        for row in verification['comparisons']:
            p=phase/row['phase'];first=read(p/'first/observed-run.json');resumed=read(p/'resumed/observed-run.json');direct=read(p/'direct.json')['run']
            low.check_inputs(p/'first/input')
            low.compare_prefix(resumed,direct)
            if low.deterministic(resumed)!=low.deterministic(direct) or first['last_committed']['stage_index']!=row['stage_index'] or not 0<row['starting_progress']<1:
                raise ValueError('production partial restart differs')
        if verification['passed'] or len(verification['comparisons'])!=4:
            raise ValueError('partial restart scope differs')
        result['production_phase_restarts']=dict(verification=verification,remeasured_entire_run_except_elapsed_seconds=True,completion_verified=False)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--evidence',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.write_text(json.dumps(measure(a.evidence),indent=2,allow_nan=False)+'\n')
