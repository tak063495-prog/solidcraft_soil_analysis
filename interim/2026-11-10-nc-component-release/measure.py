"""Remeasure every saved field/system; incomplete evidence never passes."""
import argparse
import json
from pathlib import Path
import numpy as np


def rows(p):
    with p.open() as f:
        for line in f:yield json.loads(line)


def read(p):return json.loads(p.read_text())


def audit(folder):
    c=read(folder/'controls.json');r=read(folder/'report.json')
    B=np.array(c['B']);w=np.array(c['weights']);R=np.array(c['transform']);pi=np.array(c['active_gp'])
    initial=read(folder/'initial.json')['payload'];native=c['binding']['map']=='native'
    old=initial['states'];active=np.array(c['active_elements']);contact=initial['contacts']
    final=np.array(c['final_load_kn']);h=np.array(c['horizontal_holding_kn']);v=np.array(c['vertical_holding_kn'])
    if not np.array_equal(h+v,np.array(c['holding_kn'])) or not c['final_load_exact']:
        raise ValueError('load decomposition/final endpoint differs')
    maximum_force=0.;samples=0;frozen=True;plastic_changed=set();last=None
    def fields(state):
        if native:
            return np.array([state[i]['effective_stress_compression_positive_kpa'] for i in pi]),np.array([state[i]['plastic_strain_compression_positive'] for i in pi])
        a=np.array(state)[active].reshape(-1,14);return a[:,:6],a[:,8:]
    _,ep0=fields(old)
    for x in rows(folder/'samples.jsonl'):
        stress,ep=fields(x['states']);a,b=x['fractions'];load=final+(1-a)*h+(1-b)*v
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
    if any(not x.get('rejected_committed_state_unchanged',True) for x in attempts):raise ValueError('rejected trial changed state')
    if r['passed'] or r.get('static_limit_confirmed') or r.get('production_adopted') or r['native_committed_increments']:
        raise ValueError('diagnostic falsely certifies production or stability')
    release=last['fractions'];remaining=(1-release[0])*h+(1-release[1])*v
    full=release==[1.,1.] and r['full_removal_reached']
    return dict(report=r,saved_fields=samples,linear_systems=systems,maximum_reassembled_force_residual=maximum_force,
        maximum_remeasured_linear_residual_including_rejections=maximum_linear,removed_material_and_contacts_frozen=frozen,
        active_plastic_points_changed=len(plastic_changed),rejected_attempts=sum(not a['accepted'] for a in attempts),
        initial_coordinate=initial['coordinate'],last_coordinate=last['coordinate'],last_fractions=release,
        remaining_holding_norm_kn=float(np.linalg.norm(remaining)),full_release=full,
        lambda_zero_newton=next((x for x in attempts if x.get('label')=='lambda-zero'),None),
        full_release_freeze=read(folder/'full-release.json') if (folder/'full-release.json').exists() else None)


def vector(x,field,native):
    if field=='displacement':return np.array(x['displacement_m'])
    if native:
        key={'stress':'effective_stress_compression_positive_kpa','plastic':'plastic_strain_compression_positive','pc':'preconsolidation_pressure_kpa','void':'void_ratio'}[field]
        return np.array([s[key] for s in x['states']])
    a=np.array(x['states']).reshape(-1,14)
    if field=='stress':return a[:,:6]
    if field=='plastic':return a[:,8:]
    if field=='pc':return np.exp(a[:,6])
    return np.exp(a[:,7])-1


def compare(a,b):
    ca=read(a/'controls.json');cb=read(b/'controls.json')
    native=ca['binding']['map']=='native'
    native_b=cb['binding']['map']=='native'
    aa={x['coordinate']:x for x in rows(a/'samples.jsonl')};count=0;maximum={};absolute={};end=None;load_difference=0.;reduced_load_difference=0.
    Ra=np.array(ca['transform']);Rb=np.array(cb['transform'])
    def load(c,x):
        h,v=x['fractions']
        return np.array(c['final_load_kn'])+(1-h)*np.array(c['horizontal_holding_kn'])+(1-v)*np.array(c['vertical_holding_kn'])
    for y in rows(b/'samples.jsonl'):
        if y['coordinate'] not in aa:continue
        x=aa[y['coordinate']];count+=1;end=y['coordinate']
        la,lb=load(ca,x),load(cb,y)
        load_difference=max(load_difference,float(np.linalg.norm(la-lb)/max(np.linalg.norm(lb),1.)))
        # Full nodal load includes differing support/reaction conventions.
        # Compare generalized force in each model's compatible free basis too.
        reduced_load_difference=max(reduced_load_difference,float(np.linalg.norm(Ra.T@la-Rb.T@lb)/max(np.linalg.norm(Rb.T@lb),1.)))
        for f in ['displacement','stress','plastic','pc','void']:
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
        for order in ['horizontal','vertical','proportional']:
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
                    n+=1
                result['restart'][kind+'-'+order]=dict(exact_common_full_fields=n,completion_verified=False)
    for order in ['horizontal','vertical','proportional']:
        for step in ['128','256']:
            a=evidence/('native-'+order+step);b=evidence/('independent-'+order+step)
            if a in paths and b in paths:result['cross_map'][order+step]=compare(a,b)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--evidence',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.write_text(json.dumps(measure(a.evidence),indent=2,allow_nan=False)+'\n')
