"""Reassemble every saved trace field and actual linear system, without new solves."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
from restart_audit import audit as restart_audit
from branch_math import component_bounds


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def rows(path):
    with gzip.open(path, 'rt', encoding='utf-8') as stream:
        for line in stream:
            yield json.loads(line)


def load(model, s):
    first, second = (1.-2*s, 1.) if s <= .5 else (0., 2.-2*s)
    return np.asarray(model['final_load_kn'])+first*np.asarray(model['first_group_kn'])+second*np.asarray(model['second_group_kn'])


def audit(folder):
    report, controls, model = [read(folder/(n+'.json')) for n in ['report', 'controls', 'model']]
    B, R, w = [np.asarray(model[n]) for n in ['B', 'transform', 'weights']]
    active = np.asarray(model['active_gp'])
    inactive = np.setdiff1d(np.arange(len(model['initial_states']) if controls['map']=='native' else np.asarray(model['initial_states']).reshape(-1,14).shape[0]), active)
    initial = model['initial_states']
    events = list(rows(folder/'events.jsonl.gz'))
    accepted = {(x['step'],x['attempt']) for x in events if x['accepted']}
    maximum_residual, maximum_strain, maximum_linear = 0., 0., 0.
    minimum_augmented, minimum_singular, minimum_alignment = float('inf'), float('inf'), float('inf')
    maximum_lr_all, linears, fields = 0., 0, 0
    turns = []
    diagnostics = {}
    best_coordinate=read(folder/'anchor-check.json')['coordinate']
    minimum_holding=float(np.linalg.norm(load(model,best_coordinate)-np.asarray(model['final_load_kn'])))
    for p in rows(folder/'samples.jsonl.gz'):
        states = p['states']
        if controls['map']=='native':
            stress = np.array([states[i]['effective_stress_compression_positive_kpa'] for i in active])
            frozen = all(states[i]==initial[i] for i in inactive)
        else:
            full = np.asarray(states).reshape(-1,14)
            stress = full[active,:6]
            frozen = np.array_equal(full[inactive],np.asarray(initial).reshape(-1,14)[inactive])
        if not frozen or not p['removed_material_frozen'] or p['contacts']!=model['contacts']:
            raise ValueError('removed material/contact history changed')
        target = load(model,p['coordinate'])
        component_bounds(model['holding_kn'], target-np.asarray(model['final_load_kn']))
        holding_norm=float(np.linalg.norm(target-np.asarray(model['final_load_kn'])))
        if holding_norm<minimum_holding:minimum_holding,best_coordinate=holding_norm,p['coordinate']
        f = -np.einsum('qai,qa,q->i',B,stress,w)
        rel = float(np.linalg.norm(R.T@target-f)/max(np.linalg.norm(target),1.))
        strain = float(np.abs(np.einsum('qai,i->qa',B,p['reduced_increment_m'])).sum(axis=1).max())
        if rel >= 1e-8 or abs(p['arc_residual_m']) >= 1e-12 or strain > .001:
            raise ValueError('accepted field fails original equilibrium/arc/strain gates')
        maximum_residual, maximum_strain = max(maximum_residual,rel), max(maximum_strain,strain)
        d = p['branch_diagnostic']
        diagnostics[p['step']] = d
        minimum_augmented = min(minimum_augmented,d['augmented_smallest_singular_kn_m'])
        minimum_singular = min(minimum_singular,d['smallest_singular_kn_m'])
        minimum_alignment = min(minimum_alignment,d['weak_left_load_alignment'])
        if d['turning_candidate']:
            turns.append(dict(step=p['step'],coordinate=p['coordinate'],K_sigma_min=d['smallest_singular_kn_m'],
                 augmented_sigma_min=d['augmented_smallest_singular_kn_m'],weak_left_load_alignment=d['weak_left_load_alignment']))
        fields += 1
    tangent_count = 0
    maximum_singular_pair_residual = 0.
    for row in rows(folder/'tangents.jsonl.gz'):
        k = np.asarray(row['matrix'])
        d = row['diagnostic']
        if d != diagnostics[row['step']]:
            raise ValueError('field/tangent diagnostic differs')
        left,right = np.asarray(d['weak_left_mode']),np.asarray(d['weak_right_mode'])
        sigma = d['smallest_singular_kn_m']
        residual = max(np.linalg.norm(k@right-sigma*left),np.linalg.norm(k.T@left-sigma*right))
        maximum_singular_pair_residual = max(maximum_singular_pair_residual,float(residual))
        measured = np.linalg.svd(k,compute_uv=False)[-1]
        if abs(measured-sigma) > 1e-10*max(np.linalg.norm(k),1.):
            raise ValueError('saved smallest singular value differs')
        tangent_count += 1
    if tangent_count != fields:
        raise ValueError('accepted tangent inventory differs')
    accepted_linear_count = 0
    for row in rows(folder/'linear-systems.jsonl.gz'):
        a,b,x = [np.asarray(row[n]) for n in ['matrix','rhs','solution']]
        lr = float(np.linalg.norm(a@x-b)/max(np.linalg.norm(b),1e-30))
        if lr != row['true_relative_residual']:
            raise ValueError('stored actual linear residual differs')
        maximum_lr_all = max(maximum_lr_all,lr)
        if (row['step'],row['attempt']) in accepted:
            if lr > 1e-10:
                raise ValueError('accepted attempt has failed actual linear solve')
            maximum_linear = max(maximum_linear,lr)
            accepted_linear_count += 1
        linears += 1
    if fields!=report['accepted_increments'] or fields!=len(accepted):
        raise ValueError('accepted history inventory mismatch')
    if any(not x.get('rejected_committed_state_unchanged',False) for x in events if not x['accepted']):
        raise ValueError('rejected committed material changed')
    return dict(report=report, fields_reassembled=fields, maximum_force_residual=maximum_residual,
        maximum_increment_strain_l1=maximum_strain, accepted_tangents=tangent_count,
        maximum_singular_pair_residual=maximum_singular_pair_residual, actual_linear_systems=linears, accepted_actual_linear_systems=accepted_linear_count,
        maximum_accepted_true_linear_residual=maximum_linear, maximum_all_attempts_true_linear_residual=maximum_lr_all,
        minimum_K_singular=minimum_singular if fields else None, minimum_augmented_singular=minimum_augmented if fields else None,
        minimum_weak_left_load_alignment=minimum_alignment if fields else None, turning_candidates=turns,
        all_original_component_bounds=True, removed_material_and_contact_frozen=True,
        minimum_remaining_holding_norm_kn=minimum_holding, coordinate_at_minimum_holding=best_coordinate,
        bifurcation_confirmed=False, numerical_acceptance=False)


def audit_reorder(folder):
    report,controls,model=[read(folder/(name+'.json')) for name in ['report','controls','model']]
    kind=controls['binding']['map']
    B,R,w=[np.asarray(model[name]) for name in ['B','transform','weights']]
    active=np.asarray(model['active_gp'])
    initial=model['initial_states']
    total=len(initial) if kind=='native' else np.asarray(initial).reshape(-1,14).shape[0]
    inactive=np.setdiff1d(np.arange(total),active)
    count=0;maximum_force=0.;maximum_strain=0.;last=None
    for p in rows(folder/'samples.jsonl.gz'):
        states=p['states']
        if kind=='native':
            stress=np.array([states[i]['effective_stress_compression_positive_kpa'] for i in active])
            frozen=all(states[i]==initial[i] for i in inactive)
        else:
            full=np.asarray(states).reshape(-1,14)
            stress=full[active,:6]
            frozen=np.array_equal(full[inactive],np.asarray(initial).reshape(-1,14)[inactive])
        if not frozen or p['contacts']!=model['contacts']:raise ValueError('reorder frozen history changed')
        target=load(model,p['coordinate'])
        component_bounds(model['holding_kn'],target-np.asarray(model['final_load_kn']))
        force=-np.einsum('qai,qa,q->i',B,stress,w)
        residual=float(np.linalg.norm(R.T@target-force)/max(np.linalg.norm(target),1.))
        strain=float(np.abs(np.einsum('qai,i->qa',B,p.get('reduced_increment_m',np.zeros(R.shape[1])))).sum(axis=1).max())
        if residual>=1e-8 or strain>.001:raise ValueError('reorder accepted field fails original gates')
        maximum_force=max(maximum_force,residual);maximum_strain=max(maximum_strain,strain)
        count+=1;last=p
    rejected=0
    for attempt in rows(folder/'attempts.jsonl.gz'):
        if not attempt['accepted']:
            if not attempt.get('rejected_committed_state_unchanged',False):raise ValueError('reorder rejected state changed')
            rejected+=1
    if count!=report['accepted_increments']+1:raise ValueError('reorder accepted inventory differs')
    return dict(report=report,fields_reassembled=count,maximum_force_residual=maximum_force,
        maximum_increment_strain_l1=maximum_strain,rejected_attempts=rejected,
        remaining_holding_norm_kn=float(np.linalg.norm(load(model,last['coordinate'])-np.asarray(model['final_load_kn']))),
        original_component_bounds_verified=True,removed_material_and_contact_frozen=True,
        actual_matrix_scope='original last-128 matrix codec; all solve residual summaries retained',numerical_acceptance=False)


def run(source):
    records = {}
    for folder in sorted(source.iterdir()):
        if folder.is_dir() and (folder/'model.json').exists() and (folder/'report.json').exists():
            records[folder.name] = audit(folder) if (folder/'events.jsonl.gz').exists() else audit_reorder(folder)
    if not records:
        raise ValueError('no completed diagnostic traces')
    return dict(passed=False, full_release_acceptance=False, static_limit_confirmed=False,
                independent_completion_verified=False, restart_to_completion_verified=False, partial_restart=restart_audit(source), traces=records)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    a.output.write_text(json.dumps(run(a.source),indent=2,allow_nan=False)+'\n',encoding='utf-8')
