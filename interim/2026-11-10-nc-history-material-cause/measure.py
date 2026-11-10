"""Reassemble saved states and diagnose MCC stress paths; never certify release."""
import argparse
import copy
from collections import Counter
import gzip
import json
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'tools'))
import diagnose_nc_material_path as material
import fem_reference as fem
from joint_reference import System as Joint
from probe_low_preload_nc import force_check, deterministic
from probe_nc_midpoint import inputs
from verify_high_order import prepare, internal
from low_preload_nc_reference import ClosedStickSystem
import mcc_rate_fem as rate
from diagnose_low_preload_nc import matrix


def read(p):
    with (gzip.open(p, 'rt', encoding='utf-8') if p.suffix == '.gz' else p.open(encoding='utf-8')) as stream:
        return json.load(stream)


def rows(p):
    with (gzip.open(p, 'rt', encoding='utf-8') if p.suffix == '.gz' else p.open(encoding='utf-8')) as stream:
        return [json.loads(line) for line in stream]


def parameters(mesh, ground, counts):
    materials = {m['id']: m['mcc'] for m in ground['materials']}
    return {k: np.repeat([materials[i][k] for i in mesh['material_ids']], counts)
            for k in ('critical_state_m', 'lambda', 'kappa', 'poisson')}


def constitutive(states, kind, params, pi):
    t = material.continuum_terms(states, kind, params)
    eig = np.linalg.eigvalsh(t['plastic_loading'][pi])[:, 0]
    elastic = np.linalg.eigvalsh(t['elastic'][pi])[:, 0]
    ndn = t['normal_elastic_normal'][pi]
    return dict(hypothetical_plastic_loading_negative_eigenvalue_points=int(np.sum(eig < 0)),
                plastic_loading_minimum_eigenvalue_kpa=material.metrics(eig),
                elastic_minimum_eigenvalue_kpa=material.metrics(elastic),
                consistency_denominator_over_elastic_term=material.metrics(t['denominator'][pi]/ndn),
                hardening_term=material.metrics(t['hardening'][pi]),
                scope='existing continuum equations at saved stress; no discrete tangent or branch/stability certificate')


def global_continuum(mesh, ground, case, run):
    state = run['last_committed']; model = ClosedStickSystem(mesh, ground, state['active_elements'], case)
    stress, pc, void, plastic_strain = material.unpack(state['material_states'], 'native')
    all_state = np.column_stack((stress, np.log(pc), np.log1p(void), plastic_strain))
    active_gp = np.repeat(state['active_elements'], model.b.shape[1])
    observed = all_state[active_gp]
    yf, elastic, normal, dn, den = rate.terms(observed, model.parameters)
    holding = np.asarray(run['excavation_audits'][0]['initial_holding_load_kn'])
    selected = model.select(observed, -holding, yf >= -2e-9)
    tangent = elastic.copy()
    tangent[selected] -= np.einsum('qa,qb,q->qab', dn[selected], dn[selected], 1/den[selected])
    k, ke = matrix(model, tangent), matrix(model, elastic)
    left, singular, vt = np.linalg.svd(k)
    mode = vt[-1]
    coordinates = np.zeros(model.xyz.size); coordinates[model.free] = mode
    displacement = model.lift(coordinates)
    strain = -np.einsum('eqai,ei->eqa', model.b, displacement[model.dofs]).reshape(-1, 6)
    contribution = model.weight.ravel()*np.einsum('qa,qab,qb->q', strain, tangent, strain)
    direct = material.continuum_terms(state['material_states'], 'native', parameters(mesh, ground, np.diff(prepare(mesh)[1])))
    reconstruction = elastic.copy()
    reconstruction[selected] = direct['plastic_loading'][active_gp][selected]
    return dict(active_contacts=int(model.live.sum()), selected_plastic_loading_points=int(selected.sum()),
        condition_2=float(singular[0]/singular[-1]), singular_values_kn_per_m=singular.tolist(),
        elastic_singular_values_kn_per_m=np.linalg.svd(ke, compute_uv=False).tolist(),
        asymmetry_relative=float(np.linalg.norm(k-k.T)/np.linalg.norm(k)),
        right_mode_gp_contribution_kn_per_m=contribution.tolist(),
        right_mode_contribution_sum_kn_per_m=float(contribution.sum()),
        right_mode_contribution_absolute_sum_kn_per_m=float(np.abs(contribution).sum()),
        right_mode_elastic_energy_coefficient_kn_per_m=float(mode@ke@mode),
        negative_gp_contributions=int(np.sum(contribution < 0)),
        mode_reassembly_absolute_error_kn_per_m=abs(float(contribution.sum()-mode@k@mode)),
        direct_equation_vs_existing_rate_tangent_relative=float(np.linalg.norm(reconstruction-tangent)/np.linalg.norm(tangent)),
        scope='same saved field and release direction, existing continuum branch selection; not actual discrete Newton matrix, independent solution, or physical failure verdict')


def fresh(folder, n):
    ip = folder/f'input{n}'; path = folder/f'trace{n}'
    mesh, ground, case = [read(ip/(name+'.json')) for name in ('mesh', 'ground', 'stages')]
    expected = inputs('backward_euler'); expected[2]['stages'][2]['initial_increment'] = 1/n
    if [mesh, ground, case] != list(expected):
        raise ValueError('not the unchanged physical 5 kPa NC input with declared controls')
    run = read(path/'run.json')
    raw = path/'accepted.jsonl'
    samples = rows(raw if raw.exists() else path/'accepted.jsonl.gz')
    q, offset = prepare(mesh); counts = np.diff(offset)
    weights = np.concatenate([p['weight'] for p in q]); params = parameters(mesh, ground, counts)
    joint = Joint(mesh, case['joint_system'], run['joint_integration_points'])
    previous = run['geostatic']['seed']; prevlive = np.ones(len(joint.points), dtype=bool)
    ledger = []; transitions = Counter(); frozen_steps = 0; maximum_force = 0.
    for index, sample in enumerate(samples, 1):
        if sample['accepted_count'] != index:
            raise ValueError('missing or repeated accepted state')
        state = sample['state']; active = np.asarray(state['active_elements'])
        pi = np.flatnonzero(np.repeat(active, counts)); dead = np.flatnonzero(~np.repeat(active, counts))
        live = np.array([active[p['owner_a']] and active[p['owner_b']] for p in joint.points])
        points = state['joints']['points']; oldpoints = previous['joints']['points']
        for i, (p, old) in enumerate(zip(points, oldpoints)):
            if live[i] and prevlive[i]:
                transitions[f"{old['status']}->{p['status']}"] += 1
            if not live[i] and p != old:
                raise ValueError('inactive contact history changed')
        if any(state['material_states'][i] != previous['material_states'][i] for i in dead):
            raise ValueError('inactive material history changed')
        if len(dead): frozen_steps += 1
        force = force_check(mesh, ground, case, state, q, offset, joint)
        maximum_force = max(maximum_force, force['relative_force_residual'])
        desc = material.describe(state['material_states'], 'native', params['critical_state_m'],
                                 pi, weights[pi], previous['material_states'])
        gap = np.array([np.dot(p['jump_m'], p['normal']) for p in points])
        ledger.append(dict(accepted_count=index, geostatic_complete=sample['geostatic_complete'],
            stage_index=state['stage_index'], progress=state['stage_progress'],
            active_contacts=int(live.sum()), active_contact_statuses=dict(Counter(p['status'] for i,p in enumerate(points) if live[i])),
            stored_contact_statuses=dict(Counter(p['status'] for p in points)),
            active_gap_m=material.metrics(gap[live]) if live.any() else None,
            stored_traction_l2_kpa=float(np.linalg.norm([p['traction_kpa'] for p in points])),
            current_contact_nodal_force_norm_kn=float(np.linalg.norm(joint.audit(state))),
            force=force, material=desc))
        previous = state; prevlive = live
    if samples[-1]['state'] != run['last_committed']:
        raise ValueError('trace final state differs')
    attempts = run['geostatic']['attempts']+run['attempts']
    if any(not a['accepted'] and not a['committed_state_unchanged_on_rejection'] for a in attempts):
        raise ValueError('rejection changed committed state')
    accepted = sum(a['accepted'] for a in attempts)
    if accepted != len(samples): raise ValueError('accepted state inventory differs')
    terminal = run['last_committed']; pi = np.flatnonzero(np.repeat(terminal['active_elements'], counts))
    # Deliberate read-only counterfactual on an unsealed COPY, never a solver
    # state or new load path. Show which stored histories enter reassembly.
    dead = np.flatnonzero(~np.repeat(terminal['active_elements'], counts))
    cleared = copy.deepcopy(terminal)
    for i in dead:
        cleared['material_states'][i]['effective_stress_compression_positive_kpa'] = [0.]*6
    for p in cleared['joints']['points']:
        for key in ('jump_m', 'traction_kpa', 'plastic_slip_m'): p[key] = [0.]*3
        for key in ('accumulated_plastic_slip_m','plastic_dissipation_kpa_m','numerical_dissipation_kpa_m','endpoint_work_kpa_m'): p[key] = 0.
        p['status'] = 'open'
    def assembled(s):
        return internal(mesh,ground,s['material_states'],s['active_elements'],q,offset)[0]+joint.audit(s)
    force_original, force_cleared = assembled(terminal), assembled(cleared)
    if not np.array_equal(force_original, force_cleared):
        raise ValueError('inactive-history clearing changed reconstructed current force')
    # This altered active stress is not a constitutively admissible solve;
    # evaluate only its nodal force analytically, not internal()'s yield audit.
    ei = next(e for e in range(len(q)) if offset[e] <= pi[0] < offset[e+1])
    qi = int(pi[0]-offset[ei]); response = -q[ei]['b'][qi,2]*q[ei]['weight'][qi]
    return dict(passed=False, accepted_states=len(samples), ledger=ledger,
        accepted_contact_transitions_including_seed=dict(transitions), inactive_frozen_steps=frozen_steps,
        maximum_reassembled_force_residual=maximum_force, terminal_progress=terminal['stage_progress'],
        complete=run['complete'], failure_reason=run['failure_reason'], last_attempts=run['attempts'][-5:],
        excavation_audits=run['excavation_audits'], terminal_constitutive=constitutive(terminal['material_states'], 'native', params, pi),
        terminal_global_continuum=global_continuum(mesh, ground, case, run),
        rejected_committed_states_unchanged=True,
        inactive_history_read_only_intervention=dict(removed_material_points=len(dead),
            unsealed_copy_only=True, checkpoint_edited_or_resealed=False,
            all_current_nodal_forces_bitwise_unchanged=True,
            active_one_kpa_gp_stress_probe_force_norm_kn=float(np.linalg.norm(response)),
            scope='Python reassembly dependency test; no Rust trial, independent solution, or altered physical history'),
        contact_scope='seed and every accepted increment; no assertion about every internal Newton contact candidate')


def prior(folder, mesh, ground):
    q, offset = prepare(mesh); params = parameters(mesh, ground, np.diff(offset))
    results = {}
    for p in sorted((folder/'nc/paths').iterdir()):
        model = read(p/'model.json'); report = read(p/'report.json'); samples = rows(p/'samples.jsonl.gz')
        kind = model['map']; pi = np.asarray(model['active_gp']); weights = np.asarray(model['weights'])
        if weights.shape != pi.shape or report['passed']:
            raise ValueError('prior path inventory or acceptance mismatch')
        B, R = np.asarray(model['B']), np.asarray(model['transform'])
        previous = model['initial_states']; timeline = []; maximum_force = 0.
        dead = np.setdiff1d(np.arange(len(previous)), pi)
        for sample in samples:
            f = material.unpack(sample['states'], kind)
            if kind == 'native':
                frozen = all(sample['states'][i] == model['initial_states'][i] for i in dead)
            else:
                frozen = np.array_equal(np.asarray(sample['states'])[dead], np.asarray(model['initial_states'])[dead])
            if not frozen or sample['contacts'] != model['contacts']:
                raise ValueError('prior inactive history changed')
            a,b = sample['fractions']
            load = np.asarray(model['final_load_kn'])+(1-a)*np.asarray(model['first_group_kn'])+(1-b)*np.asarray(model['second_group_kn'])
            force = -np.einsum('qai,qa,q->i', B, f[0][pi], weights)
            residual = float(np.linalg.norm(R.T@load-force)/max(np.linalg.norm(load),1.))
            if residual >= 1e-8: raise ValueError('prior saved force equilibrium failed')
            maximum_force = max(maximum_force, residual)
            timeline.append(dict(coordinate=sample['coordinate'], material=material.describe(sample['states'], kind,
                params['critical_state_m'], pi, weights, previous)))
            previous = sample['states']
        attempts = rows(p/'attempts.jsonl.gz')
        zero = [a for a in attempts if a.get('label') == 'lambda-zero']
        results[p.name] = dict(passed=False, report=report, saved_states=len(samples), timeline=timeline,
            maximum_reassembled_force_residual=maximum_force, removed_history_frozen=True,
            lambda_zero_actual_newton=zero, terminal_constitutive=constitutive(samples[-1]['states'], kind, params, pi))
    return results


def measure(data):
    observed = {str(n): fresh(data, n) for n in (128,256)}
    mesh, ground = [read(data/'input128'/(n+'.json')) for n in ('mesh','ground')]
    histories = []
    for n in (128,256):
        p = data/f'trace{n}'/'accepted.jsonl'
        samples = rows(p if p.exists() else p.with_suffix('.jsonl.gz'))
        histories.append({(r['geostatic_complete'], r['state']['stage_index'], r['state']['stage_progress']): r['state'] for r in samples})
    common = sorted(histories[0].keys() & histories[1].keys())
    differences = {}
    for key, column in [('displacement_m', None), ('stress', 0), ('pc', 1), ('void_ratio', 2), ('plastic_strain', 3)]:
        pairs = []
        for t in common:
            a,b = [h[t] for h in histories]
            x,y = ([np.asarray(s['displacement_m']) for s in (a,b)] if column is None else
                   [material.unpack(s['material_states'], 'native')[column] for s in (a,b)])
            pairs.append((float(np.max(np.abs(x-y))), float(np.linalg.norm(x-y)/max(np.linalg.norm(y),1e-30))))
        differences[key] = dict(maximum_absolute=max(p[0] for p in pairs), maximum_relative_l2=max(p[1] for p in pairs))
    return dict(schema='soil-fem-nc-history-material-cause/1', passed=False,
        starting_commit='744fce031c551fe2e1eee9dd69e02f39c2785a9c',
        physical_input=dict(pc_kpa=100., additional_preload_kpa=5., force_kn=[1.,.5,0.]),
        fresh_native=observed, archived_oblique=prior(data/'prior', mesh, ground),
        fresh_native_step_sensitivity=dict(common_accepted_coordinates=len(common), fields=differences, complete_history_converged=False),
        full_release_observed_in_both_native_runs=all(v['complete'] for v in observed.values()),
        full_release_acceptance=False, full_release_freeze=False, independent_plastic_full_path=False,
        both_implementations_completion_step_sensitivity=False, completion_restart=False,
        physical_failure_confirmed=False, p3b_p4_complete=False, p6_complete=False,
        note='保存場の再組立と材料応答の原因診断。新しい独立求解・破壊・全解放認定ではない。')


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--data',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    with args.output.open('x',encoding='utf-8') as out:
        json.dump(measure(args.data.resolve()),out,indent=2,ensure_ascii=False,allow_nan=False);out.write('\n')
