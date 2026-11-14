"""Reassemble forces, work and spectra from saved fields, without solving MCC."""
import argparse
import hashlib
import json
import itertools
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'2026-11-09-nc-beyond-fold'))
sys.path.insert(0, str(HERE.parent/'2026-11-09-nc-pseudo-arclength'))
from analyze import native_geometry, native_force, independent_geometry, diff, accepted_count
from continue_path import digest, read, write


def stream(path):
    if path.exists():
        with path.open() as f:
            for line in f:
                yield json.loads(line)


def fields(row, native):
    s = row['states']
    if native:
        return dict(stress=[p['effective_stress_compression_positive_kpa'] for p in s],
                    plastic_strain=[p['plastic_strain_compression_positive'] for p in s],
                    pc=[p['preconsolidation_pressure_kpa'] for p in s], void_ratio=[p['void_ratio'] for p in s])
    s = np.array(s).reshape(-1, 14)
    return dict(stress=s[:, :6], plastic_strain=s[:, 8:], pc=np.exp(s[:, 6]), void_ratio=np.exp(s[:, 7])-1)


def compare(a, b, native):
    fa, fb = fields(a, native), fields(b, native)
    return dict(progress_absolute_difference=abs(a['progress']-b['progress']),
                displacement=diff(a['displacement_m'], b['displacement_m']),
                material={key: diff(fa[key], fb[key]) for key in fa})


def measure(evidence, out):
    out.mkdir(parents=True, exist_ok=False)
    g = native_geometry(evidence/'origin')
    if accepted_count(g.run) != 266 or g.lam0 != .9514000415802002:
        raise ValueError('production history differs from authorized base')
    s, active, before, delta = independent_geometry(evidence/'origin/reference')
    assert not s.live.any() and not np.any(s.kjoint)
    R = np.zeros((s.xyz.size, len(s.free)))
    for j, dof in enumerate(s.free):
        v = np.zeros(s.xyz.size); v[dof] = 1.; R[:, j] = s.lift(v)
    Bi = np.einsum('eqai,eij->eqaj', s.b, R[s.dofs]).reshape(-1, 6, len(s.free))
    pinv = np.linalg.pinv(g.R)
    paths, endpoints = {}, {}
    for folder in sorted(evidence.glob('extended-*')):
        if not folder.is_dir() or not (folder/'report.json').exists():
            continue
        native = read(folder/'controls.json')['binding']['map'] == 'native'
        initial = read(folder/'initial.json')['payload']
        expected_contact_count = len(initial['contacts']['points' if native else 'jump_m'])
        if expected_contact_count != 32:
            raise ValueError('contact inventory differs')
        frozen_initial = np.array(initial['states'], dtype=object) if native else np.array(initial['states'])
        previous = initial
        rr, audits, works, progress, increments = [], [], [], [], []
        count = 0
        minimum_active_mean_pressure = float('inf')
        changed_plastic_gp = np.zeros(162, dtype=bool)
        removed_frozen = contact_frozen = True
        endpoint = initial
        for row in stream(folder/'samples.jsonl'):
            if np.array(row['displacement_m']).shape != (54, 3) or np.array(fields(row, native)['stress']).shape != (324, 6):
                raise ValueError('full-field inventory differs')
            if not np.isfinite(np.array(row['displacement_m'])).all() or not all(np.isfinite(v).all() for v in fields(row, native).values()):
                raise ValueError('nonfinite saved field')
            lam = row['progress']
            if native:
                target = g.load0-(lam-g.lam0)*g.holding
                residual = g.R.T@target.ravel()-native_force(g, row['states'])
                stress = np.array([row['states'][i]['effective_stress_compression_positive_kpa'] for i in g.pi])
                previous_stress = np.array([previous['states'][i]['effective_stress_compression_positive_kpa'] for i in g.pi])
                dep = np.array([np.array(row['states'][i]['plastic_strain_compression_positive'])-previous['states'][i]['plastic_strain_compression_positive'] for i in g.pi])
                du = pinv@(np.array(row['displacement_m'])-previous['displacement_m']).ravel()
                de = -np.einsum('qai,i->qa', g.B, du)
                weight = g.weight
                removed_frozen &= all(row['states'][i] == initial['states'][i] for i in range(len(row['states'])) if i not in g.pi)
                target0 = g.load0-(previous['progress']-g.lam0)*g.holding
            else:
                full = np.array(row['states']); old = np.array(previous['states'])
                target = before+lam*delta
                residual = s.reduce(target-s.total_force(full[active].reshape(-1, 14), np.array(row['displacement_m']).ravel()))
                stress = full[active].reshape(-1, 14)[:, :6]
                previous_stress = old[active].reshape(-1, 14)[:, :6]
                dep = (full-old)[active].reshape(-1, 14)[:, 8:]
                du = (np.array(row['displacement_m'])-previous['displacement_m']).ravel()[s.free]
                de = -np.einsum('qai,i->qa', Bi, du)
                weight = s.weight.ravel()
                removed_frozen &= np.array_equal(full[~active], frozen_initial[~active])
                target0 = before+previous['progress']*delta
            minimum_active_mean_pressure = min(minimum_active_mean_pressure, float(stress[:, :3].mean(axis=1).min()))
            changed_plastic_gp |= np.any(dep != 0., axis=1)
            contact_frozen &= row['contacts'] == initial['contacts']
            rr.append(float(np.linalg.norm(residual)/max(np.linalg.norm(target), 1.)))
            avg = (stress+previous_stress)/2
            iw = float(np.sum(np.einsum('qa,qa->q', avg, de)*weight))
            pw = float(np.sum(np.einsum('qa,qa->q', avg, dep)*weight))
            ew = float(((np.array(target0)+target)/2).ravel()@(np.array(row['displacement_m'])-previous['displacement_m']).ravel())
            works.append((iw, ew, pw))
            audits.append(row['audit']); progress.append(lam); increments.append(row['progress_increment'])
            count += 1; endpoint = row; previous = row
        endpoints[folder.name] = endpoint
        maxlinear = 0.; linearcount = 0; above = 0; spectra = []; eigdifference = 0.
        for eq in stream(folder/'linear-systems.jsonl'):
            A, b, x = map(np.array, (eq['matrix'], eq['rhs'], eq['solution']))
            lr = float(np.linalg.norm(A@x-b)/max(np.linalg.norm(b), 1e-30))
            maxlinear = max(maxlinear, lr); linearcount += 1; above += lr > 1e-10
            if eq['kind'] == 'tangent':
                K = A[:-1, :-1]
                spectra.append(float(np.linalg.eigvalsh((K+K.T)/2)[0]))
        for eig, row in zip(spectra, audits):
            eigdifference = max(eigdifference, abs(eig-row['symmetric_eigen_min_kn_m']))
        report = read(folder/'report.json')
        assert count == report['equilibrium_samples'] and len(spectra) == count
        predictor = None
        if (folder/'seed-tangent-system.json').exists():
            eq = read(folder/'seed-tangent-system.json')
            A, b, x = map(np.array, (eq['matrix'], eq['rhs'], eq['solution']))
            predictor = dict(recomputed_true_relative_residual=float(np.linalg.norm(A@x-b)/max(np.linalg.norm(b), 1e-30)),
                             stability_evidence=False)
        paths[folder.name] = dict(report, recounted_samples=count,
            diagnostic_verdict='inconclusive', zero_increment_predictor_system=predictor,
            active_material_gp=162, removed_material_gp=162, frozen_contact_gp=32,
            active_gp_with_changed_plastic_strain=int(changed_plastic_gp.sum()),
            minimum_active_mean_pressure_kpa=minimum_active_mean_pressure,
            maximum_sampled_progress=max([initial['progress']]+progress), minimum_sampled_progress=min([initial['progress']]+progress),
            positive_progress_increments=sum(x > 0 for x in increments), negative_progress_increments=sum(x < 0 for x in increments),
            maximum_reassembled_relative_force_residual=max(rr, default=None), removed_material_frozen=bool(removed_frozen),
            contact_history_frozen=bool(contact_frozen), active_contact_points=0,
            saved_linear_systems=linearcount, maximum_recomputed_true_linear_residual=maxlinear, linear_systems_above_1e_10=int(above),
            minimum_symmetric_eigenvalue_kn_m=min(spectra, default=None), maximum_symmetric_eigenvalue_min_kn_m=max(spectra, default=None),
            positive_to_negative_eigen_crossings=sum(a > 0 and b < 0 for a, b in zip(spectra, spectra[1:])),
            maximum_eigen_audit_difference=eigdifference,
            negative_directional_work_samples=sum(a['directional_incremental_work_kn_m'] < 0 for a in audits),
            cumulative_internal_work_kn_m=sum(x[0] for x in works), cumulative_external_work_kn_m=sum(x[1] for x in works),
            cumulative_plastic_stress_work_kn_m=sum(x[2] for x in works),
            maximum_absolute_work_imbalance_kn_m=max((abs(x[0]-x[1]) for x in works), default=None),
            saved_fields_chain_sha256=digest([digest(row) for row in stream(folder/'samples.jsonl')]),
            energy_interpretation='stress work only; no specified hardening free energy; no physical stability verdict')
    # No constitutive or equilibrium solve: only stored-field reassembly.
    sensitivity = {}
    restarts = {}
    captures = read(evidence/'checkpoint-captures.json')
    for kind in ['native', 'independent']:
        coarse = stream(evidence/f'extended-{kind}-next/samples.jsonl')
        fine = itertools.chain(stream(evidence/f'extended-{kind}-half-a/samples.jsonl'),
                               stream(evidence/f'extended-{kind}-half-b/samples.jsonl'))
        matched = 0
        maxima = {}
        for a, b in itertools.zip_longest(coarse, itertools.islice(fine, 1, None, 2)):
            if a is None or b is None:
                raise ValueError('sensitivity inventories differ')
            c = compare(a, b, kind == 'native')
            flat = dict(progress_absolute_difference=c['progress_absolute_difference'])
            for key, value in c['displacement'].items(): flat['displacement_'+key] = value
            for field, values in c['material'].items():
                for key, value in values.items(): flat[field+'_'+key] = value
            for key, value in flat.items(): maxima[key] = max(maxima.get(key, 0.), value)
            matched += 1
        if matched != 4096:
            raise ValueError('full additional interval comparison requires 4096 points')
        parent = read(evidence/f'parents/{kind}.json')['payload']
        coarse_start = read(evidence/f'extended-{kind}-next/initial.json')['payload']
        fine_start = read(evidence/f'extended-{kind}-half-a/initial.json')['payload']
        same_fields = all(parent[key] == coarse_start[key] == fine_start[key]
                          for key in ['states', 'displacement_m', 'contacts', 'progress', 'tangent'])
        sensitivity[kind] = dict(matched_coarse_samples=matched, coarse_samples=4096, fine_samples=8192,
            same_parent_full_fields=same_fields, parent_arc_length_m=parent['arc_length_m'],
            additional_interval_m=.4096, coarse_step_m=.0001, fine_step_m=.00005,
            path_maximum_differences=maxima,
            endpoint=compare(endpoints[f'extended-{kind}-next'], endpoints[f'extended-{kind}-half-b'], kind == 'native'),
            full_history_converged=False, scope='additional 0.9216 to 1.3312 m interval only; inherited earlier errors retained')
        capture, = [c for c in captures if c['map'] == kind]
        cp_path = evidence/f'captured-{kind}.json'
        if hashlib.sha256(cp_path.read_bytes()).hexdigest() != capture['sha256']:
            raise ValueError('captured checkpoint bytes differ')
        cp = read(cp_path)
        if cp['payload_sha256'] != digest(cp['payload']) or not capture['actual_runtime_checkpoint']:
            raise ValueError('runtime checkpoint payload differs')
        resumed = list(stream(evidence/f'extended-{kind}-restart128/samples.jsonl'))
        whole = [r for r in stream(evidence/f'extended-{kind}-next/samples.jsonl')
                 if capture['step'] < r['step'] <= capture['step']+128]
        exact = resumed == whole and len(whole) == 128
        restarts[kind] = dict(capture=capture, whole_rows=len(whole), restarted_rows=len(resumed),
            all_saved_fields_and_audits_exact=exact, production_job_restart=False,
            restart_to_completion_verified=False)
    origin = g.run['continuation']['origin']
    zero_load = g.load0+g.lam0*g.holding
    native_zero = float(np.linalg.norm(g.R.T@zero_load.ravel()-native_force(g, origin['material_states']))/max(np.linalg.norm(zero_load), 1.))
    unload = read(evidence/'origin/reference/unload.json')
    izero = s.total_force(np.array(unload['material_state'])[active].reshape(-1, 14), np.array(unload['displacement_m']).ravel())
    independent_zero = float(np.linalg.norm(s.reduce(izero-before))/max(np.linalg.norm(before), 1.))
    zero_solutions = {}
    for kind in ['native', 'independent']:
        folder = evidence/f'zero-{kind}'
        report = read(folder/'report.json')
        initial = read(folder/'initial.json')
        events = list(stream(folder/'events.jsonl'))
        systems = list(stream(folder/'linear-systems.jsonl'))
        linear = [float(np.linalg.norm(np.array(eq['matrix'])@eq['solution']-eq['rhs'])/max(np.linalg.norm(eq['rhs']), 1e-30)) for eq in systems]
        summary = dict(report, initial_increment_norm_m=float(np.linalg.norm(events[0]['reduced_increment_m'])),
            initial_relative_force_residual=events[0]['relative_force_residual'], newton_evaluations=len(events),
            actual_linear_systems=len(systems), maximum_recomputed_true_linear_residual=max(linear),
            holding_target='original pre-cut-derived holding loads; unchanged',
            force_reassembly=None, removed_material_frozen=None, contact_history_frozen=None)
        if report['diagnostic_equilibrium_solved']:
            sol = read(folder/'solution.json')
            if kind == 'native':
                target = zero_load
                residual = g.R.T@target.ravel()-native_force(g, sol['states'])
                frozen = all(sol['states'][i] == initial['states'][i] for i in range(324) if i not in g.pi)
            else:
                full = np.array(sol['states'])
                target = before
                residual = s.reduce(target-s.total_force(full[active].reshape(-1, 14), np.array(sol['displacement_m']).ravel()))
                frozen = np.array_equal(full[~active], np.array(initial['states'])[~active])
            summary.update(force_reassembly=float(np.linalg.norm(residual)/max(np.linalg.norm(target), 1.)),
                removed_material_frozen=bool(frozen), contact_history_frozen=sol['contacts'] == initial['contacts'],
                returned_increment_norm_m=float(np.linalg.norm(sol['reduced_increment_m'])),
                returned_full_field_difference=compare(dict(initial, progress=0.), sol, kind == 'native'))
        zero_solutions[kind] = summary
    cross = dict(
        progress_absolute_difference=abs(endpoints['extended-native-next']['progress']-endpoints['extended-independent-next']['progress']),
        displacement=diff(endpoints['extended-native-next']['displacement_m'], endpoints['extended-independent-next']['displacement_m']),
        material={key: diff(fields(endpoints['extended-native-next'], True)[key], fields(endpoints['extended-independent-next'], False)[key])
                  for key in fields(endpoints['extended-native-next'], True)},
        comparison_coordinate='equal diagnostic arc length from separate independently obtained seeds; not equal lambda or converged history')
    result = dict(schema='soil-fem-nc-continuation-ci/1', base_commit='42473341dcfaa3c6c969a28f66fe0f5e0d312037',
        passed=False, history_converged=False, static_limit_confirmed=False, verdict='inconclusive', production_adopted=False,
        native_highest_accepted_progress=g.lam0, native_accepted_increments=accepted_count(g.run), native_new_accepted_increments=0,
        full_removal_reached=False, full_release_freeze_verified=False, restart_to_completion_verified=False,
        condition=dict(nodes=54, tet10=12, material_gp=324, contact_gp=32, preload_kpa=5, pc_kpa=100, force_kn=[1, .5, 0]),
        remaining_region_lambda_zero_native=dict(relative_force_residual=native_zero, source='original pre-cut field and original holding loads; reassembled, not re-solved'),
        remaining_region_lambda_zero_independent=dict(relative_force_residual=independent_zero,
            note='input-only pre-cut internal-force/holding-load identity; not a separately solved equilibrium or full release'),
        remaining_region_lambda_zero_perturbed_solves=zero_solutions,
        paths=paths, step_sensitivity=sensitivity, restarts=restarts, both_map_endpoint_comparison=cross,
        missing=['lambda=1 full removal', 'full release freeze', 'whole monotonic removal history convergence',
                 'independent full-removal equilibrium', 'production restart to completion', 'global branch exclusion',
                 'physical stability / specified hardening free energy'])
    write(out/'measurements.json', result)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--evidence', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    measure(a.evidence.resolve(), a.output.resolve())
