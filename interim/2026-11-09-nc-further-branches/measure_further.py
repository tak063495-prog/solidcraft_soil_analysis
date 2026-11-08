"""Reassemble forces, work and spectra from saved fields, without solving MCC."""
import argparse
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
        frozen_initial = np.array(initial['states'], dtype=object) if native else np.array(initial['states'])
        previous = initial
        rr, audits, works, progress, increments = [], [], [], [], []
        count = 0
        removed_frozen = contact_frozen = True
        endpoint = initial
        for row in stream(folder/'samples.jsonl'):
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
    # No constitutive evaluation here: independent force reassembly uses stored stresses.
    roots = {}
    for folder in sorted(evidence.glob('roots-*')):
        if not (folder/'report.json').exists():
            continue
        controls = read(folder/'controls.json')
        native = controls['binding']['parent_binding']['map'] == 'native'
        seed = read(folder/'seed.json')['payload']
        rows = read(folder/'solutions.json')
        residuals = []
        frozen = contacts = True
        for row in rows:
            lam = row['progress']
            if native:
                target = g.load0-(lam-g.lam0)*g.holding
                residual = g.R.T@target.ravel()-native_force(g, row['states'])
                frozen &= all(row['states'][i] == seed['states'][i] for i in range(len(row['states'])) if i not in g.pi)
            else:
                full = np.array(row['states'])
                target = before+lam*delta
                residual = s.reduce(target-s.total_force(full[active].reshape(-1, 14), np.array(row['displacement_m']).ravel()))
                frozen &= np.array_equal(full[~active], np.array(seed['states'])[~active])
            contacts &= row['contacts'] == seed['contacts']
            residuals.append(float(np.linalg.norm(residual)/max(np.linalg.norm(target), 1.)))
        maxlinear = 0.
        above = systems = 0
        for eq in stream(folder/'linear-systems.jsonl'):
            A, b, x = map(np.array, (eq['matrix'], eq['rhs'], eq['solution']))
            lr = float(np.linalg.norm(A@x-b)/max(np.linalg.norm(b), 1e-30))
            maxlinear = max(maxlinear, lr)
            above += lr > 1e-10
            systems += 1
        clusters = read(folder/'distinct-roots.json')
        roots[folder.name] = dict(read(folder/'report.json'), reassembled_solutions=len(rows),
            maximum_reassembled_relative_force_residual=max(residuals, default=None),
            maximum_recomputed_true_linear_residual=maxlinear, linear_systems_above_1e_10=int(above), saved_linear_systems=systems,
            removed_material_frozen=bool(frozen), contact_history_frozen=bool(contacts),
            cluster_increment_norms_m=[float(np.linalg.norm(x['reduced_increment_m'])) for x in clusters],
            note='No nonzero candidate root found inside this finite inventory; no global branch verdict')
    sensitivity = {}
    restarts = {}
    for kind in ['native', 'independent']:
        coarse = list(itertools.islice(stream(evidence/f'extended-{kind}-long/samples.jsonl'), 512))
        for label in ['localhalf', 'fine1024']:
            fine = list(stream(evidence/f'extended-{kind}-{label}/samples.jsonl'))
            matched = list(zip(coarse, fine[1::2]))
            maxima = {}
            for a, b in matched:
                result = compare(a, b, kind == 'native')
                items = dict(progress_absolute_difference=result['progress_absolute_difference'],
                             displacement_relative_l2=result['displacement']['relative_l2'],
                             **{key+'_relative_l2':v['relative_l2'] for key, v in result['material'].items()})
                for key, value in items.items():
                    if key not in maxima or value > maxima[key]['value']:
                        maxima[key] = dict(value=value, coarse_step=a['step'], progress=a['progress'])
            ca = read(evidence/f'extended-{kind}-long/initial.json')['payload']
            fa = read(evidence/f'extended-{kind}-{label}/initial.json')['payload']
            same_start = all(ca[key] == fa[key] for key in ['states', 'contacts', 'displacement_m', 'progress', 'tangent'])
            sensitivity[kind+'-'+label] = dict(matched_points=len(matched), coarse_compared_points=len(coarse), fine_points=len(fine),
                same_full_field_start=same_start, extension_arc_interval_m=[0., .0512], parent_arc_interval_m=[.1024, .1536],
                endpoint=compare(*matched[-1], kind == 'native'), path_maxima=maxima,
                interpretation='same saved full-field anchor; local step halving' if label == 'localhalf' else 'continuation of prior coarse/fine prefixes; accumulated prior discretization differences included',
                full_extension_sensitivity_verified=False, history_converged=False)
        if kind == 'native':
            a = list(stream(evidence/'extended-native-tail128/samples.jsonl'))
            b = list(itertools.chain(stream(evidence/'extended-native-restarted64/samples.jsonl'), stream(evidence/'extended-native-restarted-tail64/samples.jsonl')))
            note = '128 continuous vs 64 plus actual checkpoint restart 64, beyond arc .512 m'
        else:
            b = list(stream(evidence/'extended-independent-restarted64/samples.jsonl'))
            indexes = {x['step'] for x in b}
            a = [x for x in stream(evidence/'extended-independent-long/samples.jsonl') if x['step'] in indexes]
            note = 'Actual checkpoint captured during ongoing 4096-step run; separate process 64-step continuation'
        restarts[kind] = dict(compared_samples=len(a), resumed_samples=len(b), all_saved_fields_exact=a == b and len(a) > 0,
            note=note, production_job_restart=False, restart_to_full_removal_verified=False)
    origin = g.run['continuation']['origin']
    zero_load = g.load0+g.lam0*g.holding
    zero = dict(relative_force_residual=float(np.linalg.norm(g.R.T@zero_load.ravel()-native_force(g, origin['material_states']))/max(np.linalg.norm(zero_load), 1.)))
    unload = read(evidence/'origin/reference/unload.json')
    izero = s.total_force(np.array(unload['material_state'])[active].reshape(-1, 14), np.array(unload['displacement_m']).ravel())
    nend = endpoints['extended-native-further']
    iend = endpoints['extended-independent-further']
    cross = dict(
        progress_absolute_difference=abs(nend['progress']-iend['progress']), displacement=diff(nend['displacement_m'], iend['displacement_m']),
        material={key:diff(fields(nend, True)[key], fields(iend, False)[key]) for key in fields(nend, True)},
        comparison_coordinate='same cumulative diagnostic arc .9216 m from distinct independently generated original seeds, not equal load or converged history')
    result = dict(schema='soil-fem-nc-further-branches/1', base_commit='90d24e4181016a1ba8211719ee16b284d82b5eef',
        passed=False, history_converged=False, static_limit_confirmed=False, verdict='inconclusive', production_adopted=False,
        native_highest_accepted_progress=g.lam0, native_accepted_increments=accepted_count(g.run), native_new_accepted_increments=0,
        full_removal_reached=False, full_release_freeze_verified=False, restart_to_completion_verified=False,
        condition=dict(nodes=54, tet10=12, material_gp=324, contact_gp=32, preload_kpa=5, pc_kpa=100, force_kn=[1, .5, 0]),
        remaining_region_lambda_zero_native=zero,
        remaining_region_lambda_zero_independent=dict(relative_force_residual=float(np.linalg.norm(s.reduce(izero-before))/max(np.linalg.norm(before), 1.)),
            note='input-only holding identity; not another solved equilibrium or full-release evidence'),
        paths=paths, local_root_searches=roots, step_sensitivity=sensitivity, restarts=restarts, both_map_endpoint_comparison=cross,
        windows_ci_for_new_head='unconfirmed', physical_stability='unconfirmed',
        energy_interpretation='stress work and symmetric spectra only; no specified hardening free energy or physical stability certificate',
        missing=['lambda=1 full removal', 'production monotone removal history beyond .9514000415802002', 'full-release freeze',
            'full removal independent plastic 3D evidence', 'whole extension/material/prehistory convergence', 'production restart to completion',
            'unexplored intervals and other/disconnected branches', 'physical stability', 'new-head Windows CI'])
    write(out/'measurements.json', result)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--evidence', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    measure(a.evidence.resolve(), a.output.resolve())
