"""Read-only spectra and fixed-history replay of the existing small NC case.

Native saved states are diagnostic inputs, never an independent reference.
No matrix is regularized or symmetrized for a solve. The symmetric-part
spectrum is labelled separately from eigenvalues of the actual matrix.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT / 'tools'), str(ROOT / 'interim/2026-11-09-nc-pseudo-arclength'),
                str(ROOT / 'interim/2026-11-10-nc-mean-stress-release')]
from linear_evidence_codec import read_rows
from mcc_rate_fem import terms


def read(p):
    return json.loads(p.read_text())


def rows(p):
    with gzip.open(p, 'rt') as stream:
        for line in stream:
            yield json.loads(line)


def write(p, value):
    p.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def checksum(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def spectrum(matrix, transform=None):
    a = np.asarray(matrix, dtype=float)
    if a.ndim != 2 or a.shape[0] != a.shape[1] or not np.isfinite(a).all():
        raise ValueError('finite square tangent required')
    left, singular, right = np.linalg.svd(a)
    v, u = right[-1].copy(), left[:, -1].copy()
    if v[np.argmax(np.abs(v))] < 0:
        v = -v
        u = -u
    eig, vectors = np.linalg.eig(a)
    order = np.argsort(np.abs(eig))
    symmetric = np.linalg.eigvalsh((a + a.T) / 2)
    result = dict(singular_values=singular.tolist(), condition_2=float(singular[0] / singular[-1]) if singular[-1] else None,
                  eigenvalues_actual_matrix=[dict(real=float(eig[i].real), imaginary=float(eig[i].imag)) for i in order],
                  symmetric_part_eigenvalues=symmetric.tolist(),
                  nonsymmetry_relative=float(np.linalg.norm(a-a.T) / max(np.linalg.norm(a), 1e-30)),
                  smallest_left_vector=u.tolist(), smallest_right_vector=v.tolist(),
                  right_singular_equation_residual=float(np.linalg.norm(a @ v - singular[-1]*u) / max(singular[0], 1e-30)),
                  left_singular_equation_residual=float(np.linalg.norm(a.T @ u - singular[-1]*v) / max(singular[0], 1e-30)),
                  static_limit_certified=False)
    if transform is not None:
        result['right_mode_nodal_displacement'] = (np.asarray(transform) @ v).reshape(-1, 3).tolist()
    return result


def state_array(sample, kind, active, pi):
    if kind == 'independent':
        return np.asarray(sample['states'])[active].reshape(-1, 14)
    return np.array([s['effective_stress_compression_positive_kpa']
                     + [np.log(s['preconsolidation_pressure_kpa']), np.log(1+s['void_ratio'])]
                     + s['plastic_strain_compression_positive'] for s in (sample['states'][i] for i in pi)])


def material_fields(state, parameters, strain, indices, weights, B, transform):
    yf, elastic, n, dn, den = terms(state, parameters)
    loading = np.einsum('qa,qa->q', dn, strain)
    plastic = (yf >= -2e-9) & (loading > 1e-12 * np.maximum(np.abs(den), 1.) / np.maximum(np.linalg.norm(n, axis=1), 1e-30))
    if np.any(den[plastic] <= 0):
        raise ValueError('nonpositive active consistency denominator')
    tangent = elastic.copy()
    tangent[plastic] -= np.einsum('qa,qb,q->qab', dn[plastic], dn[plastic], 1/den[plastic])
    assembled = np.einsum('qai,qab,qbj,q->ij', B, tangent, B, weights, optimize=True)
    elastic_global = np.einsum('qai,qab,qbj,q->ij', B, elastic, B, weights, optimize=True)
    spec = spectrum(assembled, transform)
    left = np.asarray(spec['smallest_left_vector'])
    right = np.asarray(spec['smallest_right_vector'])
    el, er = B @ left, B @ right
    contributions = weights * np.einsum('qa,qab,qb->q', el, tangent, er)
    energy = weights * np.einsum('qa,qab,qb->q', er, tangent, er)
    p = state[:, :3].mean(axis=1)
    pc = np.exp(state[:, 6])
    q2 = yf*parameters[:, 0]*pc**2-parameters[:, 0]*p*(p-pc)
    local = []
    for i, gp in enumerate(indices):
        sv = np.linalg.svd(tangent[i], compute_uv=False)
        ev = np.linalg.eigvalsh((tangent[i]+tangent[i].T)/2)
        local.append(dict(gp=int(gp), element=int(gp)//27, state=state[i].tolist(), p_kpa=float(p[i]),
                          q_kpa=float(np.sqrt(max(q2[i], 0.))), pc_kpa=float(pc[i]), p_over_pc=float(p[i]/pc[i]),
                          normalized_yield=float(yf[i]), branch_from_observed_increment='plastic' if plastic[i] else 'elastic',
                          loading_numerator=float(loading[i]), consistency_denominator=float(den[i]),
                          denominator_over_elastic_term=float(den[i]/max(abs(n[i]@dn[i]), 1e-30)),
                          local_singular_values=sv.tolist(), local_symmetric_eigenvalues=ev.tolist(),
                          right_mode_engineering_strain=er[i].tolist(), left_mode_engineering_strain=el[i].tolist(),
                          left_right_mode_contribution=float(contributions[i]), right_mode_work=float(energy[i])))
    spec.update(elastic_control=spectrum(elastic_global),
                plastic_points=int(plastic.sum()), material_points=local,
                mode_contribution_sum=float(contributions.sum()),
                absolute_mode_contribution_sum=float(np.abs(contributions).sum()),
                mode_work_sum=float(energy.sum()),
                branch_selection='observed accepted increment sign; not a new global equilibrium solution')
    return spec


def native_attempts(folder):
    linear_events = []
    for attempt in rows(folder / 'attempts.jsonl.gz'):
        label = attempt.get('label', 'step-' + str(attempt.get('step')))
        for event in attempt['events']:
            if 'reduced_increment_m' in event:
                linear_events.append((label, event['iteration'], attempt, event))
    return linear_events


def replay_requests(evidence, binary, out, name, samples, attempts, actual):
    controls = read(evidence / name / 'controls.json')
    pi = controls['active_gp']
    B = np.asarray(controls['B'])
    ground = read(evidence / 'prepared/origin/prior/fresh-head-midpoint/input/ground.json')
    mesh = read(evidence / 'prepared/origin/prior/fresh-head-midpoint/input/mesh.json')
    pars = {x['id']: dict(x['mcc'], integration='midpoint') for x in ground['materials']}
    events = native_attempts(evidence / name)
    matched = events[-len(actual):]
    if len(matched) != len(actual):
        raise ValueError('missing native Newton events')
    selected = {}
    for j, (row, event) in enumerate(zip(actual, matched)):
        label, it, attempt, iteration = event
        if row['label'] != label or row['iteration'] != it:
            raise ValueError('native actual matrix / Newton event ordering differs')
        if j in (0, len(actual)-1) or (not attempt['accepted'] and it == 0):
            selected[j] = (row, attempt, iteration)
    # Keep a bounded set including the final failure and the worst retained matrix.
    worst = int(np.argmax([np.linalg.cond(np.asarray(x['matrix'])) for x in actual]))
    selected[worst] = (actual[worst], matched[worst][2], matched[worst][3])
    chosen = sorted(set([0, len(actual)-1, worst] + list(selected)[-2:]))
    by_coordinate = {s['coordinate']: s for s in samples}
    requests = []
    witnesses = []
    for j in chosen:
        row, attempt, iteration = selected[j]
        old = by_coordinate[attempt['start']]
        du = np.asarray(iteration['reduced_increment_m'])
        strain = -B @ du
        label = 'terminal-' + str(j)
        for partitions in (1, 2, 4):
            for i, gp in enumerate(pi):
                requests.append(dict(label=f'{label}-p{partitions}', point=gp,
                                     parameters=pars[mesh['material_ids'][gp//27]], state=old['states'][gp],
                                     strain=strain[i].tolist(), partitions=partitions,
                                     widths=[1e-7, 1e-9, 1e-11] if partitions == 1 else []))
        witnesses.append(dict(label=label, matrix=row, attempt=attempt, event=iteration,
                              old_sample=old, old_sample_sha256=hashlib.sha256(json.dumps(old, sort_keys=True).encode()).hexdigest()))
    ip = out / (name + '-requests.json')
    op = out / (name + '-response.json')
    write(ip, requests)
    subprocess.run([str(binary), str(ip), str(op)], check=True, timeout=120)
    response = read(op)
    results = []
    for witness in witnesses:
        base = [x for x in response['points'] if x['label'] == witness['label']+'-p1']
        errors = [x for x in base if 'error' in x]
        record = dict(label=witness['label'], attempt=witness['attempt'], event=witness['event'],
                      old_sample_sha256=witness['old_sample_sha256'], local_errors=errors)
        if errors:
            results.append(record)
            continue
        if [x['point'] for x in base] != pi or not all(x['ad_primal_exact'] and x['substep_replay_exact'] for x in base):
            raise ValueError('native material replay not exact')
        D = np.asarray([x['implicit_tangent'] for x in base])
        K = np.einsum('qai,qab,qbj,q->ij', B, D, B, np.asarray(controls['weights']), optimize=True)
        original = np.asarray(witness['matrix']['matrix'])
        record.update(replayed_tangent=spectrum(K, controls['transform']),
                      relative_to_archived_actual_tangent=float(np.linalg.norm(K-original)/np.linalg.norm(original)),
                      points=base, partition_sensitivity={})
        for count in (2, 4):
            values = [x for x in response['points'] if x['label'] == witness['label']+f'-p{count}']
            record['partition_sensitivity'][str(count)] = values
        results.append(record)
    return results


def run(evidence, out, binary=None):
    out.mkdir(parents=True, exist_ok=False)
    proof = read(evidence / 'extraction.json')
    # Bind every retained input byte again, independently of the extraction run.
    for name, expected in proof['retained_entries'].items():
        p = evidence / name
        if checksum(p) != expected['sha256'] or p.stat().st_size != expected['size']:
            raise ValueError('diagnostic input changed: ' + name)
    ground = read(evidence / 'prepared/origin/prior/fresh-head-midpoint/input/ground.json')
    mesh = read(evidence / 'prepared/origin/prior/fresh-head-midpoint/input/mesh.json')
    materials = {x['id']: x['mcc'] for x in ground['materials']}
    from probe_nc_midpoint import inputs
    case = read(evidence / 'prepared/origin/prior/fresh-head-midpoint/input/stages.json')
    if (mesh, ground, case) != inputs('midpoint'):
        raise ValueError('original 5 kPa small NC input differs')
    seed = next(x['event'] for x in read(evidence/'prepared/origin/field-records.json') if x['label'] == 'terminal-rejected:0')
    contact_points = seed['joint_report']['points']
    if len(contact_points) != 32 or any(x['active'] for x in contact_points) or np.any(seed['joint_internal']):
        raise ValueError('post-cut contact inactivity is not established')
    upper = next(x['id'] for x in mesh['regions'] if x['body'] == 'UpperClay')
    expected_active = np.asarray(mesh['region_ids']) != upper
    result = dict(passed=False, starting_commit=proof['starting_commit'], archive_sha256=proof['archive_sha256'],
                  full_removal_accepted=False, static_limit_confirmed=False, independent_full_field_completed=False,
                  paths={}, material_replays={}, production_adopted=False)
    names = sorted(p.name for p in evidence.iterdir() if (p/'controls.json').exists())
    for name in names:
        folder = evidence / name
        controls = read(folder / 'controls.json')
        samples = list(rows(folder / 'samples.jsonl.gz'))
        actual = list(read_rows(folder / 'terminal-linear-systems.jsonl.f64.gz'))
        kind = 'native' if name.startswith('native-') else 'independent'
        active = np.asarray(controls['active_elements'])
        if not np.array_equal(active, expected_active):
            raise ValueError('existing complete UpperClay removal layout differs')
        pi = controls['active_gp']
        B, weights, R = (np.asarray(controls[k]) for k in ('B', 'weights', 'transform'))
        parameters = np.array([[materials[mesh['material_ids'][i//27]][k] for k in
                                ('critical_state_m', 'lambda', 'kappa', 'poisson')] for i in pi])
        parameters[:, 0] **= 2
        frames = []
        indices = sorted(set([max(i for i, s in enumerate(samples) if s['coordinate'] <= min(x, samples[-1]['coordinate']))
                              for x in (.5, .9, .95, .97, .9744, 1.)]))
        for index in indices:
            sample = samples[index]
            if index == 0:
                continue
            prior = samples[index-1]
            displacement = np.asarray(sample['displacement_m'])-np.asarray(prior['displacement_m'])
            du = np.linalg.lstsq(R, displacement.ravel(), rcond=None)[0]
            strain = -B @ du
            fields = material_fields(state_array(sample, kind, active, pi), parameters, strain, pi, weights, B, R)
            fields.update(coordinate=sample['coordinate'], step=sample['step'],
                          saved_force_residual=sample['reassembled_relative_force_residual'],
                          removed_material_frozen=sample['removed_material_frozen'],
                          frozen_contact_history=sample['contacts'], active_contacts=0,
                          contact_activity_basis='all owners UpperClay inactive after the existing cut; saved stick/force are frozen history')
            frames.append(fields)
        retained = []
        for row in actual:
            spec = spectrum(row['matrix'], R)
            spec.update(label=row['label'], iteration=row['iteration'],
                        rhs_projection_left=float(np.asarray(spec['smallest_left_vector']) @ np.asarray(row['rhs'])),
                        recorded_true_linear_residual=row['true_relative_residual'])
            retained.append(spec)
        # Preserve full spectra/modes for the worst actual matrix and the final matrix.
        worst = max(range(len(retained)), key=lambda i: retained[i]['condition_2'] or float('inf'))
        result['paths'][name] = dict(report=read(folder/'report.json'), frames=frames,
                                    retained_actual_matrix_count=len(actual), worst_retained_matrix=retained[worst],
                                    final_retained_matrix=retained[-1],
                                    retained_matrix_summary=[{k: x[k] for k in ('label', 'iteration', 'condition_2',
                                      'recorded_true_linear_residual', 'rhs_projection_left')} for x in retained],
                                    independent_input_only_history=kind == 'independent',
                                    all_saved_contacts_identical=all(x['contacts'] == samples[0]['contacts'] for x in samples))
        if binary is not None and name == 'native-complement-refined128':
            result['material_replays'][name] = replay_requests(evidence, binary, out, name, samples,
                                                              list(rows(folder/'attempts.jsonl.gz')), actual)
        print('diagnosed ' + name, flush=True)
    write(out / 'measurements.json', result)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--evidence', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--binary', type=Path)
    a = p.parse_args()
    run(a.evidence.resolve(), a.output.resolve(), a.binary.resolve() if a.binary else None)
