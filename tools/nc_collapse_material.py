"""Extend saved GP states using the existing Rust MCC return-map observer.

Stress-ray inverse solves and strain-direction probes are diagnostic experiments.
A failed local corrector is not labelled a constitutive limit. No state or
constitutive equation is substituted in the production solver.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from nc_collapse_classify import read, write, rows


STRESS = 'effective_stress_compression_positive_kpa'


def batch(binary, out, counter, requests):
    ip = out/f'{counter:05d}-request.json'; op = out/f'{counter:05d}-result.json'
    write(ip, requests)
    result = subprocess.run([str(binary), str(ip), str(op)], capture_output=True, text=True, timeout=30)
    (out/f'{counter:05d}.log').write_text(result.stdout+result.stderr, encoding='utf-8')
    result.check_returncode()
    points = read(op)['points']
    if len(points) != len(requests):
        raise ValueError('material response inventory differs')
    return points


def request(pi, par, state, strain, label):
    return dict(label=label, point=int(pi), parameters=par, state=state, strain=np.asarray(strain).tolist(), partitions=1, widths=[])


def run(folder, binary, output, strain_steps=160, strain_increment=1e-5, stress_steps=60, stress_increment=.1):
    if (type(strain_steps) is not int or not 1 <= strain_steps <= 2000 or type(stress_steps) is not int or not 1 <= stress_steps <= 2000
            or not np.isfinite(strain_increment) or not 0 < strain_increment <= .001
            or not np.isfinite(stress_increment) or not 0 < stress_increment <= 1):
        raise ValueError('bounded positive material probe controls required')
    from verify_high_order import prepare
    # Check original files before constructing any trial states.
    manifest = read(folder/'evidence-sha256.json')
    for name, sha in manifest.items():
        if hashlib.sha256((folder/name).read_bytes()).hexdigest() != sha:
            raise ValueError('source evidence differs')
    output.mkdir(parents=True, exist_ok=False); raw = output/'raw'; raw.mkdir()
    g = read(folder/'geometry.json'); samples = rows(folder/'samples.jsonl'); tangents = rows(folder/'gp-tangents.jsonl')
    peak = max(range(len(samples)), key=lambda i: samples[i]['progress']); state = samples[peak]
    previous = samples[peak-1] if peak else read(folder/'initial.json')['payload']
    r = np.asarray(g['transform']); b = np.asarray(g['reduced_b_rows'])
    du = np.linalg.lstsq(r, (np.asarray(state['displacement_m'])-previous['displacement_m']).ravel(), rcond=None)[0]
    de = -np.einsum('qai,i->qa', b, du)
    q, offsets = prepare(g['mesh']); active = np.flatnonzero(g['active_elements'])
    ids = np.concatenate([np.arange(offsets[i], offsets[i+1]) for i in active])
    positions = np.concatenate([q[i]['position'] for i in active])
    elements = np.concatenate([np.full(offsets[i+1]-offsets[i], i) for i in active])
    materials = {x['id']: x['mcc'] for x in g['ground']['materials']}
    pars = [dict(materials[g['mesh']['material_ids'][ei]], integration=g['case']['settings']['mcc_integration']) for ei in elements]
    start = [state['states'][i] for i in ids]
    sigma0 = np.asarray([s[STRESS] for s in start])
    ds = sigma0-np.asarray([previous['states'][i][STRESS] for i in ids])
    norm = np.linalg.norm(ds, axis=1)
    if np.any(norm == 0) or np.any(np.linalg.norm(de, axis=1) == 0):
        raise ValueError('adjacent path has zero direction')
    stress_direction = ds/norm[:, None]; strain_direction = de/np.linalg.norm(de, axis=1)[:, None]
    counter = 0; t0 = time.monotonic(); records = []
    old = list(start); trajectories = [[] for _ in ids]; running = set(range(len(ids)))
    prior_projection = np.sum(sigma0*stress_direction, axis=1); peaks = [None]*len(ids); failures = [None]*len(ids)
    for step in range(strain_steps):
        js = sorted(running)
        if not js:
            break
        requests = [request(ids[j], pars[j], old[j], strain_increment*strain_direction[j], f'strain-{step}') for j in js]
        points = batch(binary, raw, counter, requests); counter += 1
        for j, p in zip(js, points):
            if 'error' in p:
                failures[j] = dict(step=step, error=p['error']); running.remove(j); continue
            if not p['ad_primal_exact'] or not p['substep_replay_exact']:
                raise ValueError('observer changed primal return map')
            s = p['updated_state']; sigma = np.asarray(s[STRESS]); c = np.asarray(p['implicit_tangent'])
            projection = float(sigma@stress_direction[j]); slope = (projection-prior_projection[j])/strain_increment
            trajectories[j].append(dict(step=step+1, strain_length=(step+1)*strain_increment,
                stress_kpa=sigma.tolist(), pc_kpa=s['preconsolidation_pressure_kpa'], void_ratio=s['void_ratio'],
                projected_stress_kpa=projection, secant_directional_stiffness_kpa=float(slope),
                tangent_directional_projection_kpa=float(stress_direction[j]@c@strain_direction[j]),
                singular_min_kpa=float(np.linalg.svd(c, compute_uv=False)[-1])))
            if peaks[j] is None and slope <= 0:
                peaks[j] = dict(strain_length_interval=[step*strain_increment, (step+1)*strain_increment],
                    projection_kpa=prior_projection[j], sample=step, criterion='first nonpositive projected stress increment along fixed adjacent strain direction',
                    absolute_constitutive_limit=False)
            prior_projection[j] = projection; old[j] = s
    # Stress-controlled inverse map. A failure remains an unresolved corrector.
    old = list(start); cprior = np.asarray(tangents[peak]['local_tangents']).copy(); running = set(range(len(ids)))
    rays = [[] for _ in ids]; ray_failures = [None]*len(ids)
    for step in range(stress_steps):
        js = sorted(running)
        if not js:
            break
        increments = {j: np.linalg.lstsq(cprior[j], stress_increment*stress_direction[j], rcond=None)[0] for j in js}
        active_newton = set(js); next_states = {}; next_tangents = {}
        for iteration in range(14):
            js = sorted(active_newton)
            if not js:
                break
            requests = [request(ids[j], pars[j], old[j], increments[j], f'stress-{step}-it-{iteration}') for j in js]
            points = batch(binary, raw, counter, requests); counter += 1
            for j, p in zip(js, points):
                if 'error' in p:
                    ray_failures[j] = dict(step=step, iteration=iteration, reason='return_map_error', error=p['error']); active_newton.remove(j); running.remove(j); continue
                if not p['ad_primal_exact'] or not p['substep_replay_exact']:
                    raise ValueError('observer changed primal return map')
                target = sigma0[j]+(step+1)*stress_increment*stress_direction[j]
                s = p['updated_state']; residual = target-np.asarray(s[STRESS]); c = np.asarray(p['implicit_tangent'])
                rel = np.linalg.norm(residual)/max(np.linalg.norm(target), 1.)
                if rel <= 1e-8:
                    next_states[j] = s; next_tangents[j] = c; active_newton.remove(j)
                    rays[j].append(dict(stress_distance_kpa=(step+1)*stress_increment, relative_stress_residual=float(rel),
                        strain_increment=increments[j].tolist(), pc_kpa=s['preconsolidation_pressure_kpa'],
                        singular_min_kpa=float(np.linalg.svd(c, compute_uv=False)[-1]),
                        state=s))
                else:
                    correction = np.linalg.lstsq(c, residual, rcond=None)[0]
                    scale = max(1., np.sum(np.abs(increments[j]+correction))/.001)
                    increments[j] = (increments[j]+correction)/scale
        for j in active_newton:
            ray_failures[j] = dict(step=step, reason='local_stress_corrector_iteration_limit'); running.remove(j)
        for j in next_states:
            old[j] = next_states[j]; cprior[j] = next_tangents[j]
    for j, pi in enumerate(ids):
        records.append(dict(point=int(pi), element=int(elements[j]), position_m=positions[j].tolist(),
            initial_state=start[j], adjacent_stress_direction=stress_direction[j].tolist(), adjacent_strain_direction=strain_direction[j].tolist(),
            directional_peak=peaks[j], strain_probe_failure=failures[j], strain_probe=trajectories[j],
            stress_ray=rays[j], stress_ray_failure=ray_failures[j],
            stress_ray_limit_proved=False, constitutive_zero_mode_proved=False,
            limit_at_fold=False if peaks[j] is not None and peaks[j]['sample']>0 else None))
    result = dict(schema='nc-collapse-material/1', lambda_fold=state['progress'], records=records,
        controls=dict(strain_steps=strain_steps, strain_increment=strain_increment, stress_steps=stress_steps, stress_increment_kpa=stress_increment,
                      local_stress_residual=1e-8, local_maximum_iterations=14),
        gp_count=len(ids), directional_peaks_observed=sum(p is not None for p in peaks),
        stress_ray_corrector_failures=sum(p is not None for p in ray_failures),
        absolute_material_limits_proved=0, limit_at_fold_count=None,
        passed=False, nc01_complete=False, elapsed_seconds=time.monotonic()-t0,
        observer_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
        source_manifest_sha256=hashlib.sha256((folder/'evidence-sha256.json').read_bytes()).hexdigest(),
        limitation='Finite stress-ray/strain-direction horizons. A projected peak on a prescribed local strain direction is not a stress cap under every loading path. Failed inverse correctors do not establish a material limit.')
    write(output/'material-limits.json', result)
    return {k: result[k] for k in ('gp_count', 'directional_peaks_observed', 'stress_ray_corrector_failures', 'elapsed_seconds')}


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--path', type=Path, required=True); p.add_argument('--binary', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--strain-steps', type=int, default=160); p.add_argument('--stress-steps', type=int, default=60)
    args = p.parse_args(); print(json.dumps(run(args.path, args.binary, args.output, strain_steps=args.strain_steps, stress_steps=args.stress_steps)))
