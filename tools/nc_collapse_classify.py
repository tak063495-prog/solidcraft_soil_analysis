"""Read-only nonuniform NC fold diagnostics; no NC-01 acceptance criterion.

Distinguishes nonsymmetric K weak vectors, rectangular equilibrium nullspace,
and bordered arc system. Hill work and acoustic data retain their branch and
finite-increment limitations. Saved fields are never edited or resealed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from nc_collapse_bounds import discrete_bound, mechanism_bound
from mcc_acoustic_diagnostic import branch, acoustic_tensor
from mcc_acoustic_exact import criterion


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def write(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def rows(path):
    return [json.loads(line) for line in path.read_text(encoding='utf-8').splitlines()]


def matrix_diagnostic(k, load_column, plane):
    k = np.asarray(k, dtype=float); h = np.asarray(load_column, dtype=float)
    if k.shape != (len(h), len(h)) or not np.isfinite(k).all() or not np.isfinite(h).all() or np.linalg.norm(h) == 0:
        raise ValueError('finite square K and nonzero load column required')
    u, singular, vh = np.linalg.svd(k)
    left = u[:, -1]; right = vh[-1]
    rectangular = np.column_stack((k, h))
    _, sr, vr = np.linalg.svd(rectangular, full_matrices=True)
    null = vr[-1]
    augmented = np.vstack((rectangular, plane))
    sa = np.linalg.svd(augmented, compute_uv=False)
    scale = max(float(singular[0]), 1e-30)
    return dict(k_singular_values=singular.tolist(), k_numerical_rank=int(np.linalg.matrix_rank(k)),
                k_right_weak_mode=right.tolist(), k_left_weak_mode=left.tolist(),
                k_relative_right_residual=float(np.linalg.norm(k@right)/scale),
                k_relative_left_residual=float(np.linalg.norm(k.T@left)/scale),
                left_load_cosine=float(left@h/np.linalg.norm(h)), right_load_cosine=float(right@h/np.linalg.norm(h)),
                rectangular_rank=int(np.linalg.matrix_rank(rectangular)), rectangular_singular_values=sr.tolist(),
                rectangular_null_vector=null.tolist(), rectangular_null_lambda=float(null[-1]),
                rectangular_null_relative_residual=float(np.linalg.norm(rectangular@null)/max(np.linalg.norm(rectangular), 1e-30)),
                augmented_rank=int(np.linalg.matrix_rank(augmented)), augmented_singular_values=sa.tolist(),
                augmented_exact_right_null_exists=bool(np.linalg.matrix_rank(augmented) < len(augmented)),
                exact_k_zero_mode_proved=False,
                interpretation='A full-rank rectangular equilibrium block always has a 1D path nullspace; it is not a singular bordered system. For nonsymmetric K the left load cosine controls solvability. Finite samples cannot prove an exact critical point.')


def acoustic_minimum(data):
    """Multi-start directional minimisation, not an all-direction certificate."""
    def normal(x):
        return np.array([np.sin(x[0])*np.cos(x[1]), np.sin(x[0])*np.sin(x[1]), np.cos(x[0])])
    def objective(x):
        return float(np.linalg.eigvalsh(acoustic_tensor(data, normal(x), plastic=True))[0])
    seeds = [(t, a) for t in (0., np.pi/4, np.pi/2) for a in (0., np.pi/2, np.pi, 3*np.pi/2)]
    results = [minimize(objective, x, method='Nelder-Mead', options={'xatol': 1e-9, 'fatol': 1e-9, 'maxiter': 200}) for x in seeds]
    best = min(results, key=lambda r: r.fun)
    return dict(minimum_found_kpa=float(best.fun), normal=normal(best.x).tolist(), starts=len(seeds),
                best_optimizer_success=bool(best.success), global_minimum_certified=False)


def run(folder, output, solve_bounds=True):
    from nc_physical_limit_audit import stress_values
    manifest = read(folder/'evidence-sha256.json')
    for name, sha in manifest.items():
        if hashlib.sha256((folder/name).read_bytes()).hexdigest() != sha:
            raise ValueError('saved evidence changed: '+name)
    g = read(folder/'geometry.json'); samples = rows(folder/'samples.jsonl'); gp = rows(folder/'gp-tangents.jsonl')
    linear = rows(folder/'linear-systems.jsonl')
    b = np.asarray(g['reduced_b_rows']); weights = np.asarray(g['weights_m3']); r = np.asarray(g['transform'])
    loads = read(folder/'original-load.json'); final = r.T@loads['final_load_kn']; holding = r.T@loads['holding_load_kn']
    peak = max(range(len(samples)), key=lambda i: samples[i]['progress']); sample = samples[peak]
    from verify_high_order import prepare
    q, offsets = prepare(g['mesh']); active = np.flatnonzero(g['active_elements'])
    ids = np.concatenate([np.arange(offsets[i], offsets[i+1]) for i in active])
    positions = np.concatenate([q[i]['position'] for i in active])
    elements = np.concatenate([np.full(offsets[i+1]-offsets[i], i) for i in active])
    materials = {x['id']: x['mcc'] for x in g['ground']['materials']}
    parameters = [materials[g['mesh']['material_ids'][ei]] for ei in elements]
    previous = samples[peak-1] if peak else read(folder/'initial.json')['payload']
    du = np.linalg.lstsq(r, (np.asarray(sample['displacement_m'])-previous['displacement_m']).ravel(), rcond=None)[0]
    de = -np.einsum('qai,i->qa', b, du)
    states = [sample['states'][i] for i in ids]
    oldstates = [previous['states'][i] for i in ids]
    sigma = np.array([s['effective_stress_compression_positive_kpa'] for s in states])
    ds = sigma-np.array([s['effective_stress_compression_positive_kpa'] for s in oldstates])
    c = np.asarray(gp[peak]['local_tangents'])
    systems = [x for x in linear if x['step'] == sample['step'] and x['kind'] == 'tangent']
    a = np.asarray(systems[-1]['matrix']); k = a[:-1, :-1]
    md = matrix_diagnostic(k, holding, a[-1])
    mode = np.asarray(md['k_right_weak_mode'])
    if mode@np.asarray(sample['tangent'])[:-1] < 0:
        mode = -mode
    mode_de = -np.einsum('qai,i->qa', b, mode)
    work = np.einsum('qa,qa->q', ds, de)
    mode_work = np.einsum('qa,qab,qb->q', mode_de, c, mode_de)
    pc = np.array([s['preconsolidation_pressure_kpa'] for s in states]); m = np.array([p['critical_state_m'] for p in parameters])
    records = []
    for j, (pi, state, par) in enumerate(zip(ids, states, parameters)):
        data = branch(sigma[j], pc[j], state['void_ratio'], par)
        elastic_loading = float(data['normal']@data['elastic']@de[j])
        plastic_loading = elastic_loading > 0 and abs(data['normalized_yield']) <= 1e-8
        continuum = data['elastic'] if not plastic_loading else data['elastic']-np.outer(data['cn'], data['cn'])/data['denominator_kpa3']
        acoustic = None
        if data['denominator_kpa3'] > 0:
            acoustic = dict(plastic_branch_all_direction_sign=criterion(data), plastic_branch_minimum=acoustic_minimum(data))
        records.append(dict(point=int(pi), element=int(elements[j]), position_m=positions[j].tolist(),
                            **stress_values(state, m[j]),
                            q_gt_Mp_ellipse_consistent=bool(stress_values(state, m[j])['q_over_Mp'] <= 1 or (sigma[j,:3].mean() < pc[j]/2 and stress_values(state,m[j])['normalized_ellipse_yield'] <= 1e-8)),
                            adjacent_stress_increment_kpa=ds[j].tolist(), adjacent_strain_increment=de[j].tolist(),
                            second_work_actual_increment_kpa=float(work[j]), second_work_weak_mode_kpa=float(mode_work[j]),
                            local_algorithmic_singular_values_kpa=np.linalg.svd(c[j], compute_uv=False).tolist(),
                            local_symmetric_eigen_min_kpa=float(np.linalg.eigvalsh((c[j]+c[j].T)/2)[0]),
                            continuum_loading_indicator=elastic_loading,
                            continuum_selected_branch='plastic_loading' if plastic_loading else 'elastic_or_unloading',
                            continuum_second_work_kpa=float(de[j]@continuum@de[j]), acoustic=acoustic))
    output.mkdir(parents=True, exist_ok=False)
    write(output/'gp-classification.json', dict(lambda_fold=sample['progress'], records=records,
          definition=dict(p='compression-positive trace(sigma)/3', q='sqrt(3/2 s:s), physical shear stresses',
                          strain='compression-positive engineering xx yy zz xy yz zx', M='constant parameter; no Lode angle dependency; no phi conversion at runtime'),
          q_gt_Mp_count=int(sum(x['q_over_Mp']>1 for x in records)),
          dry_count=int(sum(x['stress_side']=='dry' for x in records))))
    write(output/'stability.json', dict(lambda_fold=sample['progress'], matrix=md,
        actual_increment_negative_gp=int(sum(work<0)), weak_mode_negative_gp=int(sum(mode_work<0)),
        actual_increment_integral=float(weights@work), weak_mode_integral=float(weights@mode_work),
        weak_mode_global_contraction=float(mode@k@mode),
        full_mode=(r@mode).reshape(-1,3).tolist(),
        path_matrix_diagnostics=[dict(step=s['step'], lambda_value=s['progress'],
            **matrix_diagnostic(np.asarray(next(x for x in linear if x['step']==s['step'] and x['kind']=='tangent')['matrix'])[:-1,:-1],
                holding, np.asarray(next(x for x in linear if x['step']==s['step'] and x['kind']=='tangent')['matrix'])[-1])) for s in samples],
        strict_zero_mode_proved=False, classification='D', passed=False,
        limitation='Finite-increment algorithmic work and frozen continuum acoustic branch are different objects; local negative work alone does not prove collapse or exclude other equilibria.'))
    if solve_bounds:
        bounds = {kind: discrete_bound(b, weights, final, holding, pc, m, kind) for kind in ('ellipse', 'cone', 'intersection')}
        bounds['weak_mode_upper'] = mechanism_bound(b, weights, mode, final, holding, pc, m)
        if 'dual_equilibrium' in bounds['ellipse']:
            bounds['dual_mode_upper'] = mechanism_bound(b, weights, bounds['ellipse']['dual_equilibrium'], final, holding, pc, m)
        bounds['fold_inside_discrete_bounds'] = bool(bounds['ellipse'].get('lambda_LB', 1.) <= sample['progress'] <= bounds['weak_mode_upper'].get('lambda_UB', 0.))
        bounds.update(lambda_fold=sample['progress'], nc01_complete=False, passed=False)
        write(output/'bounds.json', bounds)
    write(output/'source-binding.json', dict(source_manifest_sha256=hashlib.sha256((folder/'evidence-sha256.json').read_bytes()).hexdigest(),
          diagnostic_only=True, nc01_complete=False, path=str(folder), peak_step=sample['step'], tool_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          bounds_tool_sha256=hashlib.sha256(Path(__file__).with_name('nc_collapse_bounds.py').read_bytes()).hexdigest()))
    return dict(peak=sample['progress'], gp=len(records), negative_actual=int(sum(work<0)), negative_mode=int(sum(mode_work<0)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--path', type=Path, required=True); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); print(json.dumps(run(args.path, args.output)))
