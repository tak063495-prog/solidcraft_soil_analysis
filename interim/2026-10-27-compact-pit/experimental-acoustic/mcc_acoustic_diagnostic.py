"""Read-only continuum MCC acoustic diagnostics, separate from FEM acceptance."""
import numpy as np


def branch(stress, pc, void_ratio, material):
    stress = np.asarray(stress, dtype=float)
    if stress.shape != (6,) or not np.isfinite(stress).all():
        raise ValueError('six finite compression-positive stresses required')
    m, lam, kap, nu = [float(material[k]) for k in ('critical_state_m', 'lambda', 'kappa', 'poisson')]
    if not np.isfinite([pc, void_ratio, m, lam, kap, nu]).all() or not (pc > 0 and void_ratio > 0 and m > 0 and lam > kap > 0 and -1 < nu < .5):
        raise ValueError('admissible MCC parameters, positive pc and void ratio required')
    mean = stress[:3].mean()
    if mean <= 0:
        raise ValueError('positive mean effective stress required; no pressure clipping')
    dev = stress.copy()
    dev[:3] -= mean
    q2 = 1.5 * (dev[:3] @ dev[:3]) + 3 * (dev[3:] @ dev[3:])
    bulk = (1 + void_ratio) * mean / kap
    shear = 3 * bulk * (1 - 2 * nu) / (2 * (1 + nu))
    elastic = np.zeros((6, 6))
    elastic[:3, :3] = bulk - 2 * shear / 3
    for i in range(3):
        elastic[i, i] += 2 * shear
        elastic[i + 3, i + 3] = shear
    normal = dev * np.array([3., 3., 3., 6., 6., 6.])
    trace = m * m * (2 * mean - pc)
    normal[:3] += trace / 3
    cn = elastic @ normal
    hardening = m * m * mean * pc * (1 + void_ratio) / (lam - kap) * trace
    denominator = normal @ cn + hardening
    return dict(mean_kpa=mean, normalized_yield=(q2 + m*m*mean*(mean-pc))/(m*m*pc*pc),
                bulk_kpa=bulk, shear_kpa=shear, elastic=elastic, normal=normal, cn=cn,
                consistency_H_kpa3=hardening, denominator_kpa3=denominator)


def acoustic_tensor(data, direction, *, plastic):
    n = np.asarray(direction, dtype=float)
    if n.shape != (3,) or not np.isfinite(n).all() or abs(np.linalg.norm(n) - 1) > 1e-10:
        raise ValueError('a finite unit direction required')
    bulk, shear = data['bulk_kpa'], data['shear_kpa']
    q = shear * np.eye(3) + (bulk + shear/3) * np.outer(n, n)
    if plastic:
        den = data['denominator_kpa3']
        if not np.isfinite(den) or den <= 0:
            raise ValueError('plastic consistency denominator is not positive')
        cn = data['cn']
        tensor = np.array([[cn[0], cn[3], cn[5]], [cn[3], cn[1], cn[4]], [cn[5], cn[4], cn[2]]])
        a = tensor @ n
        q -= np.outer(a, a) / den
    return q


def directions(polar_subdivisions=12, azimuth_subdivisions=48):
    if type(polar_subdivisions) is not int or type(azimuth_subdivisions) is not int or not (2 <= polar_subdivisions <= 90 and 4 <= azimuth_subdivisions <= 360):
        raise ValueError('bounded integer angular subdivisions required')
    out = [[np.sin(p)*np.cos(a), np.sin(p)*np.sin(a), np.cos(p)]
           for p in np.linspace(0, np.pi/2, polar_subdivisions+1)
           for a in np.linspace(0, 2*np.pi, azimuth_subdivisions, endpoint=False)]
    out = np.unique(np.round(out, 14), axis=0)
    return out / np.linalg.norm(out, axis=1)[:, None]


def sample(data, normals):
    normals = np.asarray(normals, dtype=float)
    if normals.ndim != 2 or normals.shape[1] != 3 or not 1 <= len(normals) <= 32768:
        raise ValueError('bounded nonempty direction inventory required')
    best = None
    for normal in normals:
        q = acoustic_tensor(data, normal, plastic=True)
        values, vectors = np.linalg.eigh(q)
        if best is None or values[0] < best['eigenvalue_kpa']:
            vector = vectors[:, 0]
            best = dict(eigenvalue_kpa=float(values[0]), normalized_by_shear=float(values[0]/data['shear_kpa']),
                        normal=normal.tolist(), polarization=vector.tolist(),
                        rayleigh_kpa=float(vector @ q @ vector))
    best.update(diagnostic_only=True, actual_loading_branch_not_determined=True,
                positive_samples_do_not_prove_all_direction_ellipticity=True,
                is_not_global_or_algorithmic_tangent=True)
    return best


def scan_batch(data, normals):
    """Vectorized directions; each row remains an independently defined material branch."""
    if not data:
        return []
    bulk = np.array([d['bulk_kpa'] for d in data])
    shear = np.array([d['shear_kpa'] for d in data])
    den = np.array([d['denominator_kpa3'] for d in data])
    valid = np.isfinite(den) & (den > 0)
    indices = np.flatnonzero(valid)
    tensor = np.array([[[d['cn'][0], d['cn'][3], d['cn'][5]],
                        [d['cn'][3], d['cn'][1], d['cn'][4]],
                        [d['cn'][5], d['cn'][4], d['cn'][2]]] for d in data])
    best = np.full(len(data), np.inf)
    normal_index = np.full(len(data), -1, dtype=int)
    normals = np.asarray(normals)
    if normals.ndim != 2 or normals.shape[1] != 3 or not 1 <= len(normals) <= 32768:
        raise ValueError('bounded nonempty direction inventory required')
    for j, normal in enumerate(normals):
        if not np.isfinite(normal).all() or abs(np.linalg.norm(normal) - 1) > 1e-10:
            raise ValueError('a finite unit direction required')
        a = tensor[valid] @ normal
        q = shear[valid, None, None] * np.eye(3) + (bulk[valid] + shear[valid]/3)[:, None, None] * np.outer(normal, normal)
        q -= a[:, :, None] * a[:, None, :] / den[valid, None, None]
        value = np.linalg.eigvalsh(q)[:, 0]
        changed = value < best[valid]
        selected = indices[changed]
        best[selected] = value[changed]
        normal_index[selected] = j
    return [None if not valid[i] else dict(eigenvalue_kpa=float(best[i]),
             normalized_by_shear=float(best[i]/shear[i]), normal=normals[normal_index[i]].tolist())
            for i in range(len(data))]


def diagnose(mesh, ground, result, *, yield_tolerance=2e-6, normals=None):
    """Effective-stress skeleton only; near-yield points are potential plastic branches."""
    if not np.isfinite(yield_tolerance) or not 0 < yield_tolerance <= 1e-3:
        raise ValueError('finite positive dimensionless yield tolerance up to1e-3 required')
    if result.get('complete') is not True or not result.get('frames'):
        raise ValueError('complete result with committed frames required')
    layout = result.get('integration_points', [])
    owners = np.array([q['element_index'] for q in layout])
    weights = np.array([q['weight_m3'] for q in layout])
    count = len(mesh['tetrahedra'])
    if not len(layout) or owners.dtype.kind not in 'iu' or (owners < 0).any() or (owners >= count).any() or not np.isfinite(weights).all() or (weights <= 0).any():
        raise ValueError('positive integration weights and valid element ownership required')
    materials = {m['id']: m.get('mcc') for m in ground['materials']}
    if not all(isinstance(value, dict) for value in materials.values()):
        raise ValueError('MCC materials required for this diagnostic')
    normals = directions() if normals is None else np.asarray(normals)
    rows = []
    for frame in result['frames']:
        state = frame['state']
        if len(state['material_states']) != len(layout) or len(state['active_elements']) != count or any(type(x) is not bool for x in state['active_elements']):
            raise ValueError('committed state and integration inventory differ')
        active = np.asarray(state['active_elements'])[owners]
        candidates, ids = [], []
        for i in np.flatnonzero(active):
            value = state['material_states'][i]
            material = materials[mesh['material_ids'][owners[i]]]
            data = branch(value['effective_stress_compression_positive_kpa'], value['preconsolidation_pressure_kpa'], value['void_ratio'], material)
            if abs(data['normalized_yield']) <= yield_tolerance:
                candidates.append(data)
                ids.append(int(i))
        scanned = scan_batch(candidates, normals)
        negative = [k for k, value in enumerate(scanned) if value is not None and value['normalized_by_shear'] < -1e-8]
        valid_values = [value['normalized_by_shear'] for value in scanned if value is not None]
        row = dict(stage=frame['name'], active_gp=int(active.sum()), active_volume_m3=float(weights[active].sum()),
                   near_yield_candidates=len(ids), softening_candidates=int(sum(d['consistency_H_kpa3'] < 0 for d in candidates)),
                   invalid_plastic_consistency_denominators=sum(value is None for value in scanned),
                   sampled_negative_acoustic_candidates=len(negative),
                   sampled_negative_volume_m3=float(sum(weights[ids[k]] for k in negative)),
                   minimum_sampled_eigenvalue_over_G=min(valid_values) if valid_values else None)
        if negative:
            k = min(negative, key=lambda k: scanned[k]['normalized_by_shear'])
            i = ids[k]
            q = acoustic_tensor(candidates[k], scanned[k]['normal'], plastic=True)
            values, vectors = np.linalg.eigh(q)
            v = vectors[:, 0]
            row['witness'] = dict(global_gp=i, element=int(owners[i]), position_m=layout[i]['position_m'],
                                  normal=scanned[k]['normal'], polarization=v.tolist(), eigenvalue_kpa=float(values[0]),
                                  rayleigh_kpa=float(v @ q @ v), material_state=state['material_states'][i])
        rows.append(row)
    return dict(diagnostic_only=True, frames=rows, sampled_directions=len(normals),
                normalized_yield_tolerance=yield_tolerance, consistency_H_and_denominator_units='kPa^3',
                acoustic_eigenvalue_units='kPa', stress_order=['xx','yy','zz','xy','yz','xz'],
                actual_continuing_loading_branch_not_determined=True,
                sampled_positive_directions_do_not_prove_ellipticity=True,
                is_not_global_mixed_algorithmic_or_hydraulic_tangent=True,
                constitutive_law_and_fem_result_unchanged=True, spatial_convergence_accepted=False)


def main():
    import argparse,hashlib,json
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('mesh','ground','result','output'):
        parser.add_argument(name,type=Path)
    parser.add_argument('--yield-tolerance',type=float,default=2e-6)
    args=parser.parse_args()
    if args.output.exists():
        raise ValueError('new diagnostic output file required')
    paths=[args.mesh,args.ground,args.result]
    if any(path.stat().st_size > 1024**3 for path in paths):
        raise ValueError('offline diagnostic source exceeds1GiB file budget')
    data=[]
    before={}
    for path in paths:
        raw=path.read_bytes()
        before[str(path.resolve())]=hashlib.sha256(raw).hexdigest()
        data.append(json.loads(raw))
    report=diagnose(*data,yield_tolerance=args.yield_tolerance)
    if before != {str(path.resolve()):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}:
        raise ValueError('diagnostic inputs changed during computation')
    report['inputs_sha256']=before
    report['implementation_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    serialized=json.dumps(report,indent=2,allow_nan=False)
    with args.output.open('x',encoding='utf-8') as stream:
        stream.write(serialized)
    print(json.dumps(dict(diagnostic_only=True,frames=len(report['frames']),output=str(args.output))))


if __name__ == '__main__':
    main()
