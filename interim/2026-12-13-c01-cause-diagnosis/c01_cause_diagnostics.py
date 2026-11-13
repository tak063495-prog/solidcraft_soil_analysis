"""C-01 postprocessing and explicit small reference experiments; no production changes."""
import argparse
from collections import Counter
from decimal import Decimal, localcontext
import hashlib
import importlib.util
import json
from math import ceil, fsum, sqrt
from pathlib import Path
import shutil
import sys
import tempfile
import time
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PRIOR = ROOT / 'interim/2026-12-12-c01-switch-event-localization'
BASE = 'e8133f045a852c6ca374a552070526067f11dbe9'
STAGES = ('geostatic', 'shear', 'open', 'recontact', 'reverse', 'reload', 'return')
FIELDS = ('plastic_slip', 'plastic_dissipation', 'elastic_energy')


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def dependencies():
    # Same import arrangement as the retained studies; no source file is rewritten.
    sys.path.insert(0, str(ROOT / 'interim/2026-11-15-wet-joint-convergence'))
    old = module('cause_old_study', ROOT / 'interim/2026-11-15-wet-joint-convergence/study.py')
    ref = module('cause_reference', PRIOR / 'reference.py')
    return old, ref


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def write(path, value):
    with path.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')


def sha(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            result.update(block)
    return result.hexdigest()


def dot(a, b):
    return fsum(float(x) * float(y) for x, y in zip(a, b))


def norm(a):
    return sqrt(dot(a, a))


def energy(q, n, kn, kt):
    qn = dot(q, n)
    tau = [q[i] - qn * n[i] for i in range(3)]
    return qn * qn / (2 * kn) + dot(tau, tau) / (2 * kt)


def local_replay(start, qold, target, n, kn, kt, mu, count):
    """Independent scalar Coulomb return, with plane crossing and local endpoint work.

    Saved jump endpoints are inputs, never new global equilibrium solutions.
    Zero birth offset is restricted to the fixed C-01 specimen.
    """
    jump = list(start)
    q = list(qold)
    plastic = [0., 0., 0.]
    diss = work = numerical = 0.

    def closed(end):
        nonlocal jump, q, diss, work, numerical, plastic
        delta = [end[i] - jump[i] for i in range(3)]
        dn = dot(delta, n)
        qn0 = dot(q, n)
        trial = [q[i] - qn0 * n[i] + kt * (delta[i] - dn * n[i]) for i in range(3)]
        qn = kn * min(dot(end, n), 0.)
        limit = mu * -qn
        length = norm(trial)
        scale = min(1., limit / length) if length else 1.
        tau = [x * scale for x in trial]
        dp = [(trial[i] - tau[i]) / kt for i in range(3)]
        next_q = [qn * n[i] + tau[i] for i in range(3)]
        de = energy(next_q, n, kn, kt) - energy(q, n, kn, kt)
        dw = dot(next_q, delta)
        dd = max(dot(tau, dp), 0.)
        numerical += max(dw - de - dd, 0.)
        diss += dd
        work += dw
        plastic = [plastic[i] + dp[i] for i in range(3)]
        jump, q = list(end), next_q

    if list(start) == list(target):
        return q, plastic, diss, work, numerical
    for k in range(1, count + 1):
        end = list(target) if k == count else [start[i] + (target[i] - start[i]) * k / count for i in range(3)]
        g0, g1 = dot(jump, n), dot(end, n)
        if g0 <= 0 and g1 <= 0:
            closed(end)
        elif g0 > 0 and g1 > 0:
            jump, q = end, [0., 0., 0.]
        else:
            crossing = [jump[i] + (end[i] - jump[i]) * g0 / (g0 - g1) for i in range(3)]
            residual = dot(crossing, n)
            crossing = [crossing[i] - residual * n[i] for i in range(3)]
            if g0 <= 0:
                closed(crossing)
                jump, q = end, [0., 0., 0.]
            else:
                jump, q = crossing, [0., 0., 0.]
                closed(end)
    return q, plastic, diss, work, numerical


def metric(a, b, weights=None):
    """Fixed C-order traversal; 80-digit sqrt. No BLAS reduction or extrapolation."""
    import numpy as np
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('incompatible or nonfinite comparison')
    d = a - b
    weight = np.ones(a.shape) if weights is None else np.broadcast_to(weights, a.shape)
    ds = fsum(float(x) * float(x) * float(w) for x, w in zip(d.flat, weight.flat))
    rs = fsum(float(x) * float(x) * float(w) for x, w in zip(b.flat, weight.flat))
    with localcontext() as context:
        context.prec = 80
        relative = float(Decimal.from_float(ds).sqrt() / max(Decimal.from_float(rs).sqrt(), Decimal.from_float(1e-12)))
    return dict(difference_squared=ds, reference_squared=rs, pooled_relative_l2=relative)


def extract_selected(source, output):
    """Authenticate existing complete ZIP bytes, extract only diagnostic inputs/fields."""
    index_path = source / 'archive-index.json'
    index = read(index_path)
    acceptance = read(source / 'acceptance.json')
    if acceptance['archive_index_sha256'] != sha(index_path) or acceptance['passed'] is not False:
        raise ValueError('retained acceptance/index changed')
    output.mkdir(parents=True, exist_ok=False)
    selected = []
    with tempfile.TemporaryDirectory(dir=output.parent) as tmp:
        archive = Path(tmp) / 'evidence.zip'
        with archive.open('xb') as target:
            for part in index['parts']:
                path = source / part['file']
                if path.stat().st_size != part['bytes'] or sha(path) != part['sha256']:
                    raise ValueError('retained part changed')
                with path.open('rb') as stream:
                    shutil.copyfileobj(stream, target)
        if sha(archive) != index['sha256']:
            raise ValueError('retained ZIP changed')
        expected = {entry['path']: entry for entry in index['entries']}
        with zipfile.ZipFile(archive) as z:
            if len(z.namelist()) != len(expected) or set(z.namelist()) != set(expected):
                raise ValueError('retained ZIP inventory changed')
            for name in z.namelist():
                path = Path(name)
                if path.is_absolute() or '..' in path.parts:
                    raise ValueError('unsafe archive path')
                parts = path.parts
                fields = len(parts) >= 4 and parts[0] == 'paths' and parts[2] in ('native', 'independent')
                inputs = name in ('paths/mesh-h0.json', 'paths/ground.json', 'paths/stages.json')
                result = len(parts) == 3 and parts[0] == 'paths' and parts[2] == 'independent-result.json'
                if not (fields or inputs or result):
                    continue
                dest = output / path
                dest.parent.mkdir(parents=True, exist_ok=True)
                with z.open(name) as stream, dest.open('xb') as target:
                    shutil.copyfileobj(stream, target)
                if sha(dest) != expected[name]['sha256'] or dest.stat().st_size != expected[name]['bytes']:
                    raise ValueError('retained entry changed: ' + name)
                selected.append(name)
    write(output / 'selection.json', dict(index_sha256=sha(index_path), archive_sha256=index['sha256'], selected_entries=selected, all_selected_entry_sha256_verified=True))


def postprocess(folder, old, native_counts=None):
    import numpy as np
    mesh, ground, case = [read(folder / name) for name in ('mesh-h0.json', 'ground.json', 'stages.json')]
    from joint_reference import System
    system = System(mesh, case['joint_system'])
    layout_path = HERE/'geometry.json'
    if layout_path.exists():
        # Quadrature is part of the saved discretization. Reuse its sealed
        # doubles rather than recomputing eigenvalue-based Jacobi weights.
        layout = read(layout_path)['points']
        if len(layout) != len(system.points):
            raise ValueError('sealed quadrature point count differs')
        for actual, fixed in zip(system.points, layout):
            if actual['face_index'] != fixed['face_index'] or fixed['weight_m2'] <= 0:
                raise ValueError('sealed quadrature layout differs')
            actual['weight_m2'] = fixed['weight_m2']
            actual['normal'] = fixed['normal']
        system.weight = np.asarray([p['weight_m2'] for p in layout])
        system.normal = np.asarray([p['normal'] for p in layout])
    npoints = len(system.points)
    initial = dict(jump=np.zeros((npoints, 3)), traction=np.zeros((npoints, 3)), plastic_slip=np.zeros((npoints, 3)), plastic_dissipation=np.zeros(npoints))
    series, audits = {}, {}
    for kind in ('native', 'independent'):
        for steps in (16, 32, 64, 128, 256):
            key = f'{kind}/{steps}'
            path = folder / f'g{steps}' / kind
            saved_counts = None
            if kind == 'native' and native_counts is not None:
                captured = read(native_counts / f'n{steps}/counts.json')
                if captured['counts_bound_to_saved_old_and_new_displacement_bytes'] is not True or captured['accepted_coordinates_exact'] is not True:
                    raise ValueError('unbound native counts')
                saved_counts = {(r['stage'],r['progress']):r['contact_partition_counts'] for r in captured['rows']}
            if kind == 'independent':
                result = read(path.parent / 'independent-result.json')
                saved_counts = {(r['stage'], r['progress']): r['contact_partition_counts'] for r in result['path']}
            cumulative = {k: np.zeros(npoints) for k in ('replay', 'endpoint', 'trapezoid')}
            values, coordinates = [], []
            maxima = Counter()
            counts_hist = Counter()
            previous = initial
            negative_trapezoids = 0
            for row, fields in old.states(path):
                coord = (row['stage'], row['progress'])
                coordinates.append(coord)
                recalc_e = np.array([energy(q, n, kn, kt) for q, n, kn, kt in zip(fields['traction'], system.normal, system.kn, system.kt)])
                balance = fields['endpoint_work'] - recalc_e - fields['numerical_dissipation']
                maxima['elastic_energy'] = max(maxima['elastic_energy'], float(np.max(abs(recalc_e-fields['elastic_energy']))))
                maxima['balance_dissipation'] = max(maxima['balance_dissipation'], float(np.max(abs(balance-fields['plastic_dissipation']))))
                for i in range(npoints):
                    n, kn, kt, mu = system.normal[i], system.kn[i], system.kt[i], system.mu[i]
                    q0, q1 = previous['traction'][i], fields['traction'][i]
                    dp = fields['plastic_slip'][i] - previous['plastic_slip'][i]
                    tau0 = q0 - dot(q0, n) * n
                    tau1 = q1 - dot(q1, n) * n
                    ep = dot(tau1, dp)
                    trap = dot((tau0+tau1)*.5, dp)
                    cumulative['endpoint'][i] += ep
                    cumulative['trapezoid'][i] += trap
                    negative_trapezoids += int(trap < 0.)
                    count = saved_counts[coord][i] if saved_counts is not None else max(1, ceil(norm(fields['jump'][i]-previous['jump'][i])/5e-5))
                    counts_hist[count] += 1
                    q, plastic, dd, _, _ = local_replay(previous['jump'][i], q0, fields['jump'][i], n, kn, kt, mu, count)
                    cumulative['replay'][i] += dd
                    maxima['replay_traction'] = max(maxima['replay_traction'], max(abs(q[j]-q1[j]) for j in range(3)))
                    maxima['replay_plastic_increment'] = max(maxima['replay_plastic_increment'], max(abs(plastic[j]-dp[j]) for j in range(3)))
                maxima['replay_dissipation'] = max(maxima['replay_dissipation'], float(np.max(abs(cumulative['replay']-fields['plastic_dissipation']))))
                values.append(dict(plastic_slip=fields['plastic_slip'], plastic_dissipation=fields['plastic_dissipation'], elastic_energy=fields['elastic_energy'], balance=balance, **{k:v.copy() for k,v in cumulative.items()}))
                previous = fields
            if len(coordinates) != len(set(coordinates)):
                raise ValueError('duplicate actual coordinate')
            series[key] = (dict(zip(coordinates, range(len(coordinates)))), {k:np.stack([v[k] for v in values]) for k in values[0]})
            audits[key] = dict(states=len(coordinates), point_increment_samples=len(coordinates)*npoints, partition_counts_source=('saved independent path' if kind == 'independent' else 'observed native; old/new accepted u bytes matched') if saved_counts is not None else 'lower bound ceil(saved jump change / maximum); native frozen counts not archived for all phases', partition_histogram={str(k):v for k,v in sorted(counts_hist.items())}, maximum_absolute_recalculation_errors=dict(maxima), negative_global_trapezoid_increment_samples=negative_trapezoids)
    pairs = {}
    for kind in ('native', 'independent'):
        for a, b in zip((16, 32, 64, 128), (32, 64, 128, 256)):
            coords = [(s, i/a) for s in STAGES for i in range(1, a+1)]
            ia, fa = series[f'{kind}/{a}']; ib, fb = series[f'{kind}/{b}']
            ids_a, ids_b = [ia[c] for c in coords], [ib[c] for c in coords]
            fields = {}
            for key in fa:
                x, y = fa[key][ids_a], fb[key][ids_b]
                spatial_weight = system.weight[:, None] if key == 'plastic_slip' else system.weight
                tx, ty, tw = [], [], []
                for j, stage in enumerate(STAGES):
                    for i in range(a+1):
                        if i == 0 and j == 0:
                            vx, vy = np.zeros_like(x[0]), np.zeros_like(y[0])
                        else:
                            coord = (STAGES[j-1],1.) if i == 0 else (stage,i/a)
                            vx, vy = fa[key][ia[coord]], fb[key][ib[coord]]
                        tx.append(vx); ty.append(vy)
                        tw.append((.5 if i in (0,a) else 1.)/a)
                time_weight = np.asarray(tw).reshape((-1,)+(1,)*(x.ndim-1))*spatial_weight
                fields[key] = dict(original_unweighted=metric(x, y), area_weighted=metric(x, y, spatial_weight), progress_trapezoid_area_weighted=metric(np.stack(tx),np.stack(ty),time_weight))
            pairs[f'{kind}/{a}__{b}'] = fields
    return dict(audits=audits, pairs=pairs), system


def geometry(system, mesh):
    import numpy as np
    points = []
    for i, p in enumerate(system.points):
        mismatch = max(norm(np.asarray(mesh['nodes'][a])-np.asarray(mesh['nodes'][b])) for a, b in zip(p['nodes_a'], p['nodes_b']))
        points.append(dict(index=i, face_index=p['face_index'], position_m=p['position_m'], normal=p['normal'], weight_m2=p['weight_m2'], maximum_paired_node_mismatch_m=mismatch))
    return dict(contact_gp_count=len(points), face_count=len({p['face_index'] for p in system.points}), area_m2=fsum(float(w) for w in system.weight), maximum_paired_node_mismatch_m=max(p['maximum_paired_node_mismatch_m'] for p in points), geometry_kind='matched quadratic paired faces; fixed reference normals/global vector history; no facet search or moving tangent basis', points=points)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--extract', choices=('switch', 'path'))
    p.add_argument('--input', type=Path)
    p.add_argument('--native-counts', type=Path)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if args.extract:
        source = ROOT / ('interim/2026-12-11-c01-global-step-control' if args.extract == 'switch' else 'interim/2026-12-12-c01-switch-event-localization')
        extract_selected(source, args.output)
        return
    if args.input is None:
        p.error('--input required for postprocessing')
    args.output.mkdir(parents=True, exist_ok=False)
    old, _ = dependencies()
    started = time.monotonic()
    report, system = postprocess(args.input, old, args.native_counts)
    report.update(schema='c01-cause-measurement/1', starting_commit=BASE, passed=False, c01_decision='未達', numerical_solve=False, original_acceptance_retained=True, extrapolation_used=False, elapsed_seconds=time.monotonic()-started)
    write(args.output / 'measurement.json', report)
    write(args.output / 'geometry.json', geometry(system, read(args.input / 'mesh-h0.json')))


if __name__ == '__main__':
    main()
