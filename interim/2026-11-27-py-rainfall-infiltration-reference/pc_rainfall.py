"""Independent Cartesian 2D/3D VG rainfall FV for explicit user-PC runs.

3D is not run in the cloud. 2D is a strip analogue, not the square-patch case.
This module calls only this folder's VG constitutive function and SciPy/NumPy.
"""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from reference import vg


def solve(dim=2, n=8, dt=100., end=10800., user_pc=False):
    if dim not in (2, 3) or n < 4 or n % 4 or n > 64 or dt <= 0 or not 0 < end <= 10800:
        raise ValueError('dim=2/3, n=4..64 multiple of 4, positive dt, end<=10800 required')
    if dim == 3 and not user_pc:
        raise ValueError('3D requires explicit --user-pc; do not run a cloud series')
    started = time.perf_counter()
    shape, dx = (n,)*dim, 1/n
    coordinates = np.indices(shape)
    ids = np.arange(n**dim).reshape(shape)
    psi = np.broadcast_to((np.arange(n)+.5).reshape((n,)+(1,)*(dim-1))*dx-1, shape).copy()
    initial = vg(psi)[2]
    patch = np.ones(shape[1:], dtype=bool)
    for axis in range(dim-1):
        coordinate = (np.arange(n)+.5)*dx
        broadcast = [1]*(dim-1)
        broadcast[axis] = n
        patch &= ((coordinate >= .25) & (coordinate <= .75)).reshape(broadcast)
    # Cross stencil colouring: +/- axis weights are distinct mod (2*dim+1).
    colors = sum((axis+1)*coordinates[axis] for axis in range(dim)) % (2*dim+1)
    rows, cols = [ids.ravel()], [ids.ravel()]
    for axis in range(dim):
        left, right = [slice(None)]*dim, [slice(None)]*dim
        left[axis], right[axis] = slice(None, -1), slice(1, None)
        a, b = ids[tuple(left)].ravel(), ids[tuple(right)].ravel()
        rows.extend([a, b]); cols.extend([b, a])
    row, col = np.concatenate(rows), np.concatenate(cols)
    flat_colors = colors.ravel()
    masks = [flat_colors == color for color in range(2*dim+1)]
    entries = [flat_colors[col] == color for color in range(2*dim+1)]

    def fields(h, rain):
        _, k, storage = vg(h)
        divergence = np.zeros(shape)
        top = bottom = None
        for axis in range(dim):
            face_shape = list(shape); face_shape[axis] += 1
            flux = np.zeros(face_shape)
            left, right, interior = [slice(None)]*dim, [slice(None)]*dim, [slice(None)]*dim
            left[axis], right[axis], interior[axis] = slice(None, -1), slice(1, None), slice(1, -1)
            kl, kr = k[tuple(left)], k[tuple(right)]
            harmonic = 2*kl*kr/(kl+kr)
            flux[tuple(interior)] = harmonic*((1. if axis == 0 else 0.)-np.diff(h, axis=axis)/dx)
            if axis == 0:
                capacity = 2*k[0]*1e-6/(k[0]+1e-6)*(1-h[0]/(dx/2))
                top = np.minimum(rain*patch, capacity)
                bottom = 2*k[-1]*1e-6/(k[-1]+1e-6)*(1+h[-1]/(dx/2))
                flux[0], flux[-1] = top, bottom
            divergence += np.diff(flux, axis=axis)/dx
        return storage, divergence, top, bottom

    t, cumulative_top, cumulative_bottom, max_mass = 0., 0., 0., 0.
    frames, steps = [], []
    while t < end-1e-8:
        edge = min([a for a in (3600., 7200., 10800.) if a > t+1e-8]+[end])
        step = min(dt, edge-t, end-t)
        rain = 3e-6 if t < 3600.-1e-8 or t >= 7200.-1e-8 else 0.
        old = vg(psi)[2]

        def residual(flat):
            storage, divergence, _, _ = fields(flat.reshape(shape), rain)
            return ((storage-old)/step+divergence).ravel()

        h = psi.ravel().copy()
        for iteration in range(35):
            r = residual(h)
            if np.max(np.abs(r)) < 2e-15:
                break
            values = np.empty(len(row))
            epsilon = 1e-6
            for mask, selected in zip(masks, entries):
                plus, minus = h.copy(), h.copy()
                plus[mask] += epsilon; minus[mask] -= epsilon
                derivative = (residual(plus)-residual(minus))/(2*epsilon)
                values[selected] = derivative[row[selected]]
            jacobian = coo_matrix((values, (row, col)), shape=(h.size, h.size)).tocsc()
            change = spsolve(jacobian, -r)
            if not np.all(np.isfinite(change)):
                raise RuntimeError('nonfinite Newton update; not accepted')
            fraction = 1.
            while fraction >= 2**-16:
                trial = h+fraction*change
                if np.linalg.norm(residual(trial)) < np.linalg.norm(r):
                    h = trial
                    break
                fraction /= 2
            else:
                raise RuntimeError('Newton line search failed; not accepted')
        else:
            raise RuntimeError('Newton iteration limit; not accepted')
        psi = h.reshape(shape)
        storage, _, top, bottom = fields(psi, rain)
        cumulative_top += step*float(top.sum())*dx**(dim-1)
        cumulative_bottom += step*float(bottom.sum())*dx**(dim-1)
        mass = float((storage-initial).sum()*dx**dim)-cumulative_top+cumulative_bottom
        max_mass = max(max_mass, abs(mass))
        t += step
        steps.append(dict(time_s=t, dt_s=step, iterations=iteration,
                          mass_defect=mass, max_residual_per_s=float(np.max(np.abs(residual(h))))))
        if any(abs(t-a) < 1e-7 for a in (3600., 7200., 10800., end)):
            z = 1-(np.arange(n)+.5).reshape((n,)+(1,)*(dim-1))*dx
            frames.append(dict(time_s=t, head_m=(psi+z).tolist(),
                               storage_offset=float(storage.sum()*dx**dim),
                               cumulative_top=cumulative_top, cumulative_bottom=cumulative_bottom))
    return dict(dimension=dim, n=n, dt_s=dt, shape=list(shape),
                axis_order='depth_then_transverse_axes; array first index top to bottom',
                volume_unit='m3' if dim == 3 else 'm3 per m strip width',
                rain_patch_area=float(patch.sum()*dx**(dim-1)),
                original_square_geometry=dim == 3, complete=True, accepted=None,
                maximum_mass_defect=max_mass, frames=frames, steps=steps,
                elapsed_s=time.perf_counter()-started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dim', type=int, choices=(2, 3), default=2)
    parser.add_argument('--n', type=int, default=8)
    parser.add_argument('--dt', type=float, default=100.)
    parser.add_argument('--end', type=float, default=10800.)
    parser.add_argument('--user-pc', action='store_true')
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    path = Path(args.out)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() or path.with_suffix('.failure.json').exists():
        raise FileExistsError('use a fresh output path; previous evidence is preserved')
    try:
        report = solve(args.dim, args.n, args.dt, args.end, args.user_pc)
    except Exception as error:
        path.with_suffix('.failure.json').write_text(json.dumps(
            dict(complete=False, accepted=None, error=repr(error)), indent=2)+'\n', encoding='utf-8')
        raise
    path.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    print(json.dumps({k: report[k] for k in ('dimension', 'n', 'maximum_mass_defect', 'elapsed_s')}))
