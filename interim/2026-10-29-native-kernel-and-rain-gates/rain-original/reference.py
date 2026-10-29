"""Input-only, cell-centred 3D transient Richards finite volumes.

Unit cube, H=0 at base, impermeable sides, no-ponding atmospheric top.
The optional rain patch is [.25,.75] x [.25,.75]. No Rust output is read.
Seven-colour numerical Newton + independently checked BiCGSTAB, harmonic K.
This is a homogeneous/isotropic hydraulic benchmark, not a mechanical model.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
from pathlib import Path
import copy
import json
import sys
import time
import numpy as np

WORKSPACE = Path(__file__).resolve().parents[2]
ROOT = Path(__file__).resolve().parent / 'runtime'
if not ROOT.is_dir():
    ROOT = WORKSPACE
sys.path.insert(0, str(ROOT / 'tools'))
from richards_column_reference import retention
from richards_dam_reference import bicgstab


class Grid:
    def __init__(self, case, n, local=True):
        assert 2 <= n <= 64 and (not local or n % 4 == 0)
        self.n = n
        self.h = 1. / n
        self.z = (np.arange(n) + .5) / n
        self.height = self.z[:, None, None]
        mat = case['materials'][0]
        self.law = copy.deepcopy(mat['retention'])
        self.ks = float(mat['saturated_conductivity_m_s'][0][0])
        assert self.ks > 0
        for m in case['materials']:
            assert m['retention'] == self.law
            np.testing.assert_array_equal(m['saturated_conductivity_m_s'], np.eye(3) * self.ks)
            assert m.get('volume_source_per_s', 0) == 0
        assert case['initial_head'] == {'head_m': 0.}
        assert case['analysis']['kind'] == 'transient'
        a = (self.z > .25) & (self.z < .75)
        self.patch = a[:, None] & a[None, :] if local else np.ones((n, n), bool)
        z, y, x = np.indices((n, n, n))
        self.colors = ((x + 2*y + 3*z) % 7).ravel()
        indices = np.arange(n**3).reshape(n, n, n)
        self.columns = np.tile(indices.ravel()[:, None], (1, 7))
        self.valid = np.zeros_like(self.columns, dtype=bool)
        self.valid[:, 0] = True
        for j, (dz, dy, dx) in enumerate([(0,0,-1),(0,0,1),(0,-1,0),(0,1,0),(-1,0,0),(1,0,0)], 1):
            mask = (z+dz>=0)&(z+dz<n)&(y+dy>=0)&(y+dy<n)&(x+dx>=0)&(x+dx<n)
            self.valid[:, j] = mask.ravel()
            self.columns[mask.ravel(), j] = indices[(z+dz)[mask], (y+dy)[mask], (x+dx)[mask]]

    def field(self, head, rain):
        theta, kr, storage = retention(head - self.height, self.law)
        k = self.ks * kr
        harmonic = lambda a, b: 2*a*b/(a+b)
        n, h = self.n, self.h
        fx = np.zeros((n, n, n+1))
        fy = np.zeros((n, n+1, n))
        fz = np.zeros((n+1, n, n))
        fx[:, :, 1:-1] = -harmonic(k[:, :, :-1], k[:, :, 1:]) * np.diff(head, axis=2) / h
        fy[:, 1:-1, :] = -harmonic(k[:, :-1, :], k[:, 1:, :]) * np.diff(head, axis=1) / h
        fz[1:-1, :, :] = -harmonic(k[:-1, :, :], k[1:, :, :]) * np.diff(head, axis=0) / h
        # Bottom pressure head is zero; total head H=0. Positive z flux injects.
        fz[0] = -harmonic(self.ks, k[0]) * head[0] / (h/2)
        capacity = -harmonic(self.ks, k[-1]) * (1. - head[-1]) / (h/2)
        fz[-1] = np.maximum(-rain * self.patch, capacity)
        outward = (np.diff(fx, axis=2) + np.diff(fy, axis=1) + np.diff(fz, axis=0)) * h*h
        return theta, storage, outward, fx, fy, fz

    def residual(self, head, old_storage, dt, rain):
        _, storage, outward, *_ = self.field(head, rain)
        return (storage - old_storage) * self.h**3 / dt + outward

    def jacobian(self, head, old_storage, dt, rain):
        coefficients = np.zeros_like(self.columns, dtype=float)
        delta = 1e-6 * (1 + np.abs(head.ravel()))
        for color in range(7):
            mask = self.colors == color
            plus, minus = head.copy(), head.copy()
            plus.ravel()[mask] += delta[mask]
            minus.ravel()[mask] -= delta[mask]
            difference = (self.residual(plus, old_storage, dt, rain) - self.residual(minus, old_storage, dt, rain)).ravel()
            for j in range(7):
                use = self.valid[:, j] & (self.colors[self.columns[:, j]] == color)
                coefficients[use, j] = difference[use] / (2 * delta[self.columns[use, j]])
        return coefficients


def solve(case, n, dt_max, local=True, stop=None):
    assert np.isfinite(dt_max) and dt_max > 0
    grid = Grid(case, n, local)
    h = np.zeros((n, n, n))
    initial = grid.field(h, 0.)[1].sum() * grid.h**3
    clock = 0.
    boundary = rain_total = routed = 0.
    frames, steps = [], []
    maximum_mass_error = 0.
    started = time.monotonic()
    for period in case['analysis']['periods']:
        bc = {b['face_set']: b['condition'] for b in period['boundaries']}
        assert bc['bottom'] == {'kind': 'head', 'value': {'head_m': 0.}}
        assert bc['top']['kind'] == 'atmospheric_no_ponding'
        rain = float(bc['top']['rainfall_vertical_m_s'])
        end = float(period['end_time_s'])
        assert end > clock and rain >= 0
        while clock < end:
            if stop and stop.exists():
                raise InterruptedError('cooperative reference stop before next time step')
            dt = min(dt_max, end - clock)
            old = grid.field(h, rain)[1]
            trial = h.copy()
            linear_max = 0.
            for it in range(80):
                residual = grid.residual(trial, old, dt, rain)
                norm = float(np.max(np.abs(residual)))
                if norm < 1e-17:
                    break
                coeff = grid.jacobian(trial, old, dt, rain)
                increment, _, true = bicgstab(coeff, grid.columns, -residual.ravel(), rtol=2e-11)
                linear_max = max(linear_max, true)
                for half in range(24):
                    candidate = trial + 2.**(-half) * increment.reshape(h.shape)
                    new = grid.residual(candidate, old, dt, rain)
                    if np.isfinite(candidate).all() and np.max(np.abs(new)) < norm * (1-1e-4*2.**(-half)):
                        trial = candidate
                        break
                else:
                    raise AssertionError(('FV line search failed', clock, dt, norm))
            else:
                raise AssertionError(('FV Newton budget exhausted', clock, dt))
            h = trial
            theta, storage, _, fx, fy, fz = grid.field(h, rain)
            clock += dt
            boundary += float(fz[-1].sum()-fz[0].sum()) * grid.h**2 * dt
            rain_total += rain * grid.patch.sum() * grid.h**2 * dt
            routed += float((rain*grid.patch+fz[-1]).sum()) * grid.h**2 * dt
            error = abs(float(storage.sum()*grid.h**3) - initial + boundary)
            maximum_mass_error = max(maximum_mass_error, error)
            assert error < 1e-9, ('integrated FV mass error', error)
            steps.append(dict(time_s=clock, dt_s=dt, newton_iterations=it,
                maximum_cell_residual_m3_s=norm, maximum_linear_true_residual_l2=linear_max, integrated_mass_error_m3=error))
        transverse = float(np.sqrt((np.mean(fx[:,:,1:-1]**2)+np.mean(fy[:,1:-1,:]**2))/2))
        frames.append(dict(name=period['name'], time_s=clock, head_m=h.tolist(), water_content=theta.tolist(),
            storage_offset_m3=float(storage.sum()*grid.h**3), cumulative_boundary_outflow_m3=boundary,
            cumulative_precipitation_m3=float(rain_total), cumulative_surface_routed_water_m3=routed,
            top_outward_flux_m3_s=float(fz[-1].sum()*grid.h**2), bottom_upward_flux_m3_s=float(fz[0].sum()*grid.h**2),
            transverse_flux_rms_m_s=transverse))
    return dict(schema='soil-fem-independent-local-rain-3d-fv/1', cells_per_axis=n,
        maximum_time_step_s=dt_max, rain_patch_area_m2=float(grid.patch.sum()*grid.h**2),
        maximum_integrated_mass_error_m3=maximum_mass_error, elapsed_seconds=time.monotonic()-started,
        steps=steps, frames=frames, reads_native_results=False, mechanics_not_included=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('case', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--n', type=int, required=True)
    parser.add_argument('--dt', type=float, required=True)
    parser.add_argument('--full-width', action='store_true')
    args = parser.parse_args()
    value = solve(json.loads(args.case.read_text('utf-8-sig')), args.n, args.dt, not args.full_width, Path(__file__).parent/'pause-requested')
    with args.output.open('x', encoding='utf-8') as f:
        json.dump(value, f, allow_nan=False)
    print('PASS independent3D FV', args.n, args.dt, value['maximum_integrated_mass_error_m3'], flush=True)
