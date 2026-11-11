"""Linear saturated limit of existing Richards, against an independent series.

Extension of PR #19's linear diffusion eigenfunction reference to equal,
piecewise-constant head changes on both ends of a unit-height column. All
coefficients, switches and units come from the input. This control separates
algebra agreement on RT0 from error relative to a continuum PDE solution.
"""
import argparse
import copy
import importlib.util
from pathlib import Path
import time
import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('c04_saturated_study', HERE/'study.py')
s = importlib.util.module_from_spec(spec); spec.loader.exec_module(s)
spec = importlib.util.spec_from_file_location('c04_saturated_measure', HERE/'measure.py')
measure = importlib.util.module_from_spec(spec); spec.loader.exec_module(measure)


class Series:
    def __init__(self, case):
        mats = case['materials']; first = mats[0]
        self.k = first['saturated_conductivity_m_s'][0][0]
        self.ss = first['retention']['saturated_specific_storage_per_m']
        if self.ss <= 0 or case['initial_head'].get('gradient', [0, 0, 0]) != [0, 0, 0]:
            raise ValueError('constant initial head and positive saturated storage required')
        for mat in mats:
            if (mat['retention'] != first['retention'] or
                    not np.array_equal(mat['saturated_conductivity_m_s'], np.eye(3)*self.k) or
                    mat.get('volume_source_per_s', 0) != 0):
                raise ValueError('homogeneous isotropic source-free series only')
        self.initial = case['initial_head']['head_m']; prior = self.initial; start = 0.; self.switches = []
        for period in case['analysis']['periods']:
            boundaries = {x['face_set']: x['condition'] for x in period['boundaries']}
            if set(boundaries) != {'bottom', 'top'} or any(x['kind'] != 'head' for x in boundaries.values()):
                raise ValueError('equal top and bottom head boundaries required')
            values = [x['value'] for x in boundaries.values()]
            if values[0] != values[1] or values[0].get('gradient', [0, 0, 0]) != [0, 0, 0]:
                raise ValueError('equal spatially constant head boundaries required')
            head = values[0]['head_m']
            if head <= 1 or self.initial <= 1: raise ValueError('strict saturation over unit height required')
            if head != prior: self.switches.append((start, head-prior))
            prior = head; start = period['end_time_s']
        self.diffusivity = self.k/self.ss

    def fields(self, z, t, terms=256):
        z = np.asarray(z, dtype=float)
        if np.any(z < 0) or np.any(z > 1) or t < 0: raise ValueError('unit-height domain and nonnegative time')
        n = 2*np.arange(terms)+1; wave = n*np.pi
        head = np.full(z.shape, self.initial); flux = np.zeros(z.shape)
        mean_delta = 0.
        for switch, delta in self.switches:
            tau = t-switch
            if tau <= 0: continue  # Input periods use their old boundary at the previous end frame.
            decay = np.exp(-self.diffusivity*wave**2*tau)
            head += delta*(1-(np.sin(z[..., None]*wave)*(4/wave*decay)).sum(axis=-1))
            flux += 4*self.k*delta*(np.cos(z[..., None]*wave)*decay).sum(axis=-1)
            mean_delta += delta*(1-np.sum(8/wave**2*decay))
        return head, flux, self.ss*mean_delta


def run(out, steps=(4, 8, 16), timeout=900):
    folder = out/'saturated-controls'; folder.mkdir(exist_ok=False)
    mesh = s.read(out/'hydraulic-mesh.json'); reports = {}; audit = {}; processes = []; start = time.monotonic()
    for n in steps:
        path = folder/f'h{n}'; path.mkdir()
        case = s.hydraulic_case(n); case['initial_head'] = dict(head_m=2.)
        for p in case['analysis']['periods']:
            phase = int(p['name'].split('-')[0][1:])
            value = (2., 2.2, 1.8, 2.1, 1.9, 2.15, 2.)[phase]
            for b in p['boundaries']: b['condition']['value']['head_m'] = value
        s.write(path/'case.json', case)
        proc = s.invoke([s.binary(out, 'solver'), 'solve-richards', out/'hydraulic-mesh.json', path/'case.json', path/'native.json'], path/'native.log', timeout)
        processes.append(proc)
        assert proc['exit_code'] == 0
        native = s.read(path/'native.json'); assert native['complete'] and all(a['accepted'] for a in native['attempts'])
        independent = s.hydraulics.Model(mesh, case).solve(); s.write(path/'independent.json', independent)
        series = Series(case); results = {}; rows = []; max_tail = 0.; max_mass = 0.
        native_times = list(measure.hydraulic_fields(out, path, 'native'))
        for frame, (t, discrete) in zip(native['frames'][1:], native_times):
            points = frame['points']; z = np.array([p['position_m'][2] for p in points])
            exact, flux, mass = series.fields(z, t)
            doubled = series.fields(z, t, terms=512)[0]
            max_tail = max(max_tail, float(np.max(np.abs(doubled-exact))))
            measure.add(results, 'hydraulic_gp_total_head_m', np.array([p['total_head_m'] for p in points]), exact)
            exact_flux = np.column_stack((np.zeros(len(z)), np.zeros(len(z)), flux))
            measure.add(results, 'hydraulic_gp_darcy_flux_m_s', np.array([p['darcy_flux_m_s'] for p in points]), exact_flux)
            cell_z = np.asarray(mesh['nodes'])[mesh['tetrahedra']].mean(axis=1)[:, 2]
            measure.add(results, 'cell_centroid_head_m', discrete['cell_head_m'], series.fields(cell_z, t)[0])
            face_z = np.array([f['reference_elevation_m'] for f in native['faces']])
            measure.add(results, 'face_centroid_head_m', discrete['face_head_m'], series.fields(face_z, t)[0])
            node_z = np.array(mesh['nodes'])[:, 2]
            measure.add(results, 'diagnostic_nodal_head_projection_m', discrete['nodal_head_projection_m'], series.fields(node_z, t)[0])
            native_delta = sum(frame['state']['cell_storage_offset_m3'])-sum(native['frames'][0]['state']['cell_storage_offset_m3'])
            max_mass = max(max_mass, abs(native_delta-mass))
            rows.append(dict(time_s=t, native_storage_change_m3=native_delta, exact_storage_change_m3=mass, storage_difference_m3=native_delta-mass))
        assert len(rows) == 7*n
        reports[str(n)] = dict(fields=measure.finish(results), maximum_series_256_512_head_difference_m=max_tail,
                               maximum_storage_difference_m3=max_mass, all_physical_times=rows,
                               continuum_PDE_acceptance=False,
                               note='RT0 hydraulic GP head is cellwise constant; cell/face comparisons sample the series at explicitly labelled centroids. Nodal incident-face mean is an observation only.')
        audit[str(n)] = measure.hydraulic_audit(path)
    sensitivity = {}
    for a, b in zip(steps, steps[1:]):
        for kind in ('native', 'independent'):
            sensitivity[f'{kind}/{a}__{b}'] = measure.hydraulic_pair(out, folder/f'h{a}', folder/f'h{b}', kind, kind)
    report = dict(complete=True, input_only_series=True, production_modified=False, cases=reports, mass_audits=audit,
                  time_sensitivity=sensitivity, processes=processes, elapsed_seconds=time.monotonic()-start,
                  scope='Saturated linear limit, K=.02 m/s, Ss=.02 /m, unit cube, equal top/bottom step heads. Unsaturated VG case remains separate.',
                  C04='未達', passed=False, A04_complete=False)
    s.write(folder/'measurements.json', report)
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--steps', nargs='+', type=int, default=[4, 8, 16]); p.add_argument('--timeout', type=int, default=900)
    a = p.parse_args(); run(a.output.resolve(), a.steps, a.timeout)
