"""Small independent references; no production invocation. Exit 0 means generated."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from reference import cases
from numerical import integrate_forces, series, solve_head, fitted_mesh
from compare import synthetic, diagnose, write_csv_bundle
from policies import reservoir_case, trace


def plain(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):
        return obj.item()
    raise TypeError(type(obj).__name__)


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, default=plain, allow_nan=False)+'\n', encoding='utf-8')


def run(output, pc=False):
    start = time.monotonic()
    report = dict(schema='independent-wet-soil-birth/1', a03_complete=False, accepted=None,
                  production_solver_used=False, production_comparison='未実施', cases={}, series={},
                  tolerance='未定義', pc_series=pc)
    for name, (column, flow) in cases().items():
        reference = integrate_forces(column, flow)
        z = np.unique(np.r_[np.linspace(0, column.top, 31), column.interface,
                           min(column.top, column.water), *(flow.edges if flow else [])])
        p = column.pressure(z) if flow is None else flow.pressure(z)
        fields = column.fields(z, p)
        write(output/name/'profiles.json', dict(z_m=z, **fields))
        write(output/name/'all-nodal-forces.json', reference)
        metrics = dict(reference['metrics'])
        if flow:
            zz = fitted_mesh(flow)
            solution = solve_head(flow, zz)
            numerical = integrate_forces(column, flow, numerical_result=solution)
            metrics.update(analytic_flux_m_s=flow.flux, numerical_water_defect=solution['water_defect'],
                           numerical_force_metrics=numerical['metrics'])
            write(output/name/'numerical-head.json', solution)
            write(output/name/'numerical-nodal-forces.json', numerical)
        report['cases'][name] = metrics
    ns = (16, 32, 64, 128, 256, 512) if pc else (8, 16, 32, 64)
    for name in ('darcy_upward', 'darcy_downward'):
        column, flow = cases()[name]
        report['series'][name] = series(flow, ns, column)
    # Space and quadrature series are separate. Kinks are split exactly; Q4
    # force identities should be exact to roundoff at each mesh, not fitted.
    c, f = cases()['partial_new_layer']
    report['force_space_series'] = [dict(subdivisions=n, **integrate_forces(c, f, subdivisions=n)['metrics'])
                                    for n in (2, 4, 8, 16)]
    report['force_quadrature_series'] = [dict(order=q, subdivisions=4,
        **integrate_forces(c, f, quadrature=q)['metrics']) for q in (2, 3, 4)]
    report['force_series_observed_order'] = None
    report['force_series_note'] = 'exact polynomial integrals after splitting: roundoff, no inferred convergence order'
    reservoir = reservoir_case()
    for key in ('low', 'high'):
        problem = reservoir.pop(key)
        solution = solve_head(problem, fitted_mesh(problem))
        reservoir[key+'_numerical_outward'] = solution['outward']
        reservoir[key+'_numerical_water_defect'] = solution['water_defect']
    reservoir['numerical_reservoir_out'] = float(reservoir['low_numerical_outward'][1]+reservoir['high_numerical_outward'][0])
    reservoir['numerical_water_defect'] = float(reservoir['low_numerical_outward'][0]
        +reservoir['high_numerical_outward'][1]+reservoir['numerical_reservoir_out'])
    report['prescribed_internal_reservoir'] = reservoir
    report['policy_traces'] = {kind: trace(kind, 1.2, [0., 1.], head_a=4.,
        head_b=4. if kind != 'soil' else None, soil_side='a' if kind == 'soil' else None,
        level=4. if kind == 'prescribed' else None) for kind in ('soil', 'continuous', 'prescribed')}
    sample, before = synthetic('partial_new_layer')
    write(output/'synthetic-only/export.json', sample)
    write(output/'synthetic-only/history-before.json', before)
    write_csv_bundle(sample, output/'synthetic-only/csv')
    write(output/'synthetic-only/comparison.json', diagnose(sample, before, 'partial_new_layer'))
    report['elapsed_seconds'] = time.monotonic()-start
    write(output/'measurements.json', report)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--pc-series', action='store_true')
    args = parser.parse_args()
    run(args.output, args.pc_series)
