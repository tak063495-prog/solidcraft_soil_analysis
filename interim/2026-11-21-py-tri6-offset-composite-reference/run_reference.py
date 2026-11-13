"""Run only the independent small reference; no production solver invocation."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import time
import numpy as np
import scipy
import reference as ref
import compare


def packet(result, case='eccentric-cantilever'):
    report, xy, cells, d, reaction = result
    m = compare.metadata(report['e_m'], n=report['n'])
    nu = .3 if case == 'winkler' else 0.
    if case == 'winkler':
        m.update(case=case, nu=nu, load=dict(q_kpa=1., kf_kn_m3=100., shape='sin(pi*x/2)*sin(pi*y)'),
                 boundary='w=0;beta_t=0;M_nn=0;u_p=v_p=rz=0 on all nodes for membrane',
                 geometry='L=2m,H=1m,straight-planar-tri6,nx=2n,ny=n')
    if case == 'anchor':
        m.update(case=case, load=dict(Pz_kn=1., anchor_k_kn_m=50., truss_EA_kn=100., truss_length_m=2.,
                                    attachment='uniform-unit-width-tip-distributor;bilateral'))
    nodes = [dict(id=str(i), **dict(zip(compare.NODE_FIELDS,
             map(float, [*x, 0., *u, *r])))) for i, (x, u, r) in
             enumerate(zip(xy, d.reshape(-1, 6), reaction.reshape(-1, 6)))]
    sections = [dict(id=str(i), **dict(zip(compare.SECTION_FIELDS, map(float,
                [1/3, 1/3, *ref.section(xy, cell, d, 2.e5, nu, .2, report['e_m'])]))))
                for i, cell in enumerate(cells)]
    return compare.validate(dict(metadata=m, nodes=nodes, sections=sections))


def orders(rows, key):
    for previous, row in zip(rows, rows[1:]):
        a, b = previous[key], row[key]
        row['observed_order'] = None if min(a, b) <= 0 else float(np.log(a/b)/np.log(row['n']/previous['n']))
    rows[0]['observed_order'] = None


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--h-series', nargs='+', type=int, default=[2, 4, 8])
    args = p.parse_args()
    if len(args.h_series) < 3 or any(n < 1 for n in args.h_series) or any(b <= a for a, b in zip(args.h_series, args.h_series[1:])):
        p.error('at least three strictly increasing positive h subdivisions required')
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    reports = dict(accepted=None, B_02='未達', tolerance='未定義',
                   runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
                   source_sha256={f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in Path(__file__).parent.glob('*.py')})
    reports['rigid_modes'] = [ref.rigid_diagnostics(e) for e in [-.1, 0., .1]]
    xy, cells = ref.mesh(1, 1); element_xy = xy[cells[0]]
    rng = np.random.default_rng(713)
    reports['independent_energy'] = []
    for e in [-.1, 0., .1]:
        d = rng.normal(size=36)*1.e-4
        K, _ = ref.element(element_xy, 2.e5, .3, .2, e)
        direct = ref.direct_energy(element_xy, d, 2.e5, .3, .2, e)
        reports['independent_energy'].append(dict(e_m=e, direct_thickness_kn_m=direct,
               quadratic_kn_m=.5*float(d@K@d), absolute_error_kn_m=abs(direct-.5*float(d@K@d))))
    reports['eccentric_strips'] = []
    for support in ['cantilever', 'simple']:
        for e in [-.1, 0., .1]:
            for n in [1, 2, 4]:
                for N, M in [(3., 0.), (0., .2), (3., .2)]:
                    reports['eccentric_strips'].append(ref.strip_run(n, e, N, M, support)[0])
    reports['winkler_h'] = []
    for n in args.h_series:
        result = ref.winkler_run(n); report = result[0]; reports['winkler_h'].append(report)
        (args.output/f'winkler-n{n}.json').write_text(json.dumps(packet(result, 'winkler'), indent=2, allow_nan=False)+'\n', encoding='utf-8')
        print('Winkler', n, report['w_relative_L2'], flush=True)
    orders(reports['winkler_h'], 'w_relative_L2')
    reports['anchor_h'] = []
    for n in args.h_series:
        result = ref.strip_run(n, N=0., anchor=50.); report = result[0]
        (args.output/f'anchor-n{n}.json').write_text(json.dumps(packet(result, 'anchor'), indent=2, allow_nan=False)+'\n', encoding='utf-8')
        report['tip_absolute_error_m'] = report['errors']['w_tip_absolute_m']
        reports['anchor_h'].append(report)
    orders(reports['anchor_h'], 'tip_absolute_error_m')
    # Integration sequence is separate from h, using a fixed mesh.
    reports['quadrature_fixed_h'] = [dict(order=q, **ref.winkler_run(2, order=q)[0]) for q in [4, 5, 6]]
    reports['water'] = []
    xy = np.array([[0., 0.], [1., 0.], [0., 1.], [.5, 0.], [.5, .5], [0., .5]])
    for level in [-.1, 0., .4, 1., 1.3]:
        for e in [-.1, 0., .1]:
            f = ref.pressure_nodes(xy, level, e=e)
            closed = ref.pressure_closed(level)
            # Arbitrary quadratic virtual w, independent monomial nodal force.
            d = np.zeros(36); d[2::6] = .1+xy[:, 0]**2-.3*xy[:, 0]*xy[:, 1]+.2*xy[:, 1]**2
            reports['water'].append(dict(level_m=level, e_m=e, nodal_force_kn=f[2::6].tolist(),
                  closed_nodal_force_kn=closed.tolist(), nodal_force_max_error_kn=float(np.max(abs(f[2::6]-closed))),
                  force_sum_kn=float(np.sum(f[2::6])),
                  moment_sum_kn_m=np.sum(np.cross(np.column_stack((xy, np.zeros(6))), f.reshape(-1, 6)[:, :3])+f.reshape(-1, 6)[:, 3:], axis=0).tolist(),
                  virtual_work_error_kn_m=abs(float(d@f-d[2::6]@closed))))
    reports['mcc_defined_invariants_only'] = ref.mcc_invariants(np.diag([50., 70., 100.]), 150., 1.2)
    reports['defined_stage_algebra_only'] = dict(old_external_kn=[10., -5., 6.], new_external_kn=[2., -1., 3.],
        saved_removed_soil_internal_kn=[3., 1., -2.],
        death_load_and_holding_kn={str(x): [a.tolist() for a in ref.death_load([10., -5., 6.], [2., -1., 3.], [3., 1., -2.], x)] for x in [0., .25, 1.]},
        coupled_soil_birth_mcc='仕様未定義', coupled_plastic_history='未達', actual_restart='未実行', accepted=None)
    fixture = packet(ref.strip_run(2))
    (args.output/'reference-output.json').write_text(json.dumps(fixture, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    compare.write_csv(fixture, args.output/'reference-csv')
    reports['synthetic_roundtrip'] = compare.compare(fixture, compare.read(args.output/'reference-csv'))
    reports['elapsed_seconds'] = time.monotonic()-started
    (args.output/'measurements.json').write_text(json.dumps(reports, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')
    print('diagnostic completed; B-02 acceptance remains undefined', flush=True)


if __name__ == '__main__':
    main()
