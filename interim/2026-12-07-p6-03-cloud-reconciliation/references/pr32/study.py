"""Small, bounded study and optional larger 1D series for the user's PC."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time
import numpy as np
from scipy.integrate import quad
from reference import Gardner, analytical, step_q, step_u, gardner_fv, vg_column, mechanical_fv
from compare import QUANTITIES, validate


def orders(errors):
    return [float(np.log(a/b)/np.log(2)) for a, b in zip(errors[:-1], errors[1:])]


def water_level_exact(end, depth=.5, specific_yield=.2):
    excess = sum(quad(lambda t: step_q(depth, t, jump), 0, end-at,
                      epsabs=1e-13, epsrel=1e-11)[0]
                 for at, jump in ((0., 2e-5), (300., -2e-5)) if end > at)
    return .2+excess/specific_yield


def run(out, pc=False):
    start = time.perf_counter()
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    p = Gardner()
    meshes = [16, 32, 64, 128] if not pc else [64, 128, 256, 512]
    space = []
    for n in meshes:
        r = gardner_fv(n, .5 if not pc else .125, method='CN')
        half = gardner_fv(n, r['dt']/2, method='CN')
        u, q = analytical(r['x'], 600.)
        _, faces = analytical(np.arange(n+1)*p.length/n, 600.)
        level = .2+r['history'][-1, 5]/.2
        mechanical = mechanical_fv(16, level)
        exact_mechanical = mechanical_fv(16, water_level_exact(600.))
        hydraulic_stress = mechanical['effective_stress']-exact_mechanical['effective_stress']
        space.append(dict(n=n, dt_s=r['dt'], method='CN',
                          head_rms_m=float(np.sqrt(np.mean((r['psi']-np.log(u)/p.alpha)**2))),
                          face_flux_rms_m_s=float(np.sqrt(np.mean((r['q']-faces)**2))),
                          half_dt_head_difference_m=float(np.sqrt(np.mean((r['psi']-half['psi'])**2))),
                          maximum_mass_defect_m=r['max_mass_m'],
                          storage_change_m=float(p.delta*p.length/n*np.sum(r['u']-p.initial_u)),
                          cumulative_top_infiltration_m=float(r['history'][-1, 2]),
                          cumulative_bottom_outflow_m=float(r['history'][-1, 3]),
                          final_x_m=r['x'].tolist(), final_pressure_head_m=r['psi'].tolist(),
                          analytical_pressure_head_m=(np.log(u)/p.alpha).tolist(),
                          final_face_flux_down_m_s=r['q'].tolist(),
                          analytical_face_flux_down_m_s=faces.tolist(),
                          water_level_m=level, water_level_exact_m=water_level_exact(600.),
                          hydraulic_stress_max_kpa=float(np.max(np.abs(hydraulic_stress))),
                          mechanical_stress_max_kpa=mechanical['max_stress_error'],
                          mechanical_displacement_max_m=mechanical['max_displacement_error'],
                          equilibrium_residual_kpa=mechanical['equilibrium_residual']))
    # Fixed space and successive-time differences avoid a spatial error floor.
    times = [20., 10., 5., 2.5]
    temporal_runs = [gardner_fv(128, dt) for dt in times]
    differences = [float(np.sqrt(np.mean((a['psi']-b['psi'])**2)))
                   for a, b in zip(temporal_runs[:-1], temporal_runs[1:])]
    temporal = dict(n=128, method='BE', dt_s=times,
                    successive_head_rms_m=differences, observed_orders=orders(differences),
                    maximum_mass_defect_m=[r['max_mass_m'] for r in temporal_runs])
    vg = [vg_column(n, 100.) for n in ([12, 24, 48] if not pc else [48, 96, 192])]
    result = dict(parameters=asdict(p), analytical_domain='semi-infinite',
                  pulse_jumps_m_s=[[0., 2e-5], [300., -2e-5]],
                  space=space, head_space_orders=orders([r['head_rms_m'] for r in space]),
                  flux_space_orders=orders([r['face_flux_rms_m_s'] for r in space]),
                  time=temporal, deep_boundary_u_increment=step_u(p.length, 600., 2e-5),
                  vg_uniform_top_analogue=vg,
                  original_five_accepted=False, new_comparison_accepted=None,
                  tolerance='未定義', elapsed_s=time.perf_counter()-start)
    (out/'summary.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    packet = make_packet(16)
    validate(packet)
    (out/'reference-packet.json').write_text(json.dumps(packet, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    import csv
    with (out/'reference-packet.csv').open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['quantity', 'id', 'time_s', 'coordinate_m', 'value'], lineterminator='\n')
        writer.writeheader()
        writer.writerows(packet['rows'])
    (out/'reference-packet.metadata.json').write_text(json.dumps(packet['metadata'], ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    return result


def make_packet(n=16):
    """Analytical hydraulic fields + exact recharge receptor; no body output."""
    p = Gardner()
    x = (np.arange(n)+.5)*p.length/n
    faces = np.arange(n+1)*p.length/n
    z = -np.arange(5)/4
    locations = {
        'head': [{'id': f'c{i}', 'coordinate_m': float(a)} for i, a in enumerate(x)],
        'saturation': [{'id': f'c{i}', 'coordinate_m': float(a)} for i, a in enumerate(x)],
        'flux': [{'id': f'f{i}', 'coordinate_m': float(a)} for i, a in enumerate(faces)],
        'infiltration': [{'id': 'top', 'coordinate_m': 0.}],
        'effective_stress': [{'id': f'r{i}', 'coordinate_m': float(a)} for i, a in enumerate(z)]}
    parameters = asdict(p)
    parameters.update(pulse_jumps=[[0., 2e-5], [300., -2e-5]],
                      specific_yield=.2, receptor_initial_level_m=.2, recharge_depth_m=.5,
                      receptor_modulus_kpa=10000., receptor_total_stress_kpa=50.)
    metadata = dict(contract='external-rainfall-reference/1', case_id='semi-infinite-pulse-receptor',
                    parameters=parameters, source_kind='independent_reference', complete=True,
                    clock='physical_seconds', event_trace='right_limit_boundary_continuous_interior',
                    interpolation='none', quantities=QUANTITIES,
                    times_s=[0., 300., 600.], locations=locations,
                    hydraulic_model='Gardner_unsaturated_semi_infinite; depth_down; elevation=-depth',
                    mechanical_model='one_way_recharge_bucket; receptor_elevation_up; fixed_total_Terzaghi')
    rows = []
    for t in metadata['times_s']:
        u, _ = analytical(x, t)
        _, flux = analytical(faces, t)
        # At a flux discontinuity, rows use the right limit of imposed rainfall.
        flux[0] = p.ks*p.initial_u+sum(j for at, j in ((0., 2e-5), (300., -2e-5)) if at <= t)
        cumulative = p.ks*p.initial_u*t+2e-5*min(t, 300.)
        level = water_level_exact(t)
        values = {'head': np.log(u)/p.alpha-x,
                  'saturation': (p.theta_r+p.delta*u)/p.theta_s,
                  'flux': flux, 'infiltration': [cumulative],
                  'effective_stress': 50.-9.81*(level-z)}
        for quantity, array in values.items():
            for location, value in zip(locations[quantity], array):
                rows.append(dict(quantity=quantity, id=location['id'], time_s=t,
                                 coordinate_m=location['coordinate_m'], value=float(value)))
    return dict(metadata=metadata, rows=rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True)
    parser.add_argument('--pc-series', action='store_true')
    args = parser.parse_args()
    r = run(args.out, args.pc_series)
    print(json.dumps({k: r[k] for k in ('head_space_orders', 'flux_space_orders', 'elapsed_s')}, indent=2))
