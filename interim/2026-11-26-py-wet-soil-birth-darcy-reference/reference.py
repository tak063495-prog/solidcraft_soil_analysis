"""Independent 1D free-boundary Darcy FV and vertical birth-force reference.

Only numpy/scipy/stdlib; equations and intentionally restricted physics in README.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

GW, GD, GS, TOP, JOINT, SURCHARGE = 9.81, 17., 20., 2., .7, 3.
K0, RECHARGE, BOTTOM_HEAD = 1e-5, 1e-6, 1.


def conductivity(z, alpha):
    return K0*np.exp(alpha*np.asarray(z))


def exact_head(z, alpha):
    # Integral of dH/dz = recharge/K, not used by the FV assembly.
    return BOTTOM_HEAD + RECHARGE/K0*(-np.expm1(-alpha*np.asarray(z)))/alpha


def exact_surface(alpha):
    return brentq(lambda s: float(exact_head(s, alpha))-s, .71, 1.9, xtol=1e-14)


def fixed_domain(s, n, alpha):
    """Vertex-centred FV: midpoint face K, bottom head, top incoming flux.

    Unknowns are H[1:n+1]; rows are local flux balances, including top flux.
    """
    z = np.linspace(0., s, n+1)
    dz = np.diff(z)
    t = conductivity((z[:-1]+z[1:])/2, alpha)/dz
    diagonal = np.r_[t[:-1]+t[1:], t[-1]]
    matrix = diags((-t[1:], diagonal, -t[1:]), (-1, 0, 1), format='csc')
    rhs = np.zeros(n)
    rhs[0] += t[0]*BOTTOM_HEAD
    rhs[-1] += RECHARGE
    h = np.r_[BOTTOM_HEAD, spsolve(matrix, rhs)]
    flux = -t*np.diff(h)  # z positive; m/s, A=1 m2
    return z, h, flux


def solve(n, alpha=.5):
    if type(n) is not int or not 4 <= n <= 16384 or alpha not in (.5, -.25):
        raise ValueError('supported n=4..16384 and alpha=.5 or -.25 only')
    s = brentq(lambda trial: fixed_domain(trial, n, alpha)[1][-1]-trial,
               .71, 1.9, xtol=1e-14)
    z, h, q = fixed_domain(s, n, alpha)
    if np.min(h-z) < -1e-11:
        raise ValueError('negative saturated pressure; never clip')
    return dict(surface=s, z=z, head=h, flux=q)


def birth(solution):
    """A=1 m2 P1 virtual-work forces on old/new, separate interface IDs.

    Explicit linear reconstruction on each computed hydraulic interval, split at
    free surface and birth interface. No production-output interpolation.
    """
    s, z, head = solution['surface'], solution['z'], solution['head']
    profile_p = GW*(head-z)
    samples, forces = [], []
    interface = {}
    for body, low, high in [('old', 0., JOINT), ('new', JOINT, TOP)]:
        zz = np.unique(np.r_[low, z[(z > low) & (z < high)], high, s if low < s < high else []])
        p = np.interp(zz, z, profile_p, right=0.)
        total = SURCHARGE+GD*(TOP-np.maximum(zz, s))+GS*np.maximum(s-zz, 0.)
        effective = total-p
        count = len(zz)
        fields = {key: np.zeros(count) for key in
                  ('effective_internal', 'pore_internal', 'weight', 'pore_boundary',
                   'effective_boundary', 'effective_body', 'seepage', 'buoyancy')}
        for e, dz in enumerate(np.diff(zz)):
            sl = slice(e, e+2)
            wet = (zz[e]+zz[e+1])/2 < s
            gamma = GS if wet else GD
            dp = (p[e+1]-p[e])/dz
            fields['effective_internal'][sl] += np.array([1., -1.])*np.mean(effective[sl])
            fields['pore_internal'][sl] += np.array([1., -1.])*np.mean(p[sl])
            fields['weight'][sl] += -gamma*dz/2
            fields['effective_body'][sl] += (-gamma-dp)*dz/2
            fields['buoyancy'][sl] += (GW if wet else 0.)*dz/2
            # discrete total head gradient; includes FV pressure reconstruction
            fields['seepage'][sl] += (-dp-GW if wet else 0.)*dz/2
        for j, normal in ((0, -1.), (-1, 1.)):
            fields['pore_boundary'][j] = -normal*p[j]
            fields['effective_boundary'][j] = -normal*effective[j]
        for i, position in enumerate(zz):
            sample_id = f'{body}:{i}'
            samples.append(dict(id=sample_id, body=body, z_m=float(position),
                                pressure_kpa=float(p[i]), effective_kpa=float(effective[i]),
                                total_kpa=float(total[i])))
            forces.append(dict(id=sample_id, **{key+'_kn': float(value[i]) for key, value in fields.items()}))
        interface[body] = float(p[-1] if body == 'old' else p[0])
    residual, pore_identity, effective_identity = [], [], []
    for f in forces:
        residual.append(f['effective_internal_kn']+f['pore_internal_kn']-f['weight_kn']
                        -f['pore_boundary_kn']-f['effective_boundary_kn'])
        pore_identity.append(f['pore_internal_kn']-f['pore_boundary_kn']
                             +f['effective_body_kn']-f['weight_kn'])
        effective_identity.append(f['effective_body_kn']-f['weight_kn']-f['buoyancy_kn']-f['seepage_kn'])
    # Fixed post-birth trace topology audit, not a pre-birth flow solve.
    topology = dict(before_old_exposed_kn=-interface['old'],
                    after_old_interface_kn=-interface['old'],
                    after_new_interface_kn=interface['new'],
                    stale_old_exterior_after_kn=0.,
                    paired_interface_sum_kn=interface['new']-interface['old'])
    return samples, forces, dict(nodal_residual_max_kn=max(map(abs, residual)),
                                pore_identity_max_kn=max(map(abs, pore_identity)),
                                buoyancy_seepage_identity_max_kn=max(map(abs, effective_identity)),
                                topology=topology)


def export(n=32, alpha=.5):
    sol = solve(n, alpha)
    samples, forces, force_metrics = birth(sol)
    exact_s = exact_surface(alpha)
    joint = next(row for row in samples if row['id'] == 'new:0')
    exact_joint_p = float(GW*(exact_head(JOINT, alpha)-JOINT))
    exact_joint_eff = SURCHARGE+GD*(TOP-exact_s)+GS*(exact_s-JOINT)-exact_joint_p
    head_error = sol['head']-exact_head(sol['z'], alpha)
    budget = [dict(id='bottom', outward_m3_s=float(-sol['flux'][0])),
              dict(id='recharge_at_free_surface', outward_m3_s=float(sol['flux'][-1]))]
    return dict(metadata=dict(contract='darcy-birth-reference/1', producer='INDEPENDENT_REFERENCE_ONLY',
                              complete=True, status='completed', n=n, alpha=alpha,
                              units=dict(length='m', pressure='kPa', force='kN', flow='m3/s'),
                              axis='z_up', stress_sign='compression_positive',
                              flow_sign='outward_positive', force_sign='physical_z_positive',
                              stress_kind='effective_and_total_separate', head_kind='total',
                              interpolation=False, pressure_clipped=False, pressure_filled=False,
                              area_m2=1., stage='birth_initial_equilibrium', generation=1),
                surface_m=sol['surface'],
                hydraulic=[dict(id=f'h:{i}', z_m=float(z), head_m=float(h), pressure_kpa=float(GW*(h-z)))
                           for i, (z, h) in enumerate(zip(sol['z'], sol['head']))],
                samples=samples, forces=forces, budget=budget,
                metrics=dict(surface_error_m=abs(sol['surface']-exact_s),
                             head_rms_m=float(np.sqrt(np.mean(head_error**2))),
                             joint_pressure_error_kpa=abs(joint['pressure_kpa']-exact_joint_p),
                             joint_effective_error_kpa=abs(joint['effective_kpa']-exact_joint_eff),
                             flux_max_error_m_s=float(max(abs(sol['flux']+RECHARGE))),
                             local_flux_spread_m3_s=float(np.ptp(sol['flux'])),
                             water_defect_m3_s=sum(b['outward_m3_s'] for b in budget),
                             normalized_water_defect=abs(sum(b['outward_m3_s'] for b in budget))/max(sum(abs(b['outward_m3_s']) for b in budget),1e-30),
                             free_surface_pressure_kpa=float(GW*(sol['head'][-1]-sol['surface'])),
                             **force_metrics),
                accepted=None, a03_complete=False)


def policies(pressure_kpa, head_a, head_b, prescribed_head):
    """Trace audit only; SOIL is a vertical impermeable side wall (normal x)."""
    if head_a != head_b or prescribed_head != head_a:
        raise ValueError('conflicting traces: no averaging')
    return dict(continuous=dict(head_a=head_a, head_b=head_b, normal_traction_kn=[-pressure_kpa, pressure_kpa]),
                prescribed=dict(fixed_head=prescribed_head, normal_traction_kn=[-pressure_kpa, pressure_kpa],
                                reservoir_exchange_definition='sum outward fluxes of both adjoining subdomains'),
                soil=dict(soil_side='a', structure_pore_pressure_kpa=0., normal_flux_m3_s=0., normal_axis='x',
                          normal_traction_kn=[-pressure_kpa, pressure_kpa],
                          throughflow_column_applicable=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--pc-series', action='store_true')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    rows = []
    for alpha in (.5, -.25):
        previous = None
        for n in ((32, 64, 128, 256, 512, 1024) if args.pc_series else (8, 16, 32, 64)):
            data = export(n, alpha)
            (args.output/f'alpha{alpha}-n{n}.json').write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
            row = dict(alpha=alpha, n=n, surface_m=data['surface_m'], **data['metrics'])
            for key in ('surface_error_m', 'head_rms_m', 'joint_pressure_error_kpa', 'joint_effective_error_kpa'):
                row[key+'_order'] = None if previous is None else float(np.log2(previous[key]/row[key]))
            previous = row
            rows.append(row)
    (args.output/'convergence.json').write_text(json.dumps(dict(rows=rows, accepted=None, a03_complete=False), indent=2)+'\n')
    case = export(32)
    sample = next(x for x in case['samples'] if x['id'] == 'new:0')
    pressure = sample['pressure_kpa']
    head = JOINT+pressure/GW
    audit = dict(roles=policies(pressure, head, head, head), accepted=None,
                 scope='trace roles only; no new constitutive or policy-specific birth solve')
    (args.output/'policy-audit.json').write_text(json.dumps(audit, indent=2)+'\n')


if __name__ == '__main__':
    main()
