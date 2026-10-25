"""Physical/operator checks for the research-only volume projection.

No production result or production constitutive update is used.
"""
import argparse
import hashlib
import json
from pathlib import Path
from mcc_projected_volume_diagnostic import ProjectedSystem, np, rate
import mcc_projected_volume_diagnostic as projected
from prepare_mixed_mcc_reference import column
import prepare_mixed_mcc_reference as preparation


def verify(output):
    mesh = column(2)
    xyz = np.asarray(mesh['nodes'])
    for name, axis in (('xmin', 0), ('ymin', 1)):
        mesh['node_sets'][name] = np.flatnonzero(xyz[:, axis] == 0.).tolist()
    case = dict(settings=dict(kinematics='explicit'), boundary=dict(constraints=[
        dict(node_set=name, direction=np.eye(3)[axis].tolist(), value_m=0.)
        for name, axis in (('xmin', 0), ('ymin', 1), ('bottom', 2))]))
    material = dict(critical_state_m=1.2, **{'lambda': .2, 'kappa': .04}, poisson=.3)
    ground = dict(materials=[dict(id=i, mcc=material.copy()) for i in (1, 2)])
    active = np.ones(len(mesh['tetrahedra']), dtype=bool)
    system = ProjectedSystem(mesh, ground, active, case)
    # Distinct constant pressure in two materials at the same geometric nodes.
    state = np.zeros((system.weight.size, 14))
    means = np.repeat([10. if m == 1 else 30. for m in mesh['material_ids']], system.weight.shape[1])
    state[:, :3] = means[:, None]
    pressure, _ = system.pressure(state)
    expected = np.array([10. if key[1] == 1 else 30. for key in system.keys])
    pressure_error = float(abs(pressure-expected).max())
    shared = set(n for n, m in system.keys if m == 1) & set(n for n, m in system.keys if m == 2)
    assert shared and pressure_error < 1e-10
    # Isotropic NC compression: v=N-lambda ln(p), pc=p, theta=-ln(v/v0).
    initial = np.zeros_like(state)
    initial[:, :3] = 100.
    initial[:, 6] = np.log(100.)
    v0 = 2.7
    initial[:, 7] = np.log(v0)
    target = initial.copy();target[:, :3] = 130.
    change = system.force(target)-system.force(initial)
    expected_v = v0-.2*np.log(1.3)
    theta = -np.log(expected_v/v0)
    expected_u = -theta/3*xyz
    results = []
    for count in (1000, 2000):
        current = ProjectedSystem(mesh, ground, active, case)
        final, u, diagnostics = current.advance(initial, np.zeros(xyz.size), change, count)
        v_error = float(abs(np.exp(final[:, 7])-expected_v).max())
        pc_error = float(abs(np.exp(final[:, 6])-130.).max())
        u_error = float(abs(u.reshape(-1, 3)-expected_u).max())
        ep_vol = final[:, 8:11].sum(axis=1)
        ep_error = float(abs(ep_vol-(1-.04/.2)*theta).max())
        force_error = float(np.linalg.norm(current.reduce(current.force(final)-current.force(target)))/np.linalg.norm(current.reduce(current.force(target))))
        results.append(dict(steps=count, specific_volume_maximum_error=v_error,
            pc_maximum_error_kpa=pc_error, displacement_maximum_error_m=u_error,
            plastic_volume_maximum_error=ep_error, relative_force_residual=force_error,
            operator_checks=current.checks, diagnostics=diagnostics))
        print(json.dumps(results[-1]), flush=True)
        assert force_error < 1e-10 and pc_error < 1e-6, results[-1]
    assert results[1]['specific_volume_maximum_error'] < results[0]['specific_volume_maximum_error']
    assert results[1]['displacement_maximum_error_m'] < results[0]['displacement_maximum_error_m'] < 1e-7
    # Rebuilding projection after DEATH changes its weak compatibility test
    # space. Preserve old material history; do not silently replace cumulative
    # theta by the newly projected total displacement in a future controller.
    virtual = np.random.default_rng(1251025).normal(size=xyz.size)*1e-4
    old_theta_gp = system.strain(virtual)[:, :, :3].sum(axis=2)
    removed = np.asarray(mesh['material_ids']) == 2
    remaining = ProjectedSystem(mesh, ground, ~removed, case)
    new_theta_gp = remaining.strain(virtual)[:, :, :3].sum(axis=2)
    remap_difference = float(abs(old_theta_gp[~removed]-new_theta_gp).max())
    # Removing a whole material block keeps the split projection unchanged.
    assert remap_difference < 1e-12
    removed = np.zeros(len(active), dtype=bool);removed[0] = True
    local_cut = ProjectedSystem(mesh, ground, ~removed, case)
    local_new = local_cut.strain(virtual)[:, :, :3].sum(axis=2)
    local_difference = float(abs(old_theta_gp[~removed]-local_new).max())
    assert local_difference > 1e-8
    sources = [Path(__file__), Path(projected.__file__), Path(preparation.__file__), Path(rate.__file__),
               Path(projected.ref.__file__), Path(projected.mcc_rate_midpoint.__file__),
               Path(projected.verify_high_order.__file__), Path(projected.high_order_io.__file__)]
    report = dict(schema='soil-fem-projected-volume-research-checks/1', passed=True,
        diagnostic_only=True, production_backend_accepted=False, goal_complete=False, p6_complete=False,
        material_split=dict(shared_geometric_nodes=len(shared), maximum_pressure_error_kpa=pressure_error),
        isotropic_nc_theory=dict(p0_kpa=100., p1_kpa=130., lambda_value=.2, kappa_value=.04,
            expected_specific_volume=expected_v, expected_compression_positive_volume_strain=theta, refinement=results),
        domain_change=dict(whole_material_block_projection_difference=remap_difference,
            local_death_projection_difference=local_difference,
            consequence='Local DEATH requires a frozen compatibility offset or increment-based constraint; replacing old total theta by the new total-displacement projection would alter material history'),
        sha256={str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    with output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps({k: report[k] for k in ('passed', 'material_split', 'isotropic_nc_theory', 'domain_change')}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser();p.add_argument('output', type=Path)
    verify(p.parse_args().output.resolve())
