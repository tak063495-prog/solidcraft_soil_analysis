"""Independent rate-form projected-volume construction with frozen history.

Continues an explicitly identified original-input projected initial calculation.
No Rust state, pressure, displacement or holding vector generates the solution.
The changed spatial discretization remains research-only, not P6 acceptance.
"""
import argparse
import hashlib
import json
from pathlib import Path
from mcc_projected_volume_diagnostic import ProjectedSystem, np, rate, ref, read, sha
import mcc_projected_volume_diagnostic as projected


def physical_case(value):
    out = {k: v for k, v in value.items() if k not in ('schema', 'settings', 'verification_reject_once_stages')}
    out['kinematics'] = value.get('settings', {}).get('kinematics', 'side_rollers')
    out['stages'] = [{k: v for k, v in stage.items() if k != 'initial_increment'} for stage in value['stages']]
    return out


def integrals(system, state, seed, u):
    theta = (seed[:, :, 7]-state[:, :, 7])[system.active]
    div = np.einsum('eqai,ei->eqa', system.b, u[system.dofs])[:, :, :3].sum(axis=2)
    local = np.einsum('eqa,eq,eq->ea', system.n, div+theta, system.weight)
    out = np.zeros(len(system.keys))
    np.add.at(out, system.pi.reshape(-1), local.reshape(-1))
    scale = np.zeros_like(out)
    local_scale = np.einsum('eqa,eq,eq->ea', system.n, abs(div)+abs(theta), system.weight)
    np.add.at(scale, system.pi.reshape(-1), local_scale.reshape(-1))
    return out, scale, theta


def solve(mesh_file, ground_file, case_file, identity_file, initial_file, output, steps,
          system_factory=ProjectedSystem, extra_sources=()):
    if steps < 2:
        raise ValueError('at least two rate-integration steps required')
    mesh, ground, case, identity, initial = map(read, (mesh_file, ground_file, case_file, identity_file, initial_file))
    if initial['schema'] != 'soil-fem-projected-volume-mcc-diagnostic/1' or initial['load_fraction'] != 1:
        raise ValueError('an original-input projected full initialization is required')
    backend = system_factory.__module__+'.'+system_factory.__name__
    if initial.get('reference_backend', ProjectedSystem.__module__+'.'+ProjectedSystem.__name__) != backend:
        raise ValueError('initial/reference solver backend differs')
    bindings = initial['input_and_source_sha256']
    for name, digest in bindings.items():
        assert sha(Path(name)) == digest, name
    by_name = {Path(name).name: Path(name) for name in bindings}
    for current, name in ((mesh_file, 'mesh.json'), (ground_file, 'ground.json'), (identity_file, 'mesh-identity.json')):
        assert sha(current) == sha(by_name[name]), name
    assert physical_case(case) == physical_case(read(by_name['case.json'])), 'initial physical case differs'
    assert sha(Path(projected.__file__)) == sha(by_name['mcc_projected_volume_diagnostic.py'])
    cfg = ground['initialization'].get('geostatic')
    if not cfg or cfg['initial_boundary_factor'] != 0 or ground['initialization']['water_level_m'] is not None:
        raise ValueError('diagnostic requires a dry zero-boundary-factor geostatic initial state')
    if any(case.get(k) for k in ('joint_system', 'anchor_system', 'beam_system', 'initially_inactive_bodies')):
        raise ValueError('unsupported members, contacts or initial inactive domains')
    source_files = [Path(__file__), Path(projected.__file__), Path(rate.__file__), Path(ref.__file__),
                    Path(projected.mcc_rate_midpoint.__file__), Path(projected.verify_high_order.__file__), Path(projected.high_order_io.__file__)]
    source_files += list(extra_sources)
    inputs = [mesh_file, ground_file, case_file, identity_file, initial_file]
    fingerprints = {str(p): sha(p) for p in inputs+source_files}
    active = np.ones(len(mesh['tetrahedra']), dtype=bool)
    state = np.asarray(initial['material_state']).copy()
    seed = rate.seed(mesh, ground, identity)
    u = np.asarray(initial['displacement_m']).reshape(-1).copy()
    origin = u.copy()
    system = system_factory(mesh, ground, active, case)
    offset = np.zeros(len(system.keys))
    pressure = ground['initialization']['surface_surcharge_kpa']
    factor = 0.
    frames, release_audits = [], []

    def load():
        out = rate.physical(mesh, ground, active, pressure, case.get('pressure_on_original_surface_only', False))[0]
        out += factor*ref.conditions_load(mesh, case.get('boundary') or {}, active, skip_removed_surfaces=True)
        return out.reshape(-1)

    def frame(name, previous, diagnostics):
        local = state[active].reshape(-1, 14)
        force = system.force(local)
        physical = load()
        scale = max(np.linalg.norm(physical), 1.)
        residual = float(np.linalg.norm(system.reduce(physical-force))/scale)
        defect, volume_scale, theta = integrals(system, state, seed, u)
        compatibility = float(np.linalg.norm(defect-offset)/max(np.linalg.norm(volume_scale), float(system.weight.sum())*1e-12))
        assert residual < 1e-8 and compatibility < 1e-8
        p, _ = system.pressure(local)
        nodal_theta = system.project_volume(theta)
        theta_error = float(abs(np.einsum('eqa,ea->eq', system.n, nodal_theta[system.pi])-theta).max())
        assert theta_error < 1e-10
        used = np.unique(np.asarray(mesh['tetrahedra10'])[active])
        monitors = {}
        for label, key in case.get('monitor_sets', {}).items():
            nodes = np.intersect1d(mesh['node_sets'][key], used)
            if len(nodes):
                z = u.reshape(-1, 3)[nodes, 2]
                monitors[label] = dict(nodes=len(nodes), mean_from_seed_m=float(z.mean()),
                    mean_from_construction_origin_m=float((u-origin).reshape(-1, 3)[nodes, 2].mean()),
                    mean_stage_increment_m=float((u-previous).reshape(-1, 3)[nodes, 2].mean()))
        frames.append(dict(name=name, displacement_m=u.reshape(-1, 3).tolist(), active_elements=active.tolist(),
            material_state=state.tolist(), pressure_dofs=[dict(node_index=k[0], material_id=k[1], compression_positive_kpa=float(v)) for k, v in zip(system.keys, p)],
            volume_dofs=[dict(node_index=k[0], material_id=k[1], compression_positive=float(v)) for k, v in zip(system.keys, nodal_theta)],
            compatibility_offsets=[dict(node_index=k[0], material_id=k[1], integral_m3=float(v)) for k, v in zip(system.keys, offset)],
            relative_work_conjugate_force_residual=residual, relative_frozen_volume_constraint_residual=compatibility,
            maximum_nodal_gp_volume_difference=theta_error,
            plastic_points=int(np.count_nonzero(np.linalg.norm(local[:, 8:], axis=1)>1e-10)),
            monitors=monitors, diagnostics=diagnostics))
        if hasattr(system, 'solver_diagnostics'):
            frames[-1]['sparse_solver_diagnostics'] = system.solver_diagnostics()

    frame('initial_equilibrated', u.copy(), initial['diagnostics'])
    for stage in case['stages']:
        previous = u.copy()
        old_state = state.copy()
        if stage['kind'] == 'pressure':
            pressure = stage['target_pressure_kpa']
        elif stage['kind'] == 'boundary':
            factor = stage['target_factor']
        elif stage['kind'] == 'excavate':
            remove = []
            for body in stage['bodies']:
                if body in mesh.get('excavation_sets', {}):
                    remove.extend(mesh['excavation_sets'][body])
                else:
                    regions = {r['id'] for r in mesh['regions'] if r['body'] == body}
                    remove.extend(e for e, r in enumerate(mesh['region_ids']) if r in regions)
            remove = np.unique(remove)
            if not len(remove) or not active[remove].all():
                raise ValueError('unknown, empty or overlapping excavation')
            active[remove] = False
            pressure = stage['remaining_surface_pressure_kpa']
            system = system_factory(mesh, ground, active, case)
            offset = integrals(system, state, seed, u)[0]
        else:
            raise ValueError('only pressure/boundary/excavation diagnostic stages are supported')
        before = system.force(state[active].reshape(-1, 14))
        after = load()
        if stage['kind'] == 'excavate':
            release_audits.append(dict(stage=stage['name'], removed_elements=remove.tolist(),
                holding_nodal_load_kn=(before-after).reshape(-1, 3).tolist(),
                compatibility_offsets=[dict(node_index=k[0], material_id=k[1], integral_m3=float(v)) for k, v in zip(system.keys, offset)],
                material_history_unchanged_at_domain_restriction=bool(np.array_equal(state, old_state))))
        local, u, diagnostics = system.advance(state[active].reshape(-1, 14), u, after-before, steps,
            progress=lambda t, n, e: print(f'{stage["name"]}: {t:.6f}, steps {n}, events {e}', flush=True))
        state[active] = local.reshape(state[active].shape)
        assert np.array_equal(state[~active], old_state[~active]), 'removed material history changed'
        frame(stage['name'], previous, diagnostics)
        print(json.dumps(dict(stage=stage['name'], plastic_points=frames[-1]['plastic_points'],
            force_residual=frames[-1]['relative_work_conjugate_force_residual'],
            volume_residual=frames[-1]['relative_frozen_volume_constraint_residual'], monitors=frames[-1]['monitors'])), flush=True)
    assert fingerprints == {str(p): sha(p) for p in inputs+source_files}, 'source/input changed'
    report = dict(schema='soil-fem-projected-mcc-cycle-diagnostic/1', diagnostic_only=True,
        production_backend_accepted=False, goal_complete=False, p6_complete=False, reference_steps=steps,
        state_origin='Independent projected initial reference; no Rust solved field or history',
        initial_reference_sha256=sha(initial_file), reference_displacement_m=origin.reshape(-1, 3).tolist(),
        frames=frames, release_audits=release_audits, input_and_source_sha256=fingerprints,
        limitations='No spatial-convergence acceptance; explicit rate integrator differs from Rust fixed-partition updates; no members/water/contact or production export')
    report['reference_backend'] = backend
    with output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, allow_nan=False)
    print('Independent projected-volume construction diagnostic completed', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ('mesh', 'ground', 'case', 'identity', 'initial', 'output'):
        p.add_argument(name, type=Path)
    p.add_argument('--steps', type=int, default=100)
    a = p.parse_args()
    solve(*(getattr(a, k).resolve() for k in ('mesh', 'ground', 'case', 'identity', 'initial', 'output')), a.steps)
