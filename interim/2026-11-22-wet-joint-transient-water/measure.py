"""Full accepted-field measurements; no new numerical acceptance threshold."""
import argparse
from collections import Counter
import importlib.util
import itertools
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('c04_study_measure', HERE/'study.py')
s = importlib.util.module_from_spec(spec); spec.loader.exec_module(s)
read, write = s.read, s.write
CATEGORIES = ('material_return', 'contact_status', 'open')


def fields(folder):
    with (folder/'fields.jsonl').open() as f:
        for line in f:
            v = json.loads(line)
            yield v['row'], {k: np.asarray(x) for k, x in v['fields'].items()}


def metric(a, b, category=False):
    a = np.asarray(a); b = np.asarray(b)
    assert a.shape == b.shape and a.size and np.isfinite(a).all() and np.isfinite(b).all()
    if category:
        return dict(components=a.size, mismatches=int(np.count_nonzero(a != b)))
    delta = a.astype(float)-b.astype(float)
    axes = tuple(range(1, a.ndim))
    return dict(components=a.size, difference_squared=float(np.sum(delta**2)), reference_squared=float(np.sum(b.astype(float)**2)),
                maximum_absolute=float(np.max(np.abs(delta))),
                entity_difference_squared=np.sum(delta**2, axis=axes).tolist(),
                entity_reference_squared=np.sum(b.astype(float)**2, axis=axes).tolist(),
                entity_maximum_absolute=np.max(np.abs(delta), axis=axes).tolist())


def add(target, key, a, b):
    v = metric(a, b, key in CATEGORIES)
    if key not in target:
        target[key] = v; return
    old = target[key]
    for k, x in v.items():
        if k in ('maximum_absolute', 'entity_maximum_absolute'):
            old[k] = np.maximum(old[k], x).tolist()
        elif isinstance(x, list): old[k] = (np.asarray(old[k])+x).tolist()
        else: old[k] += x


def finish(store):
    for v in store.values():
        if 'mismatches' in v:
            v['mismatch_fraction'] = v['mismatches']/v['components']
        else:
            v['absolute_l2'] = v['difference_squared']**.5
            v['relative_l2'] = v['absolute_l2']/max(v['reference_squared']**.5, 1e-12)
    return store


def pair(a, b):
    """Compare every state at common prescribed progress, without interpolation."""
    left = {(r['stage'], r['progress']): (r, f) for r, f in fields(a)}
    right = {(r['stage'], r['progress']): (r, f) for r, f in fields(b)}
    assert left and set(left) <= set(right), 'missing common accepted point'
    pooled = {}; phases = {}; coordinates = []
    for key, (ra, fa) in left.items():
        rb, fb = right[key]
        assert ra['factor'] == rb['factor'] and ra['time_s'] == rb['time_s']
        assert fa.keys() == fb.keys()
        for k in fa:
            add(pooled, k, fa[k], fb[k]); add(phases.setdefault(key[0], {}), k, fa[k], fb[k])
        coordinates.append(dict(stage=key[0], progress=key[1]))
    return dict(common_accepted_points=len(left), finer_unmatched_points=len(right)-len(left), interpolation=False,
                fields=finish(pooled), stages={k: finish(v) for k, v in phases.items()},
                common_coordinates=coordinates, numerical_acceptance=False)


def audit(out, folder):
    case = read(folder/'case.json'); model = s.mechanics.Model(read(out/'mesh.json'), read(out/'ground.json'), case)
    phases = s.schedule(case)
    expected = [(p['name'], i/p['n']) for p in phases for i in range(1, p['n']+1)]
    previous = None; count = 0; maxima = {}; rows = []; initial_water = None; minp = np.inf; maxp = -np.inf
    def defect(name, value): maxima[name] = max(maxima.get(name, 0.), float(np.max(np.abs(value), initial=0.)))
    for (row, f), coord in itertools.zip_longest(fields(folder), expected):
        assert (row['stage'], row['progress']) == coord
        assert all(np.isfinite(v).all() for v in f.values())
        model.set_level(row['level_m'])
        native = (folder/'raw').exists()
        gate = case['settings']['relative_force_tolerance'] if native else 1e-10
        force = float(np.linalg.norm(model.t.T@f['reaction'].ravel())/max(np.linalg.norm(f['applied_load']), 1.))
        assert force <= gate, ('original force gate', folder.name, row, force, gate)
        defect('maximum_reassembled_force_residual', force)
        defect('head_a_minus_b_m', f['joint_head_a']-f['joint_head_b'])
        defect('pressure_mapping_kpa', f['physical_pressure']-model.interface['pressure'])
        defect('interface_water_force_kn', f['water_load']-model.interface['load'])
        defect('bulk_water_force_kn', f['bulk_water_force']-model.pore_force.reshape(-1, 3))
        defect('exterior_water_force_kn', f['exterior_water_load']-model.exterior_water)
        defect('total_effective_stress_pore_relation_kpa', f['total_stress'][:, :3]-f['stress'][:, :3]-f['pore'][:, None])
        defect('gap_from_jump_m', f['gap']-np.einsum('qa,qa->q', f['jump'], model.system.normal))
        defect('contact_energy_work_identity_kpa_m', f['endpoint_work']-f['elastic_energy']-f['plastic_dissipation']-f['numerical_dissipation'])
        closed = f['gap'] <= 0; normal = f['effective_contact_pressure']
        defect('closed_normal_spring_kpa', normal[closed]+model.system.kn[closed]*f['gap'][closed])
        if (~closed).any(): defect('open_effective_traction_kpa', f['traction'][~closed])
        tangential = f['traction']+normal[:, None]*model.system.normal
        defect('friction_circle_excess_kpa', np.maximum(np.linalg.norm(tangential, axis=1)-model.system.mu*np.maximum(normal, 0), 0))
        assert not np.any(f['plastic_strain']) and not np.any(f['equivalent_plastic_strain']) and not np.any(f['plastic_work']), 'elastic bulk scope exceeded'
        if previous is not None:
            defect('accumulated_slip_decrease_m', np.minimum(f['accumulated']-previous['accumulated'], 0))
        if initial_water is None: initial_water = (f['physical_pressure'].copy(), f['water_load'].copy())
        defect('pressure_change_from_initial_kpa', f['physical_pressure']-initial_water[0])
        defect('interface_force_change_from_initial_kn', f['water_load']-initial_water[1])
        minp = min(minp, f['physical_pressure'].min()); maxp = max(maxp, f['physical_pressure'].max())
        rows.append(dict(stage=row['stage'], progress=row['progress'], reporting_time_s=row['time_s'], level_m=row['level_m'],
                         open_points=int(f['open'].sum()), maximum_slip_m=float(f['accumulated'].max()), relative_force_residual=force))
        previous = f; count += 1
    assert count == len(expected)
    return dict(accepted_states=count, nodes=len(model.xyz), material_gp=len(model.seeds), contact_gp=len(model.system.points), all_accepted_points=rows,
                defects=maxima, minimum_interface_pressure_kpa=float(minp), maximum_interface_pressure_kpa=float(maxp),
                hydrostatic_mass_balance=None, hydrostatic_mass_balance_reason='Prescribed equilibrated levels have no flow/storage PDE, physical flux or water-mass state. Zero balance is not inferred.',
                primary_nodal_head_unknown=False, reported_node_values='samples of committed prescribed level; no nodal water solve',
                bulk_water_force_provenance=('native diagnostic internal minus independently assembled saved constitutive/contact force' if native else 'input-only cut-volume integral'),
                exterior_force_provenance=('native applied load minus initial holding, dry/body/boundary and native interface forces' if native else 'input-only exterior hydrostatic integral'))


def local_contact(out, folder):
    """Saved-old-state replay is a local diagnostic, never a fresh 3D solution."""
    model = s.mechanics.Model(read(out/'mesh.json'), read(out/'ground.json'), read(folder/'case.json'))
    system = model.system; previous = None; maxima = {}; partitions = []; increments = 0
    names = [('traction', 'traction_kpa'), ('plastic_slip', 'plastic_slip_m'), ('accumulated', 'accumulated_plastic_slip_m'),
             ('plastic_dissipation', 'plastic_dissipation_kpa_m'), ('numerical_dissipation', 'numerical_dissipation_kpa_m'), ('endpoint_work', 'endpoint_work_kpa_m')]
    with (folder/'raw/accepted.jsonl').open() as sf, (folder/'raw/constitutive.jsonl').open() as df:
        for line, diagnostic in itertools.zip_longest(sf, df):
            assert line is not None and diagnostic is not None
            state = json.loads(line); d = json.loads(diagnostic); points = state['joints']['points']
            get = lambda ps, key: np.array([p[key] for p in ps])
            target = get(points, 'jump_m')
            start = np.zeros_like(target) if previous is None else get(previous, 'jump_m')
            oldq = np.zeros_like(target) if previous is None else get(previous, 'traction_kpa')
            counts = np.asarray(d['partitions'][len(model.seeds):]); assert len(counts) == len(target)
            partitions.extend(counts.tolist())
            for count in np.unique(counts):
                mask = counts == count
                sub = SimpleNamespace(**{k: getattr(system, k)[mask] for k in ('normal', 'kn', 'kt', 'mu')})
                q, dp, da = s.mechanics.contact(sub, start[mask], oldq[mask], target[mask], int(count))
                pd, nd, work = s.mechanics.contact_work(sub, start[mask], oldq[mask], target[mask], int(count))
                for (label, raw), value in zip(names, (q[:, 0], dp[:, 0], da[:, 0], pd, nd, work)):
                    actual = get(points, raw)[mask]
                    if label != 'traction' and previous is not None: actual -= get(previous, raw)[mask]
                    maxima[label] = max(maxima.get(label, 0.), float(np.max(np.abs(actual-value))))
            system.audit(state)
            previous = points; increments += 1
    return dict(diagnostic_only_not_fresh_3d_solve=True, increments=increments, maximum_local_difference=maxima,
                native_contact_partition_histogram=dict(sorted(Counter(partitions).items())))


def nodal_projection(mesh, face_keys, values):
    """Explicit observation operator, not an RT0 nodal unknown or JOINT map."""
    total = np.zeros(len(mesh['nodes'])); count = np.zeros(len(total))
    for key, value in zip(face_keys, values):
        total[list(key)] += value; count[list(key)] += 1
    assert np.all(count > 0)
    return total/count


def hydraulic_fields(out, folder, kind):
    mesh = read(out/'hydraulic-mesh.json'); native = read(folder/'native.json'); independent = read(folder/'independent.json')
    keys = [tuple(sorted(f['corners'])) for f in native['faces']]
    lookup = {tuple(k): i for i, k in enumerate(independent['face_keys'])}
    mapping = [lookup[k] for k in keys]
    if kind == 'native':
        for frame in native['frames'][1:]:
            v = frame['state']; case = native['case']; h = np.array(v['relative_face_head_m'])+case.get('head_datum_m', 0.)
            yield v['time_s'], dict(cell_head_m=np.array(v['relative_cell_head_m'])+case.get('head_datum_m', 0.), face_head_m=h,
                                    outflow_m3_s=np.array(v['cell_outflow_m3_s']), storage_offset_m3=np.array(v['cell_storage_offset_m3']),
                                    nodal_head_projection_m=nodal_projection(mesh, keys, h))
    else:
        for v in independent['rows']:
            h = np.asarray(v['face_head_m'])[mapping]
            yield v['time_s'], dict(cell_head_m=np.array(v['cell_head_m']), face_head_m=h, outflow_m3_s=np.array(v['outflow_m3_s']),
                                    storage_offset_m3=np.array(v['storage_offset_m3']), nodal_head_projection_m=nodal_projection(mesh, keys, h))


def hydraulic_pair(out, a, b, akind, bkind):
    left = dict(hydraulic_fields(out, a, akind)); right = dict(hydraulic_fields(out, b, bkind))
    assert set(left) <= set(right)
    store = {}
    for t, fa in left.items():
        for k in fa: add(store, k, fa[k], right[t][k])
    return dict(common_physical_times=len(left), interpolation=False, fields=finish(store), numerical_acceptance=False)


def hydraulic_audit(folder):
    a = read(folder/'native.json'); b = read(folder/'independent.json')
    initial = sum(a['frames'][0]['state']['cell_storage_offset_m3'])
    flux_integral = 0.; previous = 0.; defects = []; balances = []
    for frame in a['frames'][1:]:
        t = frame['state']['time_s']; dt = t-previous; bal = frame['balance']
        assert bal['maximum_scaled_residual'] <= 1., 'original Richards mass/flux gate'
        flux_integral += (bal['boundary_outflow_m3_s']-bal['source_m3_s'])*dt
        defect = sum(frame['state']['cell_storage_offset_m3'])-initial+flux_integral
        defects.append(abs(defect)); previous = t
        balances.append(dict(time_s=t, dt_s=dt, integrated_mass_defect_m3=defect, **bal))
    assert all(v['maximum_scaled_residual'] < .01 for v in b['rows']), 'independent stricter residual unchanged'
    return dict(native_maximum_global_mass_defect_m3=max(defects), independent_maximum_global_mass_defect_m3=max(abs(v['global_mass_defect_m3']) for v in b['rows']),
                native_maximum_scaled_residual=max(x['maximum_scaled_residual'] for x in balances),
                independent_maximum_scaled_residual=max(v['maximum_scaled_residual'] for v in b['rows']),
                native_maximum_cell_mass_residual_m3_s=max(x['maximum_absolute_cell_mass_residual_m3_s'] for x in balances),
                independent_maximum_cell_mass_residual_m3_s=max(v['maximum_absolute_cell_residual_m3_s'] for v in b['rows']),
                all_native_balance_steps=balances, native_accepted_steps=len(a['attempts']), independent_accepted_steps=len(b['rows']),
                primary_unknowns='RT0 cell P0 head and hybrid face trace. No nodal unknown. Explicit incident-face mean projection is diagnostic only.',
                hydraulic_to_joint_mapping=None, hydraulic_to_joint_mapping_reason='Existing ONE_WAY rejects JOINT; no native interface-transfer state exists.',
                storage_law='van Genuchten/Mualem with saturated specific storage; 12 Tet4, order-3 hydraulic integration (324 points).',
                unsaturated_point_states=sum(p['pressure_head_m'] < 0 for f in a['frames'][1:] for p in f['points']))


def measure(out, water_steps=(8, 16, 32), hydraulic_steps=(4, 8, 16)):
    report = dict(schema='wet-joint-transient-water-diagnosis/1', starting_commit=s.BASE,
                  C04='未達', passed=False, A04_complete=False, C01_complete=False, C02_complete=False, C03_complete=False,
                  NC01_complete=False, P6_complete=False, numerical_acceptance=False,
                  mechanical_load_steps_fixed=read(out/'identity.json')['mechanical_steps'],
                  hydrostatic_water_steps=list(water_steps), Richards_physical_steps_per_second=list(hydraulic_steps),
                  mechanical_audits={}, mechanical_native_vs_independent={}, mechanical_water_step_sensitivity={},
                  hydraulic_audits={}, hydraulic_native_vs_independent={}, hydraulic_time_step_sensitivity={})
    for n in water_steps:
        for kind in ('native', 'independent'):
            report['mechanical_audits'][f'{kind}-w{n}'] = audit(out, out/f'{kind}-w{n}')
        report['mechanical_native_vs_independent'][str(n)] = pair(out/f'native-w{n}', out/f'independent-w{n}')
    for a, b in zip(water_steps, water_steps[1:]):
        for kind in ('native', 'independent'):
            report['mechanical_water_step_sensitivity'][f'{kind}/{a}__{b}'] = pair(out/f'{kind}-w{a}', out/f'{kind}-w{b}')
    for n in hydraulic_steps:
        p = out/f'hydraulic-h{n}'
        report['hydraulic_audits'][str(n)] = hydraulic_audit(p)
        report['hydraulic_native_vs_independent'][str(n)] = hydraulic_pair(out, p, p, 'native', 'independent')
    for a, b in zip(hydraulic_steps, hydraulic_steps[1:]):
        for kind in ('native', 'independent'):
            report['hydraulic_time_step_sensitivity'][f'{kind}/{a}__{b}'] = hydraulic_pair(out, out/f'hydraulic-h{a}', out/f'hydraulic-h{b}', kind, kind)
    report['local_contact_diagnostic'] = local_contact(out, out/f'native-w{water_steps[-1]}')
    report['refusals'] = read(out/'scope-refusals/verification.json')
    report['mechanical_restarts'] = read(out/'mechanical-restarts/verification.json')
    report['hydraulic_restarts'] = read(out/'hydraulic-restarts/verification.json')
    report['observed_vs_whole_run_exact_except_elapsed'] = s.strip_clock(read(out/f'native-w{water_steps[-1]}/raw/run.json')) == s.strip_clock(read(out/'whole/raw/run.json'))
    assert report['observed_vs_whole_run_exact_except_elapsed']
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--output', type=Path, required=True); p.add_argument('--report', type=Path, required=True)
    p.add_argument('--water-steps', nargs='+', type=int, default=[8, 16, 32]); p.add_argument('--hydraulic-steps', nargs='+', type=int, default=[4, 8, 16])
    a = p.parse_args(); write(a.report, measure(a.output.resolve(), a.water_steps, a.hydraulic_steps))


if __name__ == '__main__': main()
