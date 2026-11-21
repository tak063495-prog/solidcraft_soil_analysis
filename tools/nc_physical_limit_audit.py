"""Read-only load and all-active-GP audit of an existing common-lambda NC fold.

This tool does not integrate a material, solve equilibrium, alter a saved state,
or certify NC-01. Cone formulas are diagnostic surrogates, not MCC bounds.
"""
import numpy as np


def stress_values(state, m):
    s = np.asarray(state['effective_stress_compression_positive_kpa'], dtype=float)
    if s.shape != (6,) or not np.all(np.isfinite(s)) or not np.isfinite(m) or m <= 0.:
        raise ValueError('finite six-component stress and positive M required')
    p = float(s[:3].mean())
    q = float(np.sqrt(((s[0]-s[1])**2+(s[1]-s[2])**2+(s[2]-s[0])**2)/2 + 3*np.dot(s[3:], s[3:])))
    pc = state['preconsolidation_pressure_kpa']
    plastic = np.asarray(state['plastic_strain_compression_positive'], dtype=float)
    if p <= 0. or not np.isfinite(pc) or pc <= 0. or plastic.shape != (6,) or not np.all(np.isfinite(plastic)):
        raise ValueError('positive p and pc and finite six-component plastic strain required')
    evp = float(sum(plastic[:3]))
    return dict(sigma_v_kpa=float(s[2]), sigma_h_kpa=float(s[:2].mean()),
                sigma_x_kpa=float(s[0]), sigma_y_kpa=float(s[1]),
                shear_kpa=s[3:].tolist(), p_kpa=p, q_kpa=q,
                q_over_Mp=q/(m*p), pc_kpa=pc,
                stress_side='dry' if p < pc/2 else 'wet' if p > pc/2 else 'critical',
                plastic_volume_strain_compression_positive=evp,
                normalized_ellipse_yield=((q/pc)**2+m*m*(p/pc)*(p/pc-1))/(m*m),
                critical_state_m=m)


def sign(value):
    return 'positive' if value > 0 else 'negative' if value < 0 else 'zero'


def cone_lambda(m, sh, sv0, gravity_vertical=0.):
    required = sh*(3-2*m)/(3+m)
    return dict(required_vertical_kpa=required, sigma_h_kpa=sh, sigma_v0_kpa=sv0,
                self_weight_vertical_kpa=gravity_vertical,
                lambda_uniform=1-required/sv0,
                lambda_self_weight_corrected=1-(required-gravity_vertical)/(sv0-gravity_vertical),
                exact_bound_for_mcc=False)


def read(path):
    import json
    return json.loads(path.read_text(encoding='utf-8'))


def write(path, value):
    import json
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')


def resultant(xyz, force):
    force = np.asarray(force).reshape(-1, 3)
    return dict(force_kn=force.sum(axis=0).tolist(),
                moment_about_origin_kn_m=np.cross(xyz, force).sum(axis=0).tolist(),
                norm_kn=float(np.linalg.norm(force)))


def face_label(xyz):
    for axis, name in enumerate(('x', 'y', 'z')):
        for value in (0., .5, 1.):
            if np.all(np.abs(xyz[:, axis]-value) < 1e-10):
                return f'{name}={value:g}'
    return 'other'


def load_audit(geometry, loads, origin):
    import fem_reference as fem
    from joint_reference import System
    from verify_high_order import prepare, physical, internal
    from nc_fold_path import native_loads
    mesh, ground, case = (geometry[k] for k in ('mesh', 'ground', 'case'))
    xyz = np.asarray(mesh['nodes']); active = np.asarray(geometry['active_elements'], dtype=bool)
    actual_final = np.asarray(loads['final_load_kn']).reshape(-1, 3)
    actual_h = np.asarray(loads['holding_load_kn']).reshape(-1, 3)
    final, holding = native_loads(mesh, ground, case, origin, active)
    q, offsets = prepare(mesh); joint = System(mesh, case['joint_system'])
    retained = dict(origin, active_elements=active.tolist())
    removed_force, _ = internal(mesh, ground, origin['material_states'],
                                np.asarray(origin['active_elements']) & ~active, q, offsets)
    origin_external = np.asarray(origin['applied_nodal_load_kn']).reshape(-1, 3)
    contact_lost = joint.audit(origin)-joint.audit(retained)
    components = []; body = np.zeros_like(xyz)
    materials = {m['id']: m for m in ground['materials']}
    for ei in np.flatnonzero(active):
        qp = fem.volume(mesh, ei, independent=True)
        f = np.zeros_like(xyz)
        f[qp['nodes'], 2] -= (qp['weight'] @ qp['shape'])*materials[mesh['material_ids'][ei]]['unit_weight_kn_m3']
        body += f
        components.append(dict(kind='self_weight', element=int(ei), **resultant(xyz, f)))
    surface = np.zeros_like(xyz); owners = fem.topology(mesh); face_forces = {}
    base = ground['initialization']['geostatic']['base_loads']
    for li, load in enumerate(base.get('surface_loads', [])):
        for face in mesh['face_sets'][load['face_set']]:
            key = tuple(sorted(face)); live = [ei for ei in owners[key] if active[ei]]
            if not live:
                continue
            qp = fem.face(mesh, live[0], face)
            temporary = dict(mesh, face_sets={'one': [face]})
            force = fem.conditions_load(temporary, dict(surface_loads=[dict(load, face_set='one')]), active, True)
            surface += force
            face_forces[key] = face_forces.get(key, np.zeros_like(xyz))+force
            components.append(dict(kind='surface_load', source_index=li, input_load=load,
                vertices=list(face), nodes=qp['nodes'].tolist(), plane=face_label(xyz[list(face)]),
                area_m2=float(np.linalg.norm(qp['area_vector'], axis=1).sum()),
                nodal_force_kn=[dict(node=int(n), force_kn=force[n].tolist()) for n in qp['nodes']],
                **resultant(xyz, force)))
    for li, load in enumerate(base.get('nodal_loads', [])):
        f = fem.conditions_load(mesh, dict(nodal_loads=[load]), active)
        surface += f; components.append(dict(kind='nodal_load', source_index=li, **resultant(xyz, f)))
    physical_final = physical(mesh, ground, active, 0., True, joint_system=joint)[0]
    used = np.unique(fem.cells(mesh)[active]); live_nodes = set(used.tolist())
    r = np.asarray(geometry['transform']); free = np.linalg.norm(r, axis=1).reshape(-1, 3) > 0.
    nodes = [dict(node=i, position_m=x.tolist(), active=i in live_nodes,
        free_components=free[i].tolist(), final_kn=actual_final[i].tolist(), holding_kn=actual_h[i].tolist(),
        at_lambda_zero_kn=(actual_final[i]+actual_h[i]).tolist(),
        origin_external_kn=origin_external[i].tolist(), removed_internal_kn=removed_force[i].tolist(),
        lost_contact_internal_kn=contact_lost[i].tolist(), self_weight_kn=body[i].tolist(),
        retained_surface_kn=surface[i].tolist()) for i, x in enumerate(xyz)]
    faces = []
    for face, elements in owners.items():
        live = [ei for ei in elements if active[ei]]
        if len(live) != 1:
            continue
        opposite = joint.counterparts.get(face)
        if opposite is not None and active[opposite]:
            continue
        qp = fem.face(mesh, live[0], face); ids = qp['nodes']; force = face_forces.get(face, np.zeros_like(xyz))
        faces.append(dict(vertices=list(face), nodes=ids.tolist(), plane=face_label(xyz[list(face)]),
            is_cut_seam=opposite is not None or len(elements) == 2,
            area_m2=float(np.linalg.norm(qp['area_vector'], axis=1).sum()),
            final_surface_traction_resultant_kn=force.sum(axis=0).tolist(),
            nodal_final_sum_kn=actual_final[ids].sum(axis=0).tolist(),
            nodal_holding_sum_kn=actual_h[ids].sum(axis=0).tolist()))
    groups = {}
    for face in faces:
        name = face['plane']; groups.setdefault(name, dict(area_m2=0., final_surface_traction_kn=[0., 0., 0.], faces=0))
        groups[name]['area_m2'] += face['area_m2']; groups[name]['faces'] += 1
        groups[name]['final_surface_traction_kn'] = (np.asarray(groups[name]['final_surface_traction_kn'])+face['final_surface_traction_resultant_kn']).tolist()
    return dict(schema='nc-physical-load-audit/1', nodes=nodes, faces=faces, additive_physical_components=components,
        face_groups=groups, boundary_constraints=case['boundary']['constraints'],
        kinematics=case['settings']['kinematics'], active_elements=int(active.sum()), active_nodes=len(used),
        final=resultant(xyz, actual_final), holding=resultant(xyz, actual_h),
        at_lambda_zero=resultant(xyz, actual_final+actual_h), self_weight=resultant(xyz, body),
        retained_surface=resultant(xyz, surface),
        original_holding_reconstruction_max_abs_kn=float(np.max(np.abs(actual_h.ravel()-holding))),
        original_final_reconstruction_max_abs_kn=float(np.max(np.abs(actual_final.ravel()-final))),
        physical_component_sum_max_abs_kn=float(np.max(np.abs(body+surface-physical_final))),
        physical_vs_actual_final_max_abs_kn=float(np.max(np.abs(physical_final-actual_final))),
        preload_stage_factor_at_removal=origin['boundary_factor'],
        remaining_surface_pressure_kpa=case['stages'][2]['remaining_surface_pressure_kpa'],
        surface_pressure_on_new_seam_kn=[0., 0., 0.],
        note='Face traction resultants are additive; nodal sums on faces overlap at edges and include self weight and adjacent-face loads. H is an assembled nodal vector, not a unique traction on each face.')


def gp_audit(geometry, initial, samples, origin):
    from verify_high_order import prepare
    mesh, ground = geometry['mesh'], geometry['ground']
    q, offsets = prepare(mesh); active = np.asarray(geometry['active_elements'], dtype=bool)
    params = {m['id']: m for m in ground['materials']}
    peak_index = max(range(len(samples)), key=lambda i: samples[i]['progress'])
    records = []
    for si in (peak_index,):
        sample = samples[si]; previous = samples[si-1] if si else initial
        for ei in np.flatnonzero(active):
            material = params[mesh['material_ids'][ei]]; m = material['mcc']['critical_state_m']
            for qi, position in enumerate(q[ei]['position']):
                pi = int(offsets[ei])+qi; state = sample['states'][pi]
                values = stress_values(state, m); old = previous['states'][pi]
                devp = values['plastic_volume_strain_compression_positive']-sum(old['plastic_strain_compression_positive'][:3])
                sv0 = origin['material_states'][pi]['effective_stress_compression_positive_kpa'][2]
                weight = material['unit_weight_kn_m3']*(.5-position[2])
                formula = cone_lambda(m, values['sigma_h_kpa'], sv0, weight)
                records.append(dict(sample_index=si, step=sample['step'], lambda_value=sample['progress'],
                    peak_sample=si == peak_index, point=pi, element=int(ei), local_point=qi,
                    position_m=position.tolist(), weight_m3=float(q[ei]['weight'][qi]), **values,
                    pre_removal=stress_values(origin['material_states'][pi], m),
                    plastic_volume_strain_sign=sign(values['plastic_volume_strain_compression_positive']),
                    plastic_volume_strain_increment=devp, plastic_volume_increment_sign=sign(devp),
                    formula=formula, fold_minus_uniform_formula=sample['progress']-formula['lambda_uniform'],
                    fold_minus_weight_corrected_formula=sample['progress']-formula['lambda_self_weight_corrected']))
    peak = [p for p in records if p['peak_sample']]; top = [p for p in peak if p['position_m'][2] >= .45]
    def summary(rows):
        keys = ('sigma_v_kpa', 'sigma_h_kpa', 'p_kpa', 'q_kpa', 'q_over_Mp', 'pc_kpa', 'normalized_ellipse_yield',
                'plastic_volume_strain_increment', 'fold_minus_uniform_formula', 'fold_minus_weight_corrected_formula')
        return dict(count=len(rows), ranges={k: dict(min=float(min(r[k] for r in rows)),
            median=float(np.median([r[k] for r in rows])), max=float(max(r[k] for r in rows))) for k in keys},
            stress_side_counts={k: sum(r['stress_side'] == k for r in rows) for k in ('dry', 'wet', 'critical')},
            plastic_volume_increment_sign_counts={k: sum(r['plastic_volume_increment_sign'] == k for r in rows) for k in ('positive', 'negative', 'zero')},
            q_greater_than_Mp_count=sum(r['q_over_Mp'] > 1. for r in rows))
    return dict(schema='nc-physical-gp-audit/1', all_activity_gp_records=records, peak_lambda=samples[peak_index]['progress'],
        peak_sample_index=peak_index, neighbour_progress=[r['progress'] for r in samples[max(0, peak_index-2):peak_index+3]], active_gp_count=len(peak), peak_summary=summary(peak), top_z_ge_045_summary=summary(top),
        minimum_effective_mean_pressure_kpa=min(p['p_kpa'] for p in peak),
        uniform_formula=cone_lambda(1.4, 100., 100.),
        formula_limitation='Axisymmetric uniform no-shear cone surrogate, not an MCC nonexistence theorem. GP formulas use actual pre-removal sigma_v0 and peak sigma_h; gamma*(0.5-z) is only the uniform column self-weight correction, not an exact decomposition of this nonuniform FE stress field.',
        stress_side_definition='dry p<pc/2, wet p>pc/2; not a separate material branch or a new yield constraint',
        plastic_sign_convention='compression positive; cumulative and immediately preceding accepted-increment signs are distinct',
        top_selection='z>=0.45 m is a reporting band only; every active GP is retained')


def lift_recorded_mode(geometry, sample, diagnostic):
    """Lift the recorded K singular vector; do not call it an augmented null vector."""
    import fem_reference as fem
    xyz = np.asarray(geometry['mesh']['nodes']); r = np.asarray(geometry['transform'])
    v = np.asarray(diagnostic['right_weak_mode'], dtype=float).copy()
    tangent = np.asarray(sample['tangent'], dtype=float)
    if r.shape != (xyz.size, len(v)) or tangent.shape != (len(v)+1,) or not np.all(np.isfinite(v)):
        raise ValueError('recorded mode or transform shape differs')
    if v@tangent[:-1] < 0.:
        v = -v
    full = (r@v).reshape(-1, 3)
    active = np.asarray(geometry['active_elements'], dtype=bool)
    used = np.unique(fem.cells(geometry['mesh'])[active]); top = used[np.abs(xyz[used, 2]-.5) < 1e-10]
    return dict(schema='nc-recorded-weak-mode/1', lambda_value=sample['progress'],
        mode_source='existing analyze_nc_fold_path right singular vector of actual nonsymmetric K',
        exact_null_vector=False, augmented_right_null_vector_exists=False,
        note='K at the accepted sample is full rank; the bordered system is regular. This is a nearby weak mode, not an exact zero mode or a stability certificate.',
        recorded_stiffness=diagnostic['stiffness'], recorded_augmented=diagnostic['augmented'],
        orientation='aligned to increasing recorded continuation arc length',
        normalization='reduced displacement singular vector has Euclidean norm one; directions, not physical displacement amplitudes',
        nodes=[dict(node=i, position_m=x.tolist(), displacement_direction=full[i].tolist(), active=bool(i in used)) for i, x in enumerate(xyz)],
        top_nodes=top.tolist(), top_vertical_min=float(full[top, 2].min()), top_vertical_max=float(full[top, 2].max()),
        top_vertical_all_positive=bool(np.all(full[top, 2] > 0.)),
        maximum_node=int(np.linalg.norm(full, axis=1).argmax()),
        vertical_norm=float(np.linalg.norm(full[:, 2])), horizontal_norm=float(np.linalg.norm(full[:, :2])),
        dynamic_stability_classified=False, second_order_work_audited=False, acoustic_tensor_audited=False,
        passed=False, nc01_complete=False)


def run(path, origin_checkpoint, output, origin_mesh=None):
    import hashlib, json
    from pathlib import Path
    from nc_fold_path import CONTROL, digest
    manifest = read(path/'evidence-sha256.json')
    for name, expected in manifest.items():
        if hashlib.sha256((path/name).read_bytes()).hexdigest() != expected:
            raise ValueError('fold source evidence checksum differs: '+name)
    controls = read(path/'controls.json')
    if controls['binding']['fold_control'] != CONTROL:
        raise ValueError('common lambda control required')
    geometry = read(path/'geometry.json')
    if geometry['map'] != 'native':
        raise ValueError('native saved-state audit only')
    initial = read(path/'initial.json')
    if digest(initial['payload']) != initial['payload_sha256'] or initial['binding'] != controls['binding']:
        raise ValueError('initial binding or payload differs')
    samples = [json.loads(line) for line in (path/'samples.jsonl').read_text(encoding='utf-8').splitlines()]
    if not samples:
        raise ValueError('no accepted fold samples')
    checkpoint = read(origin_checkpoint); origin = checkpoint['run']['continuation']['origin']
    from nc_fixed_inputs import fixed_inputs_match
    if not fixed_inputs_match((geometry['mesh'], geometry['ground'], geometry['case']),
                             (read(origin_mesh or origin_checkpoint.parent.parent/'input/mesh.json'), checkpoint['run']['ground'], checkpoint['run']['case'])):
        raise ValueError('origin checkpoint inputs differ from fold geometry')
    output.mkdir(parents=True, exist_ok=False)
    load = load_audit(geometry, read(path/'original-load.json'), origin)
    write(output/'load-audit.json', load)
    gp = gp_audit(geometry, initial['payload'], samples, origin)
    write(output/'gp-audit.json', gp)
    audit_path = path.parent/(path.name+'-audit')/'fold-results.json'
    fold = read(audit_path)
    selected = samples[gp['peak_sample_index']]
    diagnostic = next(p['diagnostic'] for p in fold['points'] if p['step'] == selected['step'])
    write(output/'weak-mode.json', lift_recorded_mode(geometry, selected, diagnostic))
    write(output/'source-binding.json', dict(path=str(path), original_checkpoint=str(origin_checkpoint),
        diagnostic_only=True, passed=False, nc01_complete=False,
        input_sha256={n: hashlib.sha256((path/n).read_bytes()).hexdigest() for n in ('geometry.json', 'original-load.json', 'samples.jsonl', 'initial.json', 'evidence-sha256.json')},
        origin_sha256=hashlib.sha256(origin_checkpoint.read_bytes()).hexdigest(),
        tool_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
    return dict(peak=gp['peak_lambda'], gp=gp['peak_summary'], load=load['face_groups'])


def run_partial(folder, output, accepted_rows=None):
    """Audit real accepted rows from a time-bounded trace, without inventing a Run."""
    import hashlib, json
    from probe_nc_midpoint import inputs
    from nc_fixed_inputs import fixed_inputs_match
    mesh, ground, case = [read(folder/(name+'.json')) for name in ('mesh', 'ground', 'stages')]
    expected = inputs('midpoint'); denominator = round(1/case['stages'][2]['initial_increment'])
    if denominator not in (128, 256):
        raise ValueError('original proportional 128 or 256 series required')
    expected[2]['stages'][2]['initial_increment'] = 1/denominator
    if not fixed_inputs_match((mesh, ground, case), expected):
        raise ValueError('partial trace input changes original NC physics or controls')
    path = accepted_rows or folder/'history/accepted.jsonl'
    rows = [json.loads(line) for line in path.read_text(encoding='utf-8').splitlines()]
    checkpoint_path = folder/'history/checkpoint-stage-3.json'
    checkpoint = read(checkpoint_path)
    if not fixed_inputs_match((mesh, ground, case), (mesh, checkpoint['run']['ground'], checkpoint['run']['case'])):
        raise ValueError('trace checkpoint input differs')
    origin = checkpoint['run']['continuation']['origin']
    observed = [(row['accepted_count'], row['state']) for row in rows if row['state']['stage_index'] == 3]
    if len(observed) < 2:
        raise ValueError('no pair of accepted removal states')
    active = observed[-1][1]['active_elements']
    if any(row['active_elements'] != active for _, row in observed):
        raise ValueError('partial removal trace changes active topology')
    samples = [dict(step=i, progress=row['stage_progress'], states=row['material_states']) for i, row in observed]
    geometry = dict(mesh=mesh, ground=ground, active_elements=active)
    gp = gp_audit(geometry, samples[0], samples, origin)
    gp['maximum_available_lambda'] = gp.pop('peak_lambda')
    gp['maximum_available_sample_index'] = gp.pop('peak_sample_index')
    gp['maximum_available_summary'] = gp.pop('peak_summary')
    for row in gp['all_activity_gp_records']:
        row['selected_maximum_available_point'] = row.pop('peak_sample')
        for name in ('uniform', 'weight_corrected'):
            row['sample_minus_'+name+'_formula'] = row.pop('fold_minus_'+name+'_formula')
    gp.update(observation='partial proportional history near the reported PC fold; not a newly observed fold',
              local_fold_observed=False, fold_lambda=None, nc01_complete=False, passed=False,
              process_observation=read(folder/'process-observation.json'),
              source_sha256={str(p.relative_to(folder)): hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in (path, checkpoint_path, folder/'mesh.json', folder/'ground.json', folder/'stages.json')})
    output.mkdir(parents=True, exist_ok=False); write(output/'gp-audit.json', gp)
    return dict(maximum_available_lambda=gp['maximum_available_lambda'], summary=gp['maximum_available_summary'],
                local_fold_observed=False, complete=False)


if __name__ == '__main__':
    import argparse, json
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--path', type=Path)
    mode.add_argument('--partial-history', type=Path)
    parser.add_argument('--origin-checkpoint', type=Path)
    parser.add_argument('--origin-mesh', type=Path)
    parser.add_argument('--accepted-rows', type=Path, help='exact byte slice of actual accepted JSONL rows; does not create a Run or checkpoint')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.path is not None:
        if args.origin_checkpoint is None:
            parser.error('--path requires --origin-checkpoint')
        report = run(args.path, args.origin_checkpoint, args.output, args.origin_mesh)
    else:
        report = run_partial(args.partial_history, args.output, args.accepted_rows)
    print(json.dumps(report, indent=2, ensure_ascii=False))
