"""Independent load, domain, mixed-equation and output audit for projected MCC.

This accepts a completed numerical run, not spatial convergence or a calibrated
soil model. Material histories and P1 work-conjugate stress stay distinct.
"""
import copy
import json
from pathlib import Path
import meshio
import numpy as np
import fem_reference as ref
import mixed_mcc_io
from native_workflow_snapshot import sha
from verify_high_order import physical, prepare, internal as material_internal
from verify_geostatic import verify_monitors


def require(condition, message):
    if not condition:
        raise ValueError(message)


def restricted_pressure(mesh, fields, active):
    """Independent consistent P1 mass projection over the remaining domain.

    Matrix-free diagonal-preconditioned CG avoids a quadratic-size dense matrix.
    Material stress histories are inputs; no new constitutive update is made.
    """
    material = mesh['material_ids']
    keys = sorted({(node,material[e]) for e,tet in enumerate(mesh['tetrahedra']) if active[e] for node in tet})
    lookup = {key:i for i,key in enumerate(keys)}
    owners, matrices = [], []
    rhs = np.zeros(len(keys))
    diagonal = np.zeros(len(keys))
    cursor = 0
    mean = fields['raw']['material_mean_effective_stress_compression_positive_kpa']
    for e,tet in enumerate(mesh['tetrahedra']):
        q = ref.volume(mesh,e)
        count = len(q['weight'])
        if active[e]:
            nodes = [lookup[(node,material[e])] for node in tet]
            shape = np.column_stack((1-q['natural'].sum(axis=1),q['natural']))
            mass = shape.T@(q['weight'][:,None]*shape)
            np.add.at(rhs,nodes,shape.T@(q['weight']*mean[cursor:cursor+count]))
            np.add.at(diagonal,nodes,np.diag(mass))
            owners.append(nodes)
            matrices.append(mass)
        cursor += count
    indices = np.asarray(owners)
    matrices = np.asarray(matrices)
    require(np.all(diagonal>0), 'singular restricted pressure mass matrix')
    def multiply(value):
        out = np.zeros(len(keys))
        np.add.at(out,indices,np.einsum('eab,eb->ea',matrices,value[indices]))
        return out
    x = np.zeros(len(keys))
    r = rhs.copy()
    z = r/diagonal
    direction = z.copy()
    rho = float(r@z)
    tolerance = 1e-12*max(np.linalg.norm(rhs),1.)
    for iteration in range(2000):
        if np.linalg.norm(rhs-multiply(x)) <= tolerance:
            break
        product = multiply(direction)
        denominator = float(direction@product)
        require(np.isfinite(denominator) and denominator>0, 'restricted pressure projection lost positive definiteness')
        alpha = rho/denominator
        x += alpha*direction
        r -= alpha*product
        z = r/diagonal
        next_rho = float(r@z)
        direction = z+(next_rho/rho)*direction
        rho = next_rho
    else:
        raise ValueError('restricted pressure projection failed its true residual tolerance')
    result = np.zeros(len(fields['owner']))
    cursor = 0
    for e,tet in enumerate(mesh['tetrahedra']):
        q = ref.volume(mesh,e)
        count = len(q['weight'])
        if active[e]:
            shape = np.column_stack((1-q['natural'].sum(axis=1),q['natural']))
            result[cursor:cursor+count] = shape@x[[lookup[(node,material[e])] for node in tet]]
        cursor += count
    return result


def internal(mesh, fields, active, restrict=False):
    result = np.zeros((len(mesh['nodes']), 3))
    stress = fields['raw']['work_conjugate_effective_stress_compression_positive_kpa']
    if restrict:
        pressure = restricted_pressure(mesh,fields,active)
        stress = stress.copy()
        stress[:,:3] += (pressure-stress[:,:3].mean(axis=1))[:,None]
    cursor = 0
    for element in range(len(active)):
        quadrature = ref.volume(mesh, element)
        count = len(quadrature['weight'])
        if active[element]:
            value = -np.einsum('qai,qa,q->i', quadrature['b'], stress[cursor:cursor+count], quadrature['weight']).reshape(-1, 3)
            np.add.at(result, quadrature['nodes'], value)
        cursor += count
    return result


def initial_state(mesh, ground, case, run):
    first = run['frames'][0]['state']
    require(all(first['active_elements']), 'projected initial domain must be fully active')
    initialization = ground['initialization']
    if initialization['method'] == 'horizontal_k0':
        require(run.get('geostatic') is None, 'unexpected geostatic history on K0 input')
        require(not np.asarray(first['displacement_m']).any(), 'K0 initial displacement is nonzero')
        cursor = 0
        for element in range(len(mesh['tetrahedra'])):
            q = ref.volume(mesh, element)
            stress, pore, material, layer = ref.initial_fields(mesh, ground, element, q['position'])
            require(not pore.any(), 'projected K0 input must be dry')
            count = len(q['weight'])
            points = first['material_states'][cursor:cursor+count]
            require(np.allclose([p['effective_stress_compression_positive_kpa'] for p in points], stress, rtol=1e-12, atol=1e-10), 'K0 stress differs from layer weights and K0 ratios')
            mean = stress[:,:3].mean(axis=1)
            consolidation = layer['consolidation']
            pc = consolidation['ratio']*mean if consolidation['definition'] == 'mcc_pc_over_mean_effective' else np.full(count, consolidation['preconsolidation_pressure_kpa'])
            params = material['mcc']
            void = params['specific_volume_ncl_at_reference']-params['lambda']*np.log(pc/params['reference_pressure_kpa'])+params['kappa']*np.log(pc/mean)-1
            require(np.allclose([p['preconsolidation_pressure_kpa'] for p in points], pc, rtol=1e-12, atol=1e-10), 'K0 pc differs from the explicit consolidation definition')
            require(np.allclose([p['void_ratio'] for p in points], void, rtol=1e-12, atol=1e-12), 'K0 void ratio differs from the MCC reference line')
            require(not np.asarray([p['plastic_strain_compression_positive'] for p in points]).any(), 'K0 initial plastic strain is nonzero')
            cursor += count
    else:
        geo = run.get('geostatic')
        require(geo is not None and geo['complete'] and geo['progress'] == 1. and geo['failure_reason'] is None, 'incomplete geostatic history')
        seed, last = geo['seed'], geo['last_committed']
        cfg = initialization['geostatic']
        definition = cfg['seed']
        if definition['kind'] == 'body_states':
            regions = {r['id']:r['body'] for r in mesh['regions']}
            expected = [definition['states'][regions[mesh['region_ids'][p['element_index']]]] for p in run['integration_points']]
        else:
            require(definition['mesh_sha256'] == run['provenance']['mesh_sha256'], 'geostatic seed mesh differs')
            require(len(definition['points']) == len(run['integration_points']), 'geostatic seed inventory differs')
            expected = []
            for point, layout in zip(definition['points'],run['integration_points']):
                require(point['element_index'] == layout['element_index'] and point['local_index'] == layout['local_index'] and
                    point['material_id'] == mesh['material_ids'][point['element_index']] and
                    np.allclose(point['position_m'],layout['position_m'],rtol=0.,atol=1e-8), 'geostatic seed point binding differs')
                expected.append(point['state'])
        require(expected == seed['material_states'], 'geostatic seed differs from input material history')
        require(all(seed['active_elements']) and not np.asarray(seed['displacement_m']).any(), 'geostatic seed domain/displacement differs')
        # On the unchanged seed domain, the saved P1 field is authoritative.
        # Audit its weak projection independently rather than replacing it by
        # another solve, whose finite residual changes near-zero force entries.
        seed_force = np.asarray(seed['applied_nodal_load_kn']).reshape(-1,3)
        fields = mixed_mcc_io.frame_fields(mesh,run,dict(state=seed, relative_force_residual=0.,
            support_reaction_kn=[0.,0.,0.],applied_force_kn=seed_force.sum(axis=0)))
        seed_work = internal(mesh,fields,np.ones(len(mesh['tetrahedra']),bool))
        base,_ = physical(mesh,ground,np.ones(len(mesh['tetrahedra']),bool),initialization['surface_surcharge_kpa'],case.get('pressure_on_original_surface_only',False))
        require(np.allclose(seed_work.reshape(-1),seed['applied_nodal_load_kn'],rtol=1e-10,atol=1e-8), 'geostatic seed holding force differs from independent P1 projection')
        require(np.allclose((seed_work-base).reshape(-1),geo['initial_holding_load_kn'],rtol=1e-10,atol=1e-8), 'geostatic initial release differs')
        for key in ('displacement_m','material_states','mixed_mcc','active_elements','applied_nodal_load_kn'):
            require(last[key] == first[key], 'first frame differs from geostatic committed state: '+key)
        require(not np.asarray(last['geostatic_release_load_kn']).any(), 'geostatic holding force remains at equilibrium')
        require(run['reference_displacement_m'] == first['displacement_m'], 'geostatic display origin differs')
        require(all(a['committed_state_unchanged_on_rejection'] for a in geo['attempts'] if not a['accepted']), 'geostatic rejection changed history')


def audit(mesh, ground, case, run):
    require(run.get('schema') == mixed_mcc_io.PROJECTED_SCHEMA, 'projected MCC result required')
    require(ref.normalize_ground(copy.deepcopy(ground)) == run['ground'], 'result ground differs from input')
    normalized = ref.normalize_stage_case(copy.deepcopy(case))
    require(normalized == run['case'], 'result construction differs from input')
    require(not any(normalized.get(name) for name in ('joint_system', 'beam_system', 'anchor_system', 'initially_inactive_bodies')),
        'unsupported projected construction coupling')
    mixed_mcc_io.validate(mesh, run)
    case, ground = run['case'], run['ground']
    initial_state(mesh, ground, case, run)
    quadrature, offsets = prepare(mesh,run['integration_points'])
    xyz = np.asarray(mesh['nodes'])
    reference = np.asarray(run.get('reference_displacement_m', np.zeros_like(xyz)))
    active = np.ones(len(mesh['tetrahedra']), bool)
    factor = ground['initialization'].get('geostatic', {}).get('initial_boundary_factor', 0.)
    pressure = ground['initialization']['surface_surcharge_kpa']
    regions = {entry['id']:entry['body'] for entry in mesh['regions']}
    body_elements = {}
    for index, region in enumerate(mesh['region_ids']):
        body_elements.setdefault(regions[region], []).append(index)
    prior_fields = None
    rows = []
    expected_audits = []
    require(run['frames'][0]['name'] == ('initial_equilibrated' if run.get('geostatic') is not None else 'initial_k0'), 'initial frame origin differs')
    for index, frame in enumerate(run['frames']):
        state = frame['state']
        before = active.copy()
        selected = []
        if index:
            stage = case['stages'][index-1]
            require(frame['name'] == stage['name'], 'frame name differs from construction stage')
            if stage['kind'] == 'boundary':
                factor = stage['target_factor']
            elif stage['kind'] == 'pressure':
                pressure = stage['target_pressure_kpa']
            elif stage['kind'] == 'excavate':
                for body in stage['bodies']:
                    elements = mesh.get('excavation_sets', {}).get(body, body_elements.get(body))
                    require(elements is not None and len(elements) > 0, 'unknown or empty excavation body')
                    require(all(type(e) is int and 0 <= e < len(active) and before[e] for e in elements), 'excavation overlaps removed soil')
                    require(any(active[e] for e in elements), 'excavation body is already selected')
                    selected.extend(elements)
                    active[elements] = False
                selected = sorted(set(selected))
                pressure = stage['remaining_surface_pressure_kpa']
            elif stage['kind'] != 'reset_display':
                raise ValueError('unsupported projected stage')
            if stage['kind'] == 'reset_display':
                from verify_construction import verify_display_reset
                verify_display_reset(run,index)
        require(np.array_equal(state['active_elements'], active), 'active domain differs from requested excavation')
        require(state['boundary_factor'] == factor and state['surface_pressure_kpa'] == pressure, 'committed load factor/pressure differs')
        require(state['stage_index'] == index and state['stage_progress'] == (1. if index else 0.), 'nonfinal construction frame')
        force, geometry = physical(mesh, ground, active, pressure, case.get('pressure_on_original_surface_only', False))
        if case.get('boundary') and factor:
            force += factor*ref.conditions_load(mesh, case['boundary'], active)
        require(np.allclose(force.reshape(-1), state['applied_nodal_load_kn'], rtol=1e-10, atol=1e-8), 'applied load differs from independent geometry/load assembly')
        require(np.linalg.norm(state['excavation_release_load_kn']) < 1e-9, 'completed stage retains excavation holding load')
        fields = mixed_mcc_io.frame_fields(mesh, run, frame)
        _,maximum_yield = material_internal(mesh,ground,state['material_states'],active,quadrature,offsets)
        require(np.isclose(maximum_yield,frame['maximum_normalized_yield'],rtol=1e-8,atol=1e-8), 'material yield diagnostic differs')
        work = internal(mesh, fields, active)
        conditions = {**case, 'boundary':case.get('boundary') or {'constraints':[]}}
        bases, constraint_error = ref.staged_constraint_space(mesh, conditions, state)
        reaction = np.zeros_like(xyz)
        for node, basis in bases.items():
            residual = work[node]-force[node]
            reaction[node] = residual-basis@(basis.T@residual)
        require(abs(constraint_error-frame['maximum_constraint_error_m']) < 1e-11, 'constraint diagnostic differs')
        balance = (force+reaction).sum(axis=0)
        require(np.allclose(balance, frame['force_balance_error_kn'], rtol=1e-8, atol=1e-7), 'force balance diagnostic differs')
        center = (xyz.min(axis=0)+xyz.max(axis=0))/2
        moment = np.cross(xyz-center, force+reaction).sum(axis=0)
        relative_moment = float(np.linalg.norm(moment)/max(np.linalg.norm(force, axis=1).sum()*np.linalg.norm(np.ptp(xyz, axis=0)), 1.))
        require(relative_moment < 1e-7, 'independent global moment balance failed')
        display = np.asarray(state.get('display_reference_displacement_m') or reference)
        previous = run['frames'][max(0,index-1)]
        verify_monitors(mesh, case, frame, previous, display)
        if selected:
            used = np.unique(ref.cells(mesh)[active])
            holding = internal(mesh, prior_fields, active, restrict=True)-force
            holding[np.setdiff1d(np.arange(len(xyz)), used)] = 0.
            matches = [a for a in run['excavation_audits'] if a['stage'] == frame['name']]
            require(len(matches) == 1, 'excavation audit missing or duplicated')
            recorded = matches[0]
            require(recorded['removed_element_count'] == len(selected), 'removed element count differs')
            volume = sum(ref.volume(mesh, e)['weight'].sum() for e in selected)
            require(np.isclose(volume, recorded['removed_volume_m3'], rtol=1e-10, atol=1e-10), 'removed volume differs')
            require(np.allclose(holding.reshape(-1), recorded['initial_holding_load_kn'], rtol=1e-9, atol=1e-7), 'excavation holding differs from restricted work-conjugate stress')
            require(recorded['zero_release_relative_force_residual'] < case['settings']['relative_force_tolerance'], 'initial excavation residual exceeds tolerance')
            expected_audits.append(frame['name'])
        rows.append(dict(stage=frame['name'], active_elements=int(active.sum()),
            independent_force_residual=fields['metrics']['independent_work_conjugate_force_residual'],
            relative_moment_error=relative_moment, force_balance_error_kn=balance.tolist(),
            **geometry, **fields['metrics']))
        prior_fields = fields
    require([a['stage'] for a in run['excavation_audits']] == expected_audits, 'unexpected excavation audits')
    require(all(a['committed_state_unchanged_on_rejection'] for a in run['attempts'] if not a['accepted']), 'rejected trial changed committed state')
    return rows


def verify(mesh_path, ground_path, result_path, export_dir, solver=None):
    require(not any((result_path.parent/name).exists() for name in ('stage-verification.json','accepted-stage-run.json')),
        'new projected audit receipts required; existing accepted records are never overwritten')
    paths = [mesh_path, ground_path, result_path, result_path.parent/'stages.json']
    binding = {str(p.resolve()):sha(p) for p in paths}
    mesh, ground, run, case = [json.loads(p.read_text(encoding='utf-8-sig')) for p in paths]
    rows = audit(mesh, ground, case, run)
    metadata = json.loads((export_dir/'export-verification.json').read_text(encoding='utf-8'))
    require(metadata.get('schema') == 'soil-fem-stage-export/4', 'projected VTU export required')
    require(metadata['mesh_sha256'] == sha(mesh_path) and metadata['result_sha256'] == sha(result_path), 'export source differs')
    for name, digest in metadata['sha256'].items():
        path = (export_dir/name).resolve()
        require(path.is_relative_to(export_dir.resolve()) and sha(path) == digest, 'VTU output changed')
    for index, frame in enumerate(run['frames']):
        fields = mixed_mcc_io.frame_fields(mesh, run, frame)
        gp = meshio.read(export_dir/f'integration-points-{index:03}.vtu')
        active = fields['active'][fields['owner']]
        require(np.array_equal(gp.point_data['integration_point_id'], np.flatnonzero(active)), 'exported GP identity differs')
        for name, values in fields['raw'].items():
            require(name in gp.point_data and np.array_equal(gp.point_data[name], values[active]), 'exported primary/work-conjugate/volume field differs: '+name)
        material = frame['state']['material_states']
        for name in ('effective_stress_compression_positive_kpa', 'plastic_strain_compression_positive', 'void_ratio', 'preconsolidation_pressure_kpa'):
            expected = np.asarray([s[name] for s in material])[active]
            require(np.array_equal(gp.point_data[name], expected), 'exported material history differs: '+name)
    binary = solver or Path(__file__).resolve().parents[1]/'target/release/soil-fem-3d.exe'
    require(run['provenance']['solver_binary_sha256'] == sha(binary), 'result produced by a different solver')
    require(binding == {str(p.resolve()):sha(p) for p in paths}, 'input changed during audit')
    report = dict(schema='soil-fem-stage-verification/3', frames=rows, vtu_roundtrip='passed', raw_integration_point_roundtrip='passed',
        formulation=mixed_mcc_io.PROJECTED_FORMULATION, independent_load_and_domain_audit='passed',
        spatial_convergence_verified=False)
    package = dict(schema='soil-fem-accepted-staged-run/1', inputs_and_results_sha256=binding,
        solver_binary_sha256=run['provenance']['solver_binary_sha256'])
    for name, value in [('stage-verification.json', report), ('accepted-stage-run.json', package)]:
        with (result_path.parent/name).open('x', encoding='utf-8') as stream:
            json.dump(value, stream, indent=2, allow_nan=False)
    print(json.dumps(dict(projected_stage_audit='passed', frames=len(rows))), flush=True)
    return report
