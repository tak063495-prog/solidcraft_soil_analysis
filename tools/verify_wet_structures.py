"""Small, input-generated 3D elastic wet JOINT + plate + beam reference.

This is a verification driver, not a production solver. It refuses plasticity,
opening/sliding, domain changes and unsupported water/stage kinds. All frame
unknowns are retained with the soil; no Rust tangent, solved field or placement
reference generates a reference solution. Existing independent element modules
are reused without changing their acceptance thresholds.
"""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import shutil

for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_key] = '1'
import numpy as np
import fem_reference as fem
from beam_reference import System
from joint_water_reference import Water
import darcy_reference
from verify_checkpoint import invoke
from verify_high_order import physical

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def without_time(value):
    if isinstance(value, dict):
        return {k: without_time(v) for k, v in value.items() if k != 'elapsed_seconds'}
    if isinstance(value, list):
        return [without_time(v) for v in value]
    return value


def require_closed_stick_history(state):
    """The elastic reference has no plastic contact history to certify.

    These are exact branch invariants, not fitted numerical tolerances. Work
    balance alone cannot establish a status or exclude a small false history.
    """
    for point in state['joints']['points']:
        assert point['status'] == 'stick', 'non-stick status outside elastic wet reference'
        np.testing.assert_array_equal(point['plastic_slip_m'], np.zeros(3),
                                      err_msg='plastic slip outside elastic wet reference')
        assert point['accumulated_plastic_slip_m'] == 0, 'plastic history outside elastic wet reference'
        assert point['plastic_dissipation_kpa_m'] == 0, 'plastic dissipation outside elastic wet reference'


def reference(mesh, ground, case, order=5, water_modes=('continuous',), fixed_prescribed=False):
    """Generate every reference frame exclusively from the three inputs."""
    allowed = {'install_beams', 'remove_beams', 'boundary',
               'hydrostatic_water', 'steady_seepage'}
    assert all(s['kind'] in allowed for s in case['stages']), 'unsupported stage'
    assert not case.get('anchor_system') and not case.get('initially_inactive_bodies')
    assert ground['initialization']['method'] == 'equilibrate'
    cfg = ground['initialization']['geostatic']
    assert cfg['seed']['kind'] == 'body_states' and cfg['initial_boundary_factor'] == 0
    assert case['settings']['kinematics'] == 'explicit'
    assert all(c['value_m'] == 0 and not np.asarray(c['gradient']).any()
               for c in case['boundary']['constraints'])
    assert set(water_modes) <= {'continuous', 'prescribed', 'soil'}
    assert all(j['water']['kind'] in water_modes
               for j in case['joint_system']['joints']), 'water mode requires explicit reference scope'
    if fixed_prescribed:
        assert water_modes == ('prescribed',), 'fixed head requires explicit PRESCRIBED scope'
        level = ground['initialization']['water_level_m']
        assert level is not None
        assert all(j['water'].get('level_m') == level for j in case['joint_system']['joints']), 'fixed head differs from initial reservoir'
        assert all(s['target_level_m'] == level for s in case['stages']
                   if s['kind'] == 'hydrostatic_water'), 'hydrostatic level conflicts with fixed reservoir'
    else:
        assert all(j['water'].get('level_m') is None for j in case['joint_system']['joints']
                   if j['water']['kind'] == 'prescribed'), 'fixed prescribed head requires separate reference'
    roles = set(case['joint_system'].get('structure_bodies', []))
    if roles:
        assert set(water_modes) == {'soil'}, 'structure roles require explicit SOIL scope'
        assert roles <= {r['body'] for r in mesh['regions']}, 'unknown structure body'
    else:
        assert 'soil' not in water_modes, 'SOIL requires an impermeable structure role'
    xyz = np.asarray(mesh['nodes']); ndof = xyz.size
    assert ndof <= 6500, 'small direct reference only'
    active = np.ones(len(mesh['tetrahedra']), dtype=bool)
    water = Water(mesh, ground, case); contact = water.contact
    for point in water.faces:
        policy = water.definitions[point['joint_id']]['water']
        if policy['kind'] == 'soil':
            side = policy['soil_side']; other = 'b' if side == 'a' else 'a'
            assert water.soil[point['owner_'+side]] and not water.soil[point['owner_'+other]], 'SOIL side/role mismatch'
    beams = System(mesh, case['beam_system']); state = beams.initial_state()
    materials = {m['id']: m['mohr_coulomb'] for m in ground['materials']}
    regions = {r['id']: r['body'] for r in mesh['regions']}
    operators = []; k = np.zeros((ndof, ndof))
    for ei, mid in enumerate(mesh['material_ids']):
        m = materials[mid]; d = fem.elastic_tensor(m['young_kpa'], m['poisson'])
        seed = cfg['seed']['states'][regions[mesh['region_ids'][ei]]]
        assert seed['equivalent_plastic_strain'] == 0
        q = fem.volume(mesh, ei, independent=True, independent_order=order)
        dofs = (q['nodes'][:, None]*3+np.arange(3)).ravel()
        k[np.ix_(dofs, dofs)] += np.einsum(
            'qai,ab,qbj,q->ij', q['b'], d, q['b'], q['weight'], optimize=True)
        operators.append((q, dofs, d, np.asarray(seed['effective_stress_compression_positive_kpa']), seed['void_ratio']))
    # Closed/sticking candidate: acceptance below independently checks that
    # this branch is physically admissible, rather than presuming it is so.
    constitutive = np.asarray([kt*np.eye(3)+(kn-kt)*np.outer(n, n)
                              for kt, kn, n in zip(contact.kt, contact.kn, contact.normal)])
    kj = np.einsum('qai,qab,qbj,q->ij', contact.c, constitutive,
                   contact.c, contact.weight, optimize=True)
    basis, error = fem.staged_constraint_space(mesh, case, dict(
        active_elements=active, boundary_factor=0., displacement_m=np.zeros_like(xyz)))
    assert error == 0
    count = sum(b.shape[1] for b in basis.values()); t = np.zeros((ndof, count)); start = 0
    for n, b in basis.items():
        t[3*n:3*n+3, start:start+b.shape[1]] = b; start += b.shape[1]
    u = np.zeros_like(xyz); factor = 0.; level = ground['initialization']['water_level_m']
    gamma = water.gamma; field = None; frames = []
    for stage in [None]+case['stages']:
        darcy = None
        if stage is not None:
            kind = stage['kind']
            if kind in ('install_beams', 'remove_beams'):
                state = beams.change(state, u, stage)
            elif kind == 'boundary':
                factor = stage['target_factor']
            elif kind == 'hydrostatic_water':
                level = stage['target_level_m']; field = None
            else:
                level = stage['exterior_water_level_m']
                darcy = water.steady(stage['flow'], active, mesh['material_ids'], level)
                field = darcy['field']
        interface = water.hydrostatic(active, level) if field is None else water.from_heads(active, field, level)
        seed_force = np.zeros(ndof)
        for ei, (q, dofs, _, s0, _) in enumerate(operators):
            if not water.soil[ei]:
                pore = np.zeros(len(q['weight']))
            elif field is None:
                assert level >= xyz[:, 2].max(), 'partial volume wetting requires separate quadrature convergence'
                pore = gamma*(level-q['position'][:, 2])
            else:
                nodes = q['nodes']; h = field['relative_head'][nodes]
                pore = gamma*(field['datum'][nodes[0]]-q['position'][:, 2]+h[0]+q['shape']@(h-h[0]))
                assert pore.min() >= 0, 'unsaturated reference refused'
            total = np.broadcast_to(s0, (len(pore), 6)).copy(); total[:, :3] += pore[:, None]
            seed_force[dofs] -= np.einsum('qai,qa,q->i', q['b'], total, q['weight'])
        frame_ground = copy.deepcopy(ground); frame_ground['initialization']['water_level_m'] = level
        load, _ = physical(mesh, frame_ground, active, ground['initialization']['surface_surcharge_kpa'],
                           case.get('pressure_on_original_surface_only', False), joint_system=contact)
        load += factor*fem.conditions_load(mesh, case['boundary'], active)+interface['load']
        kb, rest, fb, eff, extf, extm, info = beams.assemble(state, factor)
        offset, binding, qb, normals = beams.space(state, u, soil_basis=t)
        kk = np.block([[t.T@(k+kj)@t+binding.T@kb@binding, binding.T@kb@qb],
                       [qb.T@kb@binding, qb.T@kb@qb]])
        assert np.linalg.matrix_rank(kk) == len(kk), 'independent combined system has a mechanism'
        assert np.linalg.eigvalsh(kk).min() > 0, 'independent stiffness is not positive definite'
        rhsb = fb-rest-kb@offset
        rhs = np.r_[t.T@(load.ravel()-seed_force)+binding.T@rhsb, qb.T@rhsb]
        solution = np.linalg.solve(kk, rhs)
        residual = float(np.linalg.norm(kk@solution-rhs)/max(np.linalg.norm(rhs), 1.))
        assert residual < 1e-10
        u = (t@solution[:count]).reshape(-1, 3)
        br = beams.results(state, offset+binding@solution[:count]+qb@solution[count:],
                           kb, rest, fb, eff, extf, extm, info, normals)
        state = br['state']
        jump = np.einsum('qag,g->qa', contact.c, u.ravel())
        traction = np.einsum('qab,qb->qa', constitutive, jump)
        pressure = -np.einsum('qa,qa->q', traction, contact.normal)
        shear = traction+pressure[:, None]*contact.normal
        assert pressure.min() > 1e-6, 'open/no-pressure contact is outside this reference'
        assert np.all(np.linalg.norm(shear, axis=1) < contact.mu*pressure-1e-6), 'sliding reference refused'
        stress = []; void = []; pore_values = []
        from mc_reference import normalized_yield
        for ei, (_, _, d, s0, e0) in enumerate(operators):
            q = fem.volume(mesh, ei)
            strain = np.einsum('qai,i->qa', q['b'], u[q['nodes']].ravel())
            s = s0-strain@d.T
            assert np.max(normalized_yield(materials[mesh['material_ids'][ei]], s)) < 0
            stress.extend(s); void.extend((1+e0)*np.exp(strain[:, :3].sum(axis=1))-1)
            if not water.soil[ei]:
                pore_values.extend(np.zeros(len(q['weight'])))
            elif field is None:
                pore_values.extend(gamma*(level-q['position'][:, 2]))
            else:
                h = field['relative_head'][q['nodes']]
                pore_values.extend(gamma*(field['datum'][q['nodes'][0]]-q['position'][:, 2]+h[0]+q['shape']@(h-h[0])))
        balance = ((k+kj)@u.ravel()+seed_force-load.ravel()).reshape(-1, 3)
        for n, value in br['interface_residual'].items():
            balance[n] += value
        support = np.zeros_like(xyz)
        for n, b in basis.items():
            support[n] = balance[n]-b@(b.T@balance[n])
        moment = np.cross(xyz, load+support).sum(axis=0)+br['external_moment']+br['support_moment']
        assert np.linalg.norm(moment) < 1e-6
        frames.append(dict(u=u.copy(), beams=copy.deepcopy(state), results=br,
                           stress=np.asarray(stress), void=np.asarray(void), traction=traction,
                           pore=np.asarray(pore_values),
                           pressure=pressure, interface=interface, head=field, darcy=darcy, load=load,
                           support=support.sum(axis=0)+br['support_force'], moment=moment,
                           residual=residual))
    return frames


def compare(mesh, ground, case, native, expected):
    assert native['complete'] and native['schema'] == 'soil-fem-stage-result/14'
    assert native['ground'] == ground and native['case'] == case
    assert len(native['frames']) == len(expected)
    assert [f['name'] for f in native['frames']] == ['initial_equilibrated']+[s['name'] for s in case['stages']]
    flow_indices = {i for i, s in enumerate(case['stages']) if s['kind'] == 'steady_seepage'}
    assert len(native['flow_audits']) == len(flow_indices)
    assert {a['stage_index'] for a in native['flow_audits']} == flow_indices
    water = Water(mesh, ground, case, native['joint_integration_points']); rows = []
    previous_actual = previous_ref = None
    for frame, r in zip(native['frames'], expected, strict=True):
        s = frame['state']; actual = np.asarray(s['displacement_m'])
        require_closed_stick_history(s)
        du = float(np.linalg.norm(actual-r['u'])/max(np.linalg.norm(r['u']), 1e-12))
        assert du < 1e-7
        increment = 0.
        if previous_actual is not None:
            inc = r['u']-previous_ref
            increment = float(np.linalg.norm(actual-previous_actual-inc)/max(np.linalg.norm(inc), 1e-12))
            assert increment < 1e-7
        previous_actual = actual; previous_ref = r['u']
        for key in ('active_members', 'active_plates'):
            assert s['beams'][key] == r['beams'][key]
        for key in ('node_displacement_m', 'node_rotation_rad', 'reference_global_kinematics', 'plate_reference_global_kinematics'):
            np.testing.assert_allclose(s['beams'][key], r['beams'][key], rtol=0, atol=1e-9)
        points = s['material_states']; assert all(p['equivalent_plastic_strain'] == 0 and p['plastic_work_kpa'] == 0 for p in points)
        stress = np.asarray([p['effective_stress_compression_positive_kpa'] for p in points])
        ds = float(np.linalg.norm(stress-r['stress'])/np.linalg.norm(r['stress'])); assert ds < 1e-7
        np.testing.assert_allclose([p['void_ratio'] for p in points], r['void'], rtol=0, atol=1e-8)
        if s.get('pore_pressure_kpa') is not None:
            np.testing.assert_allclose(s['pore_pressure_kpa'], r['pore'], rtol=1e-10, atol=1e-8)
        np.testing.assert_allclose([p['traction_kpa'] for p in s['joints']['points']], r['traction'], rtol=1e-7, atol=1e-7)
        assert all(p['accumulated_plastic_slip_m'] == 0 for p in s['joints']['points'])
        water.audit(s, frame['joint_results'], r['interface'])
        applied = r['load'].copy()
        for n, value in r['results']['attached_load'].items():
            applied[n] += value
        np.testing.assert_allclose(s['applied_nodal_load_kn'], applied.ravel(), rtol=1e-9, atol=1e-8)
        np.testing.assert_allclose(frame['support_reaction_kn'], r['support'], rtol=1e-9, atol=1e-5)
        np.testing.assert_allclose(frame['applied_force_kn'], r['load'].sum(axis=0)+r['results']['external_force'], rtol=1e-9, atol=1e-8)
        np.testing.assert_allclose(frame['global_moment_balance_error_kn_m'], r['moment'], rtol=0, atol=1e-6)
        assert frame['relative_force_residual'] <= case['settings']['relative_force_tolerance']
        assert frame['joint_support']['unrestrained_modes'] == 0
        assert {m['id'] for m in frame['beam_results']['members']} == set(r['results']['members'])
        for member in frame['beam_results']['members']:
            for key, value in r['results']['members'][member['id']].items():
                np.testing.assert_allclose(member[key], value, rtol=1e-7, atol=1e-6)
        assert len(frame['beam_results'].get('plates', [])) == len(r['results']['plates'])
        assert {p['id'] for p in frame['beam_results'].get('plates', [])} == set(r['results']['plates'])
        for plate in frame['beam_results'].get('plates', []):
            predicted = r['results']['plates'][plate['id']]
            for key in ('area_m2', 'local_axes', 'membrane_energy_kn_m', 'bending_energy_kn_m', 'drilling_penalty_energy_kn_m'):
                np.testing.assert_allclose(plate[key], predicted[key], rtol=1e-7, atol=1e-8)
            for point, refpoint in zip(plate['samples'], predicted['samples'], strict=True):
                for key in refpoint:
                    np.testing.assert_allclose(point[key], refpoint[key], rtol=1e-7, atol=1e-8)
        rows.append(dict(name=frame['name'],relative_displacement_error=du,
                         relative_increment_error=increment,relative_stress_error=ds,
                         reference_true_residual=r['residual'],minimum_effective_contact_pressure_kpa=float(r['pressure'].min())))
    for flow in native['flow_audits']:
        ref = expected[flow['stage_index']+1]
        field = ref['head']
        np.testing.assert_array_equal(flow['joint_flow']['original_field']['defined_nodes'], field['defined'])
        if set(case['joint_system'].get('structure_bodies', [])) or any(
                j['water']['kind'] == 'prescribed' and j['water'].get('level_m') is not None
                for j in case['joint_system']['joints']):
            np.testing.assert_array_equal(flow['solution']['active_elements'], water.soil)
            hydraulic = ref['darcy']
            assert len(flow['solution']['points']) == sum(len(q['weight']) for q in hydraulic['solution']['quadrature'].values())
            rows[flow['stage_index']+1]['independent_darcy'] = darcy_reference.verify(
                hydraulic['hydraulic_mesh'], hydraulic['case'], flow['solution'], mesh['material_ids'])
        np.testing.assert_allclose(flow['joint_flow']['original_field']['head_m'], field['head'], rtol=0, atol=1e-8)
    return rows


def run(solver, output, water_mode=None):
    output.mkdir(parents=True, exist_ok=False)
    inputs = []
    for src, dest in [('joint-water-column-tet10-mesh.json', 'mesh.json'),
                      ('ground-wet-plate-beam.json', 'ground.json'), ('stages-wet-plate-beam.json', 'stages.json')]:
        path = output/dest; shutil.copyfile(ROOT/'examples'/src, path); inputs.append(path)
    fixed_prescribed = water_mode == 'prescribed-fixed'
    policy_mode = 'prescribed' if fixed_prescribed else ('soil' if water_mode == 'soil-source' else water_mode)
    water_modes = ('continuous',) if policy_mode is None else (policy_mode,)
    assert water_mode in (None, 'prescribed', 'prescribed-fixed', 'soil', 'soil-source')
    policy_source = None
    if water_mode is not None:
        policy_source = ROOT/'examples'/('stages-joint-water-'+policy_mode+'.json')
        source_case = read(policy_source)
        policy = source_case['joint_system']['joints'][0]['water']
        case_input = read(inputs[2])
        assert len(case_input['joint_system']['joints']) == 1
        case_input['joint_system']['joints'][0]['water'] = policy
        if fixed_prescribed:
            # Existing initial reservoir is held fixed; the incompatible rise
            # is omitted, not reinterpreted as a new hydraulic policy.
            level = read(inputs[1])['initialization']['water_level_m']
            policy['level_m'] = level
            case_input['stages'] = [s for s in case_input['stages']
                                   if s['kind'] != 'hydrostatic_water' or s['target_level_m'] == level]
        if policy_mode == 'soil':
            # Carry the existing SOIL continuum role AND its single soil-side
            # bottom-head flow boundary; the structural top is not an aquifer.
            case_input['joint_system']['structure_bodies'] = source_case['joint_system']['structure_bodies']
            flow_source = next(s for s in source_case['stages'] if s['kind'] == 'steady_seepage')
            for stage in case_input['stages']:
                if stage['kind'] == 'steady_seepage':
                    stage['flow'] = copy.deepcopy(flow_source['flow'])
                    stage['exterior_water_level_m'] = flow_source['exterior_water_level_m']
                    if water_mode == 'soil-source':
                        # Existing Darcy source law (tests/darcy.rs quadratic control),
                        # with the inherited SOIL material conductivity. No new geometry/BC.
                        for material in stage['flow']['materials']:
                            material['volume_source_per_s'] = 2*material['hydraulic_conductivity_m_s'][2][2]
        inputs[2].write_text(json.dumps(case_input, indent=2)+'\n', encoding='utf-8')
    sources = [ROOT/'tools'/n for n in ('verify_wet_structures.py', 'beam_reference.py', 'plate_reference.py',
               'joint_water_reference.py', 'joint_reference.py', 'fem_reference.py', 'darcy_reference.py',
               'hydrostatic_reference.py', 'verify_high_order.py', 'mc_reference.py', 'verify_checkpoint.py',
               'export_stages.py', 'high_order_io.py', 'beam_io.py', 'plate_io.py', 'joint_io.py')]
    if policy_source is not None:
        sources.append(policy_source)
    if water_mode == 'soil-source':
        sources.append(ROOT/'tools/verify_soil_source_column.py')
    fingerprints = {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p): sha(p) for p in [solver]+inputs+sources}
    rp = output/'stage-results.json'
    invoke([str(solver), 'solve-stages', *map(str, inputs), str(rp)], output/'solve.log')
    mesh = read(inputs[0]); ground = fem.normalize_ground(read(inputs[1])); case = fem.normalize_stage_case(read(inputs[2])); native = read(rp)
    assert native['provenance']['solver_binary_sha256'] == sha(solver)
    expected = reference(mesh, ground, case, water_modes=water_modes, fixed_prescribed=fixed_prescribed); rows = compare(mesh, ground, case, native, expected)
    from export_stages import export
    export(inputs[0], rp, output/'exports')
    exported = read(output/'exports/export-verification.json')
    assert exported['result_sha256'] == sha(rp) and exported['mesh_sha256'] == sha(inputs[0])
    assert all(exported[key] == 'passed' for key in ('vtu_roundtrip', 'integration_point_roundtrip',
               'joint_point_roundtrip', 'plate_sample_roundtrip'))
    if policy_mode == 'soil':
        import meshio
        point_soil = np.asarray([Water(mesh, ground, case).soil[p['element_index']]
                                 for p in native['integration_points']])
        for index, ref in enumerate(expected):
            gp = meshio.read(output/'exports'/f'integration-points-{index:03}.vtu').point_data
            ids = gp['integration_point_id']
            np.testing.assert_array_equal(gp['material_stress_reference'], (~point_soil[ids]).astype(np.uint8))
            np.testing.assert_allclose(gp['pore_pressure_kpa'], ref['pore'][ids], rtol=1e-10, atol=1e-8)
            total = ref['stress'][ids].copy(); total[:, :3] += ref['pore'][ids, None]
            np.testing.assert_allclose(gp['total_stress_compression_positive_kpa'], total, rtol=1e-7, atol=1e-7)
    fine = reference(mesh, ground, case, order=7, water_modes=water_modes, fixed_prescribed=fixed_prescribed)
    refinement = max(float(np.linalg.norm(a['u']-b['u'])/max(np.linalg.norm(b['u']), 1e-12)) for a, b in zip(expected, fine, strict=True))
    assert refinement < 1e-7
    # Fresh clean accepted path catches leakage from the deliberate rejections.
    clean = read(inputs[2]); selected = clean['verification_reject_once_stages']; clean['verification_reject_once_stages'] = []
    rejected = [a for a in native['attempts'] if not a['accepted']]
    assert len(rejected) == len(selected)
    assert all(a['committed_state_unchanged_on_rejection'] and a['rejected_trial_structure_changed'] and a['rejected_trial_material_changed'] for a in rejected)
    for stage in clean['stages']:
        if stage['name'] in selected:
            stage['initial_increment'] = (stage.get('initial_increment') or clean['settings']['initial_increment'])/2
    cp = output/'clean-stages.json'; cp.write_text(json.dumps(clean), encoding='utf-8')
    cr = output/'clean-result.json'
    invoke([str(solver), 'solve-stages', *map(str, inputs[:2]), str(cp), str(cr)], output/'clean.log')
    clean_run = read(cr)
    assert clean_run['complete'] and clean_run['last_committed'] == native['last_committed'] and clean_run['frames'] == native['frames']
    assert without_time(clean_run.get('flow_audits')) == without_time(native.get('flow_audits'))
    checkpoints = []; accepted = sum(a['accepted'] for a in native['geostatic']['attempts']); seen = set()
    restart_points = [('geostatic', 1)]
    for attempt in native['attempts']:
        if not attempt['accepted']: continue
        accepted += 1; name = attempt['stage']
        if name in seen: continue
        seen.add(name); restart_points.append((name, accepted))
    for index, (name, accepted) in enumerate(restart_points):
        checkpoint = output/f'checkpoint-{index}.json'; resumed = output/f'resumed-{index}.json'
        invoke([str(solver), 'checkpoint-stages', *map(str, inputs), str(accepted), str(checkpoint)], output/f'checkpoint-{index}.log')
        saved = read(checkpoint); assert not saved['run']['complete']
        progress = saved['run']['geostatic']['progress'] if name == 'geostatic' else saved['run']['last_committed']['stage_progress']
        assert 0 < progress < 1
        invoke([str(solver), 'resume-stages', *map(str, inputs), str(checkpoint), str(resumed)], output/f'resume-{index}.log')
        restored = read(resumed); assert restored['complete']
        for key in ('last_committed', 'frames', 'flow_audits', 'excavation_audits', 'transition_audits'):
            assert without_time(restored.get(key)) == without_time(native.get(key)), (name, key)
        # Only measured linear-solver timings are excluded.
        assert without_time(restored['attempts']) == without_time(native['attempts'])
        checkpoints.append(dict(stage=name, accepted_count=accepted, progress=progress, exact=True))
    assert fingerprints == {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p): sha(p) for p in [solver]+inputs+sources}
    report = dict(schema='soil-fem-wet-structures-verification/1', passed=True, frames=rows,
                  independent_gl5_to_gl7_relative_difference=refinement,
                  clean_rejected_path_exact=True, checkpoints=checkpoints, inputs_sources_solver_sha256=fingerprints,
                  result_sha256=sha(rp), output_roundtrip=exported, p6_complete=False, p3b_p4_rows_complete=False,
                  water_modes=list(water_modes), fixed_prescribed=fixed_prescribed,
                  soil_volume_source=water_mode == 'soil-source',
                  scope='12 Tet10, 2 linear plates, 1 linear beam; hydrostatic and steady Darcy '+str(water_modes)+' JOINT, elastic soil, closed sticking contact; full retained direct 3D reference. No plastic, Richards, BIRTH, crack-tip or general interface convergence acceptance.')
    if water_mode == 'soil-source':
        from verify_soil_source_column import audit
        report['analytic_soil_source'] = audit(mesh, ground, case, native)
    (output/'verification.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--water-mode', choices=('prescribed', 'prescribed-fixed', 'soil', 'soil-source'))
    a = p.parse_args()
    print(json.dumps(run(a.solver.resolve(), a.output.resolve(), a.water_mode), indent=2))
