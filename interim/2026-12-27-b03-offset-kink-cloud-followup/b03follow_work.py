"""Action-surface accounting of sealed states; zero new equilibrium solves."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys

for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_key] = '1'
import numpy as np
from b03follow_math import arm_map, dot, mv, norm, scatter, write
from b03follow_rotation import query

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
_path = sys.path[:]
try:
    sys.path.insert(0, str(ROOT/'tools'))
    import fem_reference as ref
    from plate6_offset_reference import Element, System
finally:
    sys.path[:] = _path


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def element_account(kp, a, d, reference=None, load=None):
    d = np.asarray(d)
    reference = np.zeros_like(d) if reference is None else np.asarray(reference)
    relative = d-reference
    dp, relp = mv(a, d), mv(a, relative)
    qp = mv(kp, relp)
    qs = mv(a.T, qp)
    t = np.array([i for i in range(len(d)) if i % 6 < 3])
    r = np.array([i for i in range(len(d)) if i % 6 >= 3])
    trans, intrinsic = dot(qs[t], d[t]), dot(qp[r], d[r])
    couple = dot(qs[r]-qp[r], d[r])
    twice_energy = dot(relp, qp)
    report = dict(neutral_internal_absolute_work_kn_m=dot(dp, qp),
        attachment_translation_internal_work_kn_m=trans,
        neutral_intrinsic_moment_internal_work_kn_m=intrinsic,
        attachment_arm_couple_internal_work_kn_m=couple,
        attachment_total_internal_work_kn_m=dot(d, qs),
        twice_elastic_energy_kn_m=twice_energy,
        installation_reference_internal_work_kn_m=dot(mv(a, reference), qp),
        surface_mapping_work_residual_kn_m=dot(d, qs)-dot(dp, qp),
        reference_identity_residual_kn_m=dot(d, qs)-twice_energy-dot(mv(a, reference), qp))
    if load is not None:
        load = np.asarray(load)
        # Per-panel external force IS present for these coupled snapshots.
        # Invert A^T explicitly; do not invent an allocation of the free patch's nodal seed.
        fp = load.copy()
        fp[r] -= mv((a-np.eye(len(d))).T, load)[r]
        report.update(attachment_external_work_kn_m=dot(load, d),
            neutral_external_work_kn_m=dot(fp, dp),
            attachment_external_translation_work_kn_m=dot(load[t], d[t]),
            neutral_external_intrinsic_moment_work_kn_m=dot(fp[r], d[r]),
            attachment_external_arm_couple_work_kn_m=dot(load[r]-fp[r], d[r]),
            external_surface_mapping_residual_kn_m=dot(load, d)-dot(fp, dp))
    return report, qs


def free_accounts(prior, output):
    source = prior/'b03cause-final-patches'
    requests, responses = read(source/'requests.json'), read(source/'responses.json')
    by_case = {}
    for request, response in zip(requests, responses, strict=True):
        name = request['id'].rsplit('/', 1)[0].replace('/', '-')
        by_case.setdefault(name, []).append((request, response))
    rows = []
    for side, folder in (('native', 'assemblies'), ('independent', 'independent-assemblies')):
        for path in sorted((source/folder).glob('*.json')):
            saved = read(path)
            d = np.array(saved['global_displacement'])
            f = np.array(saved['global_force'])
            parts, panels = [], []
            for ei, (element, (request, response)) in enumerate(zip(saved['elements'], by_case[path.stem], strict=True)):
                dofs = [6*n+a for n in element['nodes'] for a in range(6)]
                points = np.array(request['points'])
                normal = np.cross(points[1]-points[0], points[2]-points[0])
                normal /= norm(normal)
                a = arm_map(request['e']*normal)
                kp = np.asarray(response['neutral_stiffness']) if side == 'native' else Element(points, request['section'], request['e']).neutral.stiffness
                account, internal = element_account(kp, a, d[dofs])
                account['element_index'] = ei
                panels.append(account)
                parts.append((dofs, internal))
            internal = scatter(parts, len(d))
            twice_saved_energy = dot(d, mv(saved['stiffness'], d))
            couple = math.fsum(p['attachment_arm_couple_internal_work_kn_m'] for p in panels)
            neutral = math.fsum(p['neutral_internal_absolute_work_kn_m'] for p in panels)
            trans = dot(f.reshape(-1, 6)[:, :3], d.reshape(-1, 6)[:, :3])
            moments = dot(f.reshape(-1, 6)[:, 3:], d.reshape(-1, 6)[:, 3:])
            rows.append(dict(id=path.stem, side=side, panels=panels,
                external_attachment_translation_work_kn_m=trans,
                external_attachment_moment_work_kn_m=moments,
                external_attachment_total_work_kn_m=dot(f, d),
                twice_saved_matrix_energy_kn_m=twice_saved_energy,
                neutral_internal_total_work_kn_m=neutral,
                attachment_arm_couple_internal_work_kn_m=couple,
                force_dot_d_minus_twice_u_kn_m=dot(f, d)-twice_saved_energy,
                full_surface_account_residual_kn_m=dot(f, d)-neutral,
                omitted_arm_account_residual_kn_m=dot(f, d)-(neutral-couple),
                element_scatter_minus_saved_internal_l2=norm(internal-saved['global_internal_force']),
                soil_side=dict(status='not_applicable', reason='no soil in this saved free assembly', work_kn_m=None),
                support_reaction=dict(status='not_applicable', physical_reaction_kn_m=0.,
                    numerical_residual_work_kn_m=dot(internal-f, d)),
                external_per_panel_neutral_work_kn_m=None,
                external_allocation_reason='Manufactured global nodal seed has no saved per-panel external allocation; rejected folded arms have no unique common neutral point.',
                source_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    write(output/'free-accounts.json', rows)
    return rows


def soil_account(mesh, run, state):
    u = np.asarray(state['displacement_m']).ravel()
    cells = ref.cells(mesh)
    parts, energy_parts = [], []
    by_element = {}
    for index, point in enumerate(run['integration_points']):
        by_element.setdefault(point['element_index'], []).append((index, point))
    for ei, cell in enumerate(cells):
        if not state['active_elements'][ei]:
            continue
        dofs = [3*n+a for n in cell for a in range(3)]
        positions = np.asarray(mesh['nodes'])[cell]
        terms = [[] for _ in dofs]
        for index, point in sorted(by_element[ei], key=lambda item: item[1]['local_index']):
            material = state['material_states'][index]
            if any(material['plastic_strain_compression_positive']) or material['plastic_work_kpa'] != 0.:
                raise ValueError('elastic energy accounting unavailable for plastic history')
            stress = -np.asarray(material['effective_stress_compression_positive_kpa'])
            _, gradient = ref.basis(np.array([point['natural_coordinates']]), len(cell) == 10)
            j = (positions-positions[0]).T@gradient[0]
            g = gradient[0]@np.linalg.inv(j)
            b = np.zeros((6, len(dofs)))
            for i, (gx, gy, gz) in enumerate(g):
                b[0, 3*i] = gx; b[1, 3*i+1] = gy; b[2, 3*i+2] = gz
                b[3, 3*i] = gy; b[3, 3*i+1] = gx
                b[4, 3*i+1] = gz; b[4, 3*i+2] = gy
                b[5, 3*i] = gz; b[5, 3*i+2] = gx
            weight = point['weight_m3']
            for col in range(len(dofs)):
                terms[col].append(dot(b[:, col], stress)*weight)
            energy_parts.append(dot(stress, mv(b, u[dofs]))*weight)
        parts.append((dofs, np.array([math.fsum(v) for v in terms])))
    return scatter(parts, len(u)), math.fsum(energy_parts)


def coupled_accounts(prior, output, probe):
    rows = []
    for name in ('planar', 'kink-zero', 'normal-sign-pair', 'step-split', 'kink-split', 'height-step'):
        folder = prior/'b03cause-final-controls/cli'/name
        mesh, case, run = (read(folder/(k+'.json')) for k in ('mesh', 'stages', 'run'))
        if not run['complete'] or any(any(v) for v in (s['effective_stress_compression_positive_kpa'] for s in run['ground']['initialization']['geostatic']['seed']['states'].values())):
            raise ValueError('completed dry zero-stress-seed model required for this energy definition')
        model = System(mesh, case['beam_system'])
        native_requests, native_keys = [], []
        for frame in run['frames']:
            state = frame['state']; bs = state['beams']
            for i, (nodes, section, kernel) in enumerate(model.plate_data):
                if not bs['active_plates'][i]:
                    continue
                plate = model.plates[i]
                traction = np.asarray(plate.get('constant_traction_kpa', [0.]*3))+state['boundary_factor']*np.asarray(plate.get('boundary_traction_kpa', [0.]*3))
                traction[2] -= section['unit_weight_kn_m3']*section['thickness_m']
                native_keys.append((frame['name'], i))
                native_requests.append(dict(id=name+'/'+frame['name']+'/'+str(i), points=model.positions[nodes].tolist(),
                    section=section, e=plate.get('attachment_offset_m', 0.), traction=traction.tolist(), fields=[]))
        native_results = dict(zip(native_keys, query(probe, output/'coupled-kernels', name, native_requests), strict=True))
        for operator_side in ('native', 'independent'):
            for frame in run['frames']:
                state = frame['state']; bs = state['beams']
                u = np.asarray(state['displacement_m']).ravel()
                d = np.column_stack((bs['node_displacement_m'], bs['node_rotation_rad'])).ravel()
                internal_parts, load_parts, panels = [], [], []
                for i, (nodes, section, kernel) in enumerate(model.plate_data):
                    if not bs['active_plates'][i]:
                        continue
                    plate = model.plates[i]
                    dofs = [6*n+a for n in nodes for a in range(6)]
                    traction = np.asarray(plate.get('constant_traction_kpa', [0.]*3))+state['boundary_factor']*np.asarray(plate.get('boundary_traction_kpa', [0.]*3))
                    traction[2] -= section['unit_weight_kn_m3']*section['thickness_m']
                    result = native_results[(frame['name'], i)] if operator_side == 'native' else None
                    load = np.asarray(result['load']) if result is not None else kernel.load(traction)
                    kp = np.asarray(result['neutral_stiffness']) if result is not None else kernel.neutral.stiffness
                    account, internal = element_account(kp, kernel.map, d[dofs], bs['plate_reference_global_kinematics'][i], load)
                    account['id'] = plate['id']
                    panels.append(account)
                    internal_parts.append((dofs, internal)); load_parts.append((dofs, load))
                if model.members or case['beam_system'].get('nodal_loads'):
                    raise ValueError('this sealed six-case set is plate-only; do not omit additional members/loads')
                q, f = scatter(internal_parts, len(d)), scatter(load_parts, len(d))
                soil_internal, soil_two_u = soil_account(mesh, run, state)
                mapped_internal, mapped_external = [], []
                for i, n in enumerate(model.attach):
                    if n is not None:
                        mapped_internal.append(([3*n+a for a in range(3)], q[6*i:6*i+3]))
                        mapped_external.append(([3*n+a for a in range(3)], f[6*i:6*i+3]))
                qi, fi = scatter(mapped_internal, len(u)), scatter(mapped_external, len(u))
                soil_external = np.asarray(state['applied_nodal_load_kn'])-fi
                residual_soil = soil_internal+qi-state['applied_nodal_load_kn']
                bases, _ = ref.staged_constraint_space(mesh, case, state)
                reaction = np.zeros(len(u))
                for n, basis in sorted(bases.items()):
                    r = residual_soil[3*n:3*n+3]
                    reaction[3*n:3*n+3] = r-mv(basis, mv(basis.T, r))
                soil_work, structure_work = dot(soil_external, u), dot(f, d)
                internal_work = dot(soil_internal, u)+dot(q, d)
                structure_two_u = math.fsum(p['twice_elastic_energy_kn_m'] for p in panels)
                reference_work = math.fsum(p['installation_reference_internal_work_kn_m'] for p in panels)
                reaction_work = dot(reaction, u)
                # Internal structural work at shared translations cancels the soil-side transfer.
                interface_structure = dot(qi, u)
                rows.append(dict(id=name+'/'+frame['name'], operator_side=operator_side, panels=panels,
                    soil_external_work_kn_m=soil_work, structure_external_work_kn_m=structure_work,
                    soil_internal_work_kn_m=dot(soil_internal, u),
                    structure_internal_work_kn_m=dot(q, d),
                    soil_twice_elastic_energy_kn_m=soil_two_u,
                    structure_twice_elastic_energy_kn_m=structure_two_u,
                    installation_reference_internal_work_kn_m=reference_work,
                    soil_support_reaction_work_kn_m=reaction_work,
                    soil_support_reaction_kn=[math.fsum(float(v) for v in np.asarray(reaction).reshape(-1, 3)[:, axis]) for axis in range(3)],
                    sealed_support_reaction_kn=frame['support_reaction_kn'],
                    soil_to_structure_internal_interface_work_kn_m=interface_structure,
                    structure_to_soil_internal_interface_work_kn_m=-interface_structure,
                    attachment_constraint_work_difference_kn_m=math.fsum(dot(q[6*i:6*i+3], d[6*i:6*i+3]-u[3*n:3*n+3]) for i, n in enumerate(model.attach) if n is not None),
                    external_total_kn_m=soil_work+structure_work,
                    internal_absolute_total_kn_m=internal_work,
                    force_dot_d_minus_twice_u_kn_m=soil_work+structure_work-soil_two_u-structure_two_u,
                    full_reference_and_reaction_account_residual_kn_m=soil_work+structure_work+reaction_work-soil_two_u-structure_two_u-reference_work,
                    force_equilibrium_free_residual_l2_kn=norm(residual_soil-reaction),
                    structural_rotation_residual_l2_kn_m=norm((q-f).reshape(-1, 6)[:, 3:]),
                    cumulative_work_recomputed=False, new_solves=0,
                    source_run_sha256=hashlib.sha256((folder/'run.json').read_bytes()).hexdigest()))
    write(output/'coupled-accounts.json', rows)
    return rows


def run(prior, output, probe):
    output.mkdir(parents=True, exist_ok=False)
    free = free_accounts(prior, output)
    coupled = coupled_accounts(prior, output, probe)
    report = dict(passed=False, b03_complete=False, diagnostic_only=True, production_scope_expanded=False,
        new_equilibrium_solves=0, free_saved_native=23, free_saved_independent=23,
        coupled_saved_runs=6, coupled_unique_stage_end_states=len(coupled)//2, coupled_operator_accounts=len(coupled), cumulative_s11_work_recomputed=False,
        free_maximum_full_account_residual_kn_m=max(abs(r['full_surface_account_residual_kn_m']) for r in free),
        free_maximum_omitted_arm_account_residual_kn_m=max(abs(r['omitted_arm_account_residual_kn_m']) for r in free),
        free_maximum_mapping_residual_kn_m=max(abs(p['surface_mapping_work_residual_kn_m']) for r in free for p in r['panels']),
        coupled_maximum_naive_force_dot_d_minus_twice_u_kn_m=max(abs(r['force_dot_d_minus_twice_u_kn_m']) for r in coupled),
        coupled_maximum_reference_work_kn_m=max(abs(r['installation_reference_internal_work_kn_m']) for r in coupled),
        coupled_maximum_full_account_residual_kn_m=max(abs(r['full_reference_and_reaction_account_residual_kn_m']) for r in coupled),
        coupled_maximum_constraint_work_difference_kn_m=max(abs(r['attachment_constraint_work_difference_kn_m']) for r in coupled),
        free_external_face_allocation_identifiable=False,
        rejected_junction_physical_work_primary_cause=None,
        accepted_account_does_not_certify_undefined_junction=True)
    write(output/'summary.json', report)
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prior-extracted', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--probe', type=Path, required=True)
    a = p.parse_args()
    run(a.prior_extracted.resolve(), a.output.resolve(), a.probe.resolve())
