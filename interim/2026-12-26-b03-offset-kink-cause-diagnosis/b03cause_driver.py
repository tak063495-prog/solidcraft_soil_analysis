"""Small free structural assemblies; internal counterfactual, never acceptance."""
import argparse
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

for _key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[_key] = '1'
import numpy as np
from b03cause_math import arm_observation, cross, dot, frame, mv, norm, rigid_basis, wrench

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
_old_path = sys.path[:]
try:
    sys.path.insert(0, str(ROOT/'tools'))
    from plate6_offset_reference import Element as IndependentElement
finally:
    sys.path[:] = _old_path
SECTION = dict(id='b03cause', young_kpa=3e7, poisson=.2, thickness_m=.1,
               unit_weight_kn_m3=0., drilling_factor=1e-4)
TRACTION = [1., .5, -.25]


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def panels(angle, offsets, split=False, permutation=None, rotation=None):
    r = math.radians(angle)
    corners = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.],
                        [0., -math.cos(r), math.sin(r)]])
    triangles = [[0, 1, 2], [1, 0, 3]]
    positions, ids, elements = [], {}, []
    for k, tri in enumerate(triangles):
        local = [(str(v), corners[v]) for v in tri]
        for i, j in ((0, 1), (1, 2), (2, 0)):
            edge = sorted((tri[i], tri[j]))
            local.append(('m'+str(edge), (corners[edge[0]]+corners[edge[1]])/2))
        nodes = []
        for name, point in local:
            name = (str(k)+':'+name) if split else name
            if name not in ids:
                ids[name] = len(positions); positions.append(point.tolist())
            nodes.append(ids[name])
        e = offsets[k]
        if permutation is not None and k == 1:
            nodes = [nodes[i] for i in permutation]
        elements.append(dict(nodes=nodes, e=e))
    positions = np.asarray(positions)
    if rotation is not None:
        positions = positions@np.asarray(rotation).T
    return positions, elements


def global_report(positions, elements, results, independent_components):
    n = 6*len(positions); matrix = np.zeros((n, n)); term_lists = {}
    element_rows = []
    rigid = rigid_basis(positions)
    for element, result in zip(elements, results):
        nodes = element['nodes']; dofs = [6*i+a for i in nodes for a in range(6)]
        k = np.asarray(result['stiffness']); kp = np.asarray(result['neutral_stiffness'])
        arm = frame(positions[nodes])*element['e']; a = np.eye(36)
        for i in range(6):
            a[6*i:6*i+3, 6*i+3:6*i+6] = cross(arm)
        # Element insertion order is explicit and every scalar sum is compensated.
        for i, gi in enumerate(dofs):
            for j, gj in enumerate(dofs):
                term_lists.setdefault((gi, gj), []).append(float(k[i, j]))
        trial = np.sin(np.arange(36)+.25)*1e-3
        neutral_d = mv(a, trial)
        # Virtual work of neutral traction and its attachment force/couple.
        f = np.asarray(result['load']); fp = np.zeros(36)
        # Integral of a quadratic Tri6 shape is zero at corners and area/3 at midsides.
        area = norm(np.cross(positions[nodes[1]]-positions[nodes[0]],
                             positions[nodes[2]]-positions[nodes[0]]))/2
        for i in (3, 4, 5):
            fp[6*i:6*i+3] = np.asarray(TRACTION)*area/3
        transformed_load = np.array([dot(a[:, j], fp) for j in range(36)])
        correct_work = dot(f, trial)
        omitted_work = dot(fp, trial)
        elem_rigid = rigid[dofs]
        element_rows.append(dict(arm_m=arm.tolist(), dofs=dofs,
            stiffness_transform_max_abs=float(np.max(np.abs(k-a.T@kp@a))),
            load_transform_max_abs=float(np.max(np.abs(f-transformed_load))),
            virtual_work_residual_kn_m=correct_work-dot(fp, neutral_d),
            omitted_arm_couple_work_error_kn_m=omitted_work-correct_work,
            membrane_bending_arm_energy_difference_kn_m=
                .5*dot(trial, mv(k, trial))-.5*dot(trial, mv(kp, trial)),
            rigid_residuals=[norm(mv(k, elem_rigid[:, i])) for i in range(6)],
            kernel_fields=result['fields']))
    for (i, j), terms in sorted(term_lists.items()):
        matrix[i, j] = math.fsum(terms)
    eigenvalues, vectors = np.linalg.eigh(matrix)
    cutoff = 1e-10*float(eigenvalues[-1])
    null = eigenvalues <= cutoff
    basis = []
    for component in independent_components:
        r = rigid_basis(positions[component]); full = np.zeros((n, 6))
        for i, node in enumerate(component):
            full[6*node:6*node+6] = r[6*i:6*i+6]
        basis.append(full)
    gauge = np.concatenate(basis, axis=1)
    q, _ = np.linalg.qr(gauge, mode='complete')
    v = q[:, gauge.shape[1]:]
    force_seed = np.sin(np.arange(n)+.125)
    # Fixed self-equilibrated manufactured load, with no physical support.
    force = mv(v, mv(v.T, force_seed))
    reduced = v.T@matrix@v
    solution = mv(v, np.linalg.solve(reduced, mv(v.T, force)))
    internal = mv(matrix, solution); residual = internal-force
    trans = np.array([i for i in range(n) if i % 6 < 3]); rot = np.array([i for i in range(n) if i % 6 >= 3])
    work = dot(force, solution); twice_energy = dot(solution, internal)
    maximum_rigid = max(norm(mv(matrix, rigid[:, i]))/(norm(matrix)*norm(rigid[:, i])) for i in range(6))
    raw = dict(positions_m=positions.tolist(), elements=elements, stiffness=matrix.tolist(),
               rigid_basis=rigid.tolist(), eigenvalues=eigenvalues.tolist(),
               global_force=force.tolist(), global_displacement=solution.tolist(),
               global_internal_force=internal.tolist(), global_residual=residual.tolist(),
               gauge_basis=gauge.tolist(), element_observations=element_rows)
    summary = dict(nodes=len(positions), dofs=n, tri6=2, physical_rigid_modes=6,
        component_rigid_modes=gauge.shape[1], inherited_eigen_cutoff=cutoff,
        inherited_algebraic_zero_count=int(np.sum(null)),
        additional_modes_below_inherited_cutoff=int(np.sum(null))-gauge.shape[1],
        global_rigid_normalized_residual=maximum_rigid,
        force_residual_l2_kn=norm(residual[trans]), moment_residual_l2_kn_m=norm(residual[rot]),
        external_wrench=wrench(positions, force), internal_wrench=wrench(positions, internal),
        work_kn_m=work, twice_energy_kn_m=twice_energy, work_balance_residual_kn_m=work-twice_energy,
        maximum_omitted_arm_couple_work_error_kn_m=max(abs(x['omitted_arm_couple_work_error_kn_m']) for x in element_rows),
        original_general_joint_acceptance=None, production_scope_expanded=False)
    return raw, summary


def run(exe, output):
    output.mkdir(parents=True, exist_ok=False)
    cases = [('zero90', 90, [0., 0.], False, None),
             ('paired-sign', 0, [.05, -.05], False, (0, 2, 1, 5, 4, 3)),
             ('normal-only', 0, [.05, .05], False, (0, 2, 1, 5, 4, 3)),
             ('e-only', 0, [.05, -.05], False, None),
             ('step', 0, [.05, .04], False, None),
             ('one-zero', 90, [.05, 0.], False, None),
             ('split90', 90, [.05, .05], True, None)]
    cases += [(f'fold{angle}-{sign}', angle, [sign*.05, sign*.05], False, None)
              for angle in (15, 45, 90, 135, 179) for sign in (1, -1)]
    rotation = np.array([[0., -1., 0.], [0., 0., -1.], [1., 0., 0.]])
    rows = []; requests = []; layouts = []; assemblies = {}
    for name, angle, offsets, split, permutation in cases:
        variants = [('base', None, False)]
        if name in ('fold90-1', 'paired-sign', 'zero90'):
            variants += [('rotated', rotation, False), ('element-order', None, True)]
        for variant, rot, reverse in variants:
            positions, elements = panels(angle, offsets, split, permutation, rot)
            if reverse:
                elements.reverse()
            identity = name+'/'+variant
            for index, element in enumerate(elements):
                points = positions[element['nodes']]
                rigid = rigid_basis(points)
                fields = [dict(name='rigid'+str(i), d=rigid[:, i].tolist()) for i in range(6)]
                fields.append(dict(name='elastic-trial', d=(np.sin(np.arange(36)+.25)*1e-3).tolist()))
                requests.append(dict(id=identity+'/'+str(index), points=points.tolist(),
                    section=SECTION, e=element['e'], traction=TRACTION, fields=fields))
            components = [list(range(len(positions)))] if not split else [list(range(6)), list(range(6, 12))]
            layouts.append((identity, positions, elements, components))
    write(output/'requests.json', requests)
    started = time.monotonic()
    p = subprocess.Popen([str(exe)], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate(json.dumps(requests).encode())
    (output/'stderr.log').write_bytes(stderr); (output/'responses.json').write_bytes(stdout)
    write(output/'process.json', dict(command=[str(exe)], pid=p.pid, exit_code=p.returncode,
          elapsed_seconds=time.monotonic()-started, binary_sha256=hashlib.sha256(exe.read_bytes()).hexdigest()))
    if p.returncode:
        raise RuntimeError('probe failed; retain stderr.log')
    results = json.loads(stdout)
    repeat = subprocess.Popen([str(exe)], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    repeat_stdout, repeat_stderr = repeat.communicate(json.dumps(requests).encode())
    (output/'repeat-responses.json').write_bytes(repeat_stdout)
    (output/'repeat-stderr.log').write_bytes(repeat_stderr)
    reproducibility = dict(first_pid=p.pid, repeat_pid=repeat.pid, repeat_exit_code=repeat.returncode,
                           response_bytes_exact=stdout==repeat_stdout, repeated_kernel_requests=len(requests))
    for i, (identity, positions, elements, components) in enumerate(layouts):
        r, summary = global_report(positions, elements, results[2*i:2*i+2], components)
        independent = []
        for e in elements:
            kernel = IndependentElement(positions[e['nodes']], SECTION, e['e'])
            independent.append(dict(stiffness=kernel.stiffness.tolist(),
                neutral_stiffness=kernel.neutral.stiffness.tolist(),
                load=kernel.load(TRACTION).tolist(), fields=[]))
        ir, insummary = global_report(positions, elements, independent, components)
        k, ik = np.asarray(r['stiffness']), np.asarray(ir['stiffness'])
        summary['input_only_independent'] = dict(
            relative_stiffness_difference=norm(k-ik)/norm(ik),
            maximum_displacement_difference=float(np.max(np.abs(np.subtract(r['global_displacement'], ir['global_displacement'])))),
            global_observations=insummary)
        arms = [arm_observation(positions[e['nodes']], e['e'], .1) for e in elements]
        summary.update(id=identity, arms=arms,
                       neutral_point_gap_at_common_attachment_m=norm(np.subtract(arms[0]['arm_m'], arms[1]['arm_m'])))
        write(output/'assemblies'/(identity.replace('/', '-')+'.json'), r)
        write(output/'independent-assemblies'/(identity.replace('/', '-')+'.json'), ir)
        assemblies[identity] = r
        rows.append(summary)
    covariance = []
    for name in ('fold90-1', 'paired-sign', 'zero90'):
        base = np.asarray(assemblies[name+'/base']['stiffness'])
        rotated = np.asarray(assemblies[name+'/rotated']['stiffness'])
        mapping = np.zeros_like(base)
        for j in range(len(base)//3):
            mapping[3*j:3*j+3, 3*j:3*j+3] = rotation
        covariance.append(dict(case=name,
            rotated_stiffness_relative_difference=norm(rotated-mapping@base@mapping.T)/norm(base),
            element_order_stiffness_exact=bool(np.array_equal(base, assemblies[name+'/element-order']['stiffness'])),
            warning='rotated force seed is different; compliance is not an orientation-invariant comparison'))
    # Correct midside permutation, with signed e reversal for odd vertex order.
    order_requests = []; order_maps = []
    positions, elements = panels(90, [.05, .05])
    edge_mid = {frozenset((0, 1)):3, frozenset((1, 2)):4, frozenset((2, 0)):5}
    for side, e in enumerate(elements):
        original = positions[e['nodes']]
        for tri in itertools.permutations(range(3)):
            perm = list(tri)+[edge_mid[frozenset((tri[i], tri[j]))] for i, j in ((0, 1), (1, 2), (2, 0))]
            points = original[perm]
            signed_e = .05 if dot(frame(points), frame(original)) > 0 else -.05
            order_requests.append(dict(id=f'panel{side}-'+''.join(map(str, tri)), points=points.tolist(),
                section=SECTION, e=signed_e, traction=TRACTION, fields=[]))
            order_maps.append((side, perm))
    write(output/'node-order-requests.json', order_requests)
    ordered = subprocess.run([str(exe)], input=json.dumps(order_requests).encode(), capture_output=True, check=True)
    (output/'node-order-responses.json').write_bytes(ordered.stdout)
    order_rows = []
    # The first fold90-1/base corresponds to the same positions and local ordering.
    original_index = next(i for i, a in enumerate(layouts) if a[0] == 'fold90-1/base')
    for request, response, (side, perm) in zip(order_requests, json.loads(ordered.stdout), order_maps):
        source = np.asarray(results[2*original_index+side]['stiffness'])
        dofs = [6*i+a for i in perm for a in range(6)]
        reordered = np.zeros((36, 36)); reordered[np.ix_(dofs, dofs)] = response['stiffness']
        order_rows.append(dict(id=request['id'], permutation=perm, e_m=request['e'],
                              physical_id_stiffness_relative_difference=norm(source-reordered)/norm(source)))
    report = dict(passed=False, b03_complete=False, diagnostic_only=True,
                  production_scope_expanded=False, general_joint_acceptance=None,
                  public_kernel_unique_requests=len(requests)+len(order_requests),
                  public_kernel_total_evaluations=2*len(requests)+len(order_requests),
                  global_equilibrium_solves=2*len(rows), node_order=order_rows, reproducibility=reproducibility,
                  native_patch_solves=len(rows), input_only_independent_patch_solves=len(rows),
                  covariance=covariance, cases=rows,
                  independent_source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in [ROOT/'tools/plate6_offset_reference.py', ROOT/'tools/plate6_reference.py', ROOT/'tools/plate_reference.py']})
    write(output/'summary.json', report)
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--probe', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); print(json.dumps(run(a.probe.resolve(), a.output.resolve())))
