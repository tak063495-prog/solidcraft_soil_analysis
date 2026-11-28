"""Distorted, permitted flat patches under covariant loads; explicit use only."""
import argparse
import hashlib
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
from b03follow_math import assemble, cross, dot, mv, norm, rigid, write

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
_path = sys.path[:]
try:
    sys.path.insert(0, str(ROOT/'tools'))
    from plate6_offset_reference import Element as IndependentElement
finally:
    sys.path[:] = _path
SECTION = dict(id='b03follow', young_kpa=3e7, poisson=.2, thickness_m=.1,
               unit_weight_kn_m3=0., drilling_factor=1e-4)
RAW = ('membrane_strain', 'curvature_per_m', 'membrane_force_kn_m', 'bending_moment_kn',
       'transverse_shear_strain', 'transverse_shear_force_kn_m', 'drilling_mismatch_rad')
ENERGIES = ('membrane_energy_kn_m', 'bending_energy_kn_m',
            'transverse_shear_energy_kn_m', 'drilling_penalty_energy_kn_m')


def rotations():
    result = [('identity', np.eye(3)), ('cyclic', np.array([[0., -1., 0.], [0., 0., -1.], [1., 0., 0.]]))]
    for name, axis, degrees in [('axis123-37', [1., 2., 3.], 37.),
                                 ('axis123--83', [1., 2., 3.], -83.),
                                 ('axis2m1half-179', [2., -1., .5], 179.)]:
        a = np.array(axis)/norm(axis)
        angle = math.radians(degrees)
        c = cross(a)
        result.append((name, np.eye(3)+math.sin(angle)*c+(1-math.cos(angle))*(c@c)))
    return result


def patch(cells):
    positions = []
    for j in range(cells+1):
        for i in range(cells+1):
            x, y = i/cells, j/cells
            # Skewed quadrilateral, with deterministic interior distortion.
            dx = .13*math.sin(1.7*i+2.1*j) if 0 < i < cells and 0 < j < cells else 0.
            dy = .11*math.cos(2.3*i-.9*j) if 0 < i < cells and 0 < j < cells else 0.
            positions.append([1.3*x+.27*y+.12*x*y+dx, .8*y+.08*x*y+dy, 0.])
    triangles = []
    for j in range(cells):
        for i in range(cells):
            a = j*(cells+1)+i
            b, d = a+1, a+cells+1
            triangles.extend(((a, b, d+1), (a, d+1, d)))
    mids, elements = {}, []
    for tri in triangles:
        nodes = list(tri)
        for i, j in ((0, 1), (1, 2), (2, 0)):
            edge = tuple(sorted((tri[i], tri[j])))
            if edge not in mids:
                mids[edge] = len(positions)
                positions.append(((np.array(positions[edge[0]])+positions[edge[1]])/2).tolist())
            nodes.append(mids[edge])
        elements.append(nodes)
    return np.array(positions), elements


def polynomial_field(positions, e):
    x, y = np.asarray(positions)[:, :2].T
    epsx, epsy, gamma = 2e-5, -1e-5, 1.4e-5
    kx, ky, kxy = 3e-5, -2e-5, 1e-5
    theta = np.column_stack((-(kxy*x+ky*y), kx*x+kxy*y, np.zeros(len(x))))
    u = np.column_stack((epsx*x+gamma*y/2, epsy*y+gamma*x/2,
                         -.5*(kx*x*x+2*kxy*x*y+ky*y*y)))
    u -= theta@cross([0., 0., e]).T
    return np.column_stack((u, theta)).ravel()


def query(probe, output, label, requests):
    write(output/(label+'-requests.json'), requests)
    started = time.monotonic()
    p = subprocess.Popen([str(probe)], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate(json.dumps(requests, allow_nan=False).encode())
    (output/(label+'-responses.json')).write_bytes(stdout)
    (output/(label+'-stderr.log')).write_bytes(stderr)
    write(output/(label+'-process.json'), dict(pid=p.pid, exit_code=p.returncode,
          elapsed_seconds=time.monotonic()-started, requests=len(requests), binary_sha256=hashlib.sha256(probe.read_bytes()).hexdigest()))
    if p.returncode:
        raise RuntimeError('probe failed; retained '+label)
    rows = json.loads(stdout)
    if len(rows) != len(requests) or any(r.get('error') for r in rows) or any(a['id'] != b['id'] for a, b in zip(rows, requests)):
        raise ValueError('kernel responses differ from complete accepted request inventory')
    return rows


def solve(k, positions, force):
    r = rigid(positions)
    q, _ = np.linalg.qr(r, mode='complete')
    v = q[:, 6:]
    reduced = v.T@k@v
    d = mv(v, np.linalg.solve(reduced, mv(v.T, force)))
    internal = mv(k, d)
    residual = internal-force
    eigenvalues = np.linalg.eigvalsh(k)
    cutoff = 1e-10*float(eigenvalues[-1])
    zeros = int(np.sum(eigenvalues <= cutoff))
    report = dict(eigenvalues=eigenvalues, inherited_cutoff=cutoff, algebraic_zero_modes=zeros,
        additional_zero_modes=zeros-6,
        rigid_normalized_residuals=[norm(mv(k, r[:, i]))/(norm(k)*norm(r[:, i])) for i in range(6)],
        force_residual_l2_kn=norm(residual.reshape(-1, 6)[:, :3]),
        moment_residual_l2_kn_m=norm(residual.reshape(-1, 6)[:, 3:]),
        force_rigid_orthogonality=[dot(r[:, i], force) for i in range(6)],
        external_work_kn_m=dot(force, d), twice_energy_kn_m=dot(d, internal),
        work_balance_kn_m=dot(force, d)-dot(d, internal),
        reduced_condition_number=float(np.linalg.cond(reduced)),
        displacement=d, internal_force=internal, residual=residual, rigid_basis=r)
    return report


def stresses(recovery):
    axes = np.asarray(recovery['local_axes'])
    tensors = []
    for sample in recovery['samples']:
        n, m = np.asarray(sample['membrane_force_kn_m']), np.asarray(sample['bending_moment_kn'])
        for sign in (-1, 1):
            sx, sy, sxy = n/.1+sign*6*m/.1**2
            local = np.array([[sx, sxy, 0.], [sxy, sy, 0.], [0., 0., 0.]])
            tensors.append(axes.T@local@axes)
    return np.array(tensors)


def fields_difference(base, other, rotation):
    raw = {key: [] for key in RAW}
    stress = []
    attach = []
    energy = []
    for a, b in zip(base, other, strict=True):
        for sa, sb in zip(a['samples'], b['samples'], strict=True):
            for key in RAW:
                raw[key].append(np.asarray(sb[key])-sa[key])
            attach.append((np.asarray(sb['bending_moment_kn'])-b['offset_m']*np.asarray(sb['membrane_force_kn_m']))-
                          (np.asarray(sa['bending_moment_kn'])-a['offset_m']*np.asarray(sa['membrane_force_kn_m'])))
        stress.append(np.array([rotation.T@v@rotation for v in stresses(b)])-stresses(a))
        energy.extend(float(b[key])-float(a[key]) for key in ENERGIES)
    return dict(raw_max_abs={k: max(abs(float(x)) for row in values for x in np.asarray(row).ravel()) for k, values in raw.items()},
                attachment_moment_max_abs_kn=float(np.max(np.abs(attach))),
                derived_top_bottom_inplane_stress_max_abs_kpa=float(np.max(np.abs(stress))),
                component_energy_max_abs_kn_m=max(abs(x) for x in energy))


def run(probe, output):
    output.mkdir(parents=True, exist_ok=False)
    layouts, requests = [], []
    for cells in (1, 2, 3):
        base_positions, elements = patch(cells)
        size = 6*len(base_positions)
        q, _ = np.linalg.qr(rigid(base_positions), mode='complete')
        v = q[:, 6:]
        base_force = mv(v, mv(v.T, np.sin(np.arange(size)+.125)))
        for e in (0., .025, -.05):
            for name, rotation in rotations():
                identity = f'{2*cells*cells}-e{e}-{name}'
                positions = base_positions@rotation.T
                mapping = np.kron(np.eye(size//3), rotation)
                indices = []
                for ei, nodes in enumerate(elements):
                    indices.append(len(requests))
                    requests.append(dict(id=identity+'/'+str(ei), points=positions[nodes].tolist(),
                        section=SECTION, e=e, traction=(rotation@np.array([1., .5, -.25])).tolist(), fields=[]))
                layouts.append(dict(id=identity, cells=cells, rotation=rotation, rotation_name=name, e=e,
                    positions=positions, elements=elements, request_indices=indices, mapping=mapping,
                    force=mv(mapping, base_force), polynomial=mv(mapping, polynomial_field(base_positions, e))))
    results = query(probe, output, 'matrices', requests)
    recovery_requests, records, independent_recoveries = [], [], []
    for layout in layouts:
        native_parts, independent_parts, kernels = [], [], []
        for nodes, index in zip(layout['elements'], layout['request_indices'], strict=True):
            dofs = [6*n+a for n in nodes for a in range(6)]
            kernel = IndependentElement(layout['positions'][nodes], SECTION, layout['e'])
            native_parts.append((dofs, np.asarray(results[index]['stiffness'])))
            independent_parts.append((dofs, kernel.stiffness))
            kernels.append(kernel)
        size = len(layout['force'])
        k, ik = assemble(native_parts, size), assemble(independent_parts, size)
        n = solve(k, layout['positions'], layout['force'])
        ind = solve(ik, layout['positions'], layout['force'])
        fields = {'solved': [], 'polynomial': []}
        for nodes, index, kernel in zip(layout['elements'], layout['request_indices'], kernels, strict=True):
            dofs = [6*v+a for v in nodes for a in range(6)]
            request = dict(requests[index])
            request['fields'] = [dict(name='solved', d=n['displacement'][dofs].tolist()),
                                 dict(name='polynomial', d=layout['polynomial'][dofs].tolist())]
            recovery_requests.append(request)
            for key, d in (('solved', ind['displacement']), ('polynomial', layout['polynomial'])):
                rec = kernel.recover(d[dofs])
                rec['offset_m'] = layout['e']
                fields[key].append(rec)
        independent_recoveries.append(fields)
        record = {k: v for k, v in layout.items() if k != 'mapping'}
        record.update(native=n, independent=ind, native_stiffness=k, independent_stiffness=ik,
            relative_stiffness_difference=norm(k-ik)/norm(ik),
            native_independent_translation_max_abs_m=float(np.max(np.abs((n['displacement']-ind['displacement']).reshape(-1, 6)[:, :3]))),
            native_independent_rotation_max_abs_rad=float(np.max(np.abs((n['displacement']-ind['displacement']).reshape(-1, 6)[:, 3:]))))
        records.append(record)
    recovered = query(probe, output, 'recoveries', recovery_requests)
    offset = 0
    baselines = {}
    rows = []
    for record, layout, ind_fields in zip(records, layouts, independent_recoveries, strict=True):
        native_fields = {'solved': [], 'polynomial': []}
        for response in recovered[offset:offset+len(record['elements'])]:
            for field in response['fields']:
                rec = field['recovery']
                rec['offset_m'] = layout['e']
                native_fields[field['name']].append(rec)
        offset += len(record['elements'])
        record['native_fields'], record['independent_fields'] = native_fields, ind_fields
        record['field_independence'] = {key: fields_difference(native_fields[key], ind_fields[key], np.eye(3)) for key in native_fields}
        base_key = (layout['cells'], layout['e'])
        if layout['rotation_name'] == 'identity':
            baselines[base_key] = record
        base = baselines[base_key]
        covariance = {}
        for side, fields_key in (('native', 'native_fields'), ('independent', 'independent_fields')):
            d = record[side]['displacement'].reshape(-1, 3)@layout['rotation']
            db = np.asarray(base[side]['displacement']).reshape(-1, 3)
            delta = (d-db).reshape(-1, 6)
            kb = base[side+'_stiffness']
            mapping = layout['mapping']
            covariance[side] = dict(
                translation_max_abs_m=float(np.max(np.abs(delta[:, :3]))),
                rotation_max_abs_rad=float(np.max(np.abs(delta[:, 3:]))),
                relative_stiffness_difference=norm(record[side+'_stiffness']-mapping@kb@mapping.T)/norm(kb),
                external_work_difference_kn_m=record[side]['external_work_kn_m']-base[side]['external_work_kn_m'],
                twice_energy_difference_kn_m=record[side]['twice_energy_kn_m']-base[side]['twice_energy_kn_m'],
                fields={key: fields_difference(base[fields_key][key], record[fields_key][key], layout['rotation']) for key in native_fields})
        record['rotation_covariance'] = covariance
        write(output/'assemblies'/(record['id']+'.json'), record)
        rows.append(dict(id=record['id'], tri6=len(record['elements']), nodes=len(record['positions']),
            dofs=len(record['force']), e_m=record['e'], rotation=record['rotation_name'],
            native={k: v for k, v in record['native'].items() if k not in ('eigenvalues','displacement','internal_force','residual','rigid_basis')},
            independent={k: v for k, v in record['independent'].items() if k not in ('eigenvalues','displacement','internal_force','residual','rigid_basis')},
            relative_stiffness_difference=record['relative_stiffness_difference'],
            native_independent_translation_max_abs_m=record['native_independent_translation_max_abs_m'],
            native_independent_rotation_max_abs_rad=record['native_independent_rotation_max_abs_rad'],
            field_independence=record['field_independence'], rotation_covariance=covariance))
    report = dict(passed=False, b03_complete=False, diagnostic_only=True, production_scope_expanded=False,
        native_solves=len(rows), independent_solves=len(rows), kernel_evaluations=2*len(requests),
        maximum_tri6=max(r['tri6'] for r in rows), maximum_dofs=max(r['dofs'] for r in rows),
        physical_reference_solution_certified=False, shared_global_assembly_and_gauge_algorithm=True,
        stress_definition='N/t +/- 6*Mp/t^2, local in-plane tensor transformed to global; not soil stress',
        cases=rows)
    write(output/'summary.json', report)
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--probe', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    run(a.probe.resolve(), a.output.resolve())
