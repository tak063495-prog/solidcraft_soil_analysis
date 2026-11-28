"""Explicit, bounded internal observations; production geometry guards remain active."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from b01cause_geometry import EDGES, PARAM, dot, geometry, matvec, mesh, metric, norm, rigid, surface

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SECTION = dict(id='b01cause', young_kpa=3e7, poisson=.2, thickness_m=.1,
               unit_weight_kn_m3=0., drilling_factor=1e-4)


def request(name, points, fields=True, section=None):
    x = np.array(points, float)
    probes = []
    if fields:
        probes += [dict(name='rigid-'+str(i), d=rigid(x)[:, i].tolist()) for i in range(6)]
        d = np.zeros((6, 6)); d[:, :3] = .001*x
        probes.append(dict(name='isotropic-membrane', d=d.ravel().tolist()))
    return dict(id=name, points=x.tolist(), section=section or SECTION, fields=probes)


def cases():
    plane = surface('plane', PARAM)
    bow = plane.copy(); bow[3, 1] += .02
    bulge = plane.copy(); bulge[3:, 2] += .02
    rows = [request('plane', plane), request('inplane-bow', bow), request('outplane-bulge', bulge)]
    rows += [request(k, surface(k, PARAM)) for k in ('cylinder', 'sphere')]
    # Every vertex ordering carries its actual associated edge node.
    for kind, x in [('plane', plane), ('outplane-bulge', bulge)]:
        for p in itertools.permutations(range(3)):
            index = list(p)
            for a, b in EDGES:
                index.append(next(3+e for e, (i, j) in enumerate(EDGES) if {i, j} == {p[a], p[b]}))
            r = request(kind+'-order-'+''.join(map(str, p)), x[index], fields=False)
            r['original_node_indices'] = index; rows.append(r)
    for scale in (1e-3, 1., 1e3):
        for shift in (0., 1e6):
            x = plane*scale+[shift, 0., 0.]
            cutoff = 64*np.finfo(float).eps*max(norm(x[i]-x[j]) for i, j in EDGES)
            for direction in ('inplane', 'outplane'):
                for factor in (.5, 1., 2.):
                    y = x.copy(); y[3, 1 if direction == 'inplane' else 2] += factor*cutoff
                    r = request(f'boundary-{scale}-{shift}-{direction}-{factor}', y, fields=False)
                    r.update(requested_cutoff_multiple=factor, direction=direction); rows.append(r)
    for digits in (6, 12, 15):
        rows.append(request('sphere-round-'+str(digits), np.round(surface('sphere', PARAM), digits), fields=False))
    x = surface('sphere', PARAM); x[3, 2] = np.nextafter(x[3, 2], np.inf)
    rows.append(request('sphere-nextafter', x, fields=False))
    for thickness in (.01, .001):
        for drill in (1e-8, 1e-4, .01):
            s = dict(SECTION, thickness_m=thickness, drilling_factor=drill)
            rows.append(request(f'bulge-t{thickness}-d{drill}', bulge, section=s))
    meshes = []
    for kind in ('plane', 'cylinder', 'sphere'):
        for n in (1, 2):
            x, elements, boundary = mesh(kind, n)
            ids = []
            for e, nodes in enumerate(elements):
                name = f'mesh-{kind}-{n}-{e}'; ids.append(name)
                rows.append(request(name, x[nodes], fields=False))
                facet = x[nodes].copy()
                for j, (a, b) in enumerate(EDGES):
                    facet[3+j] = (facet[a]+facet[b])*.5
                rows.append(request(name+'-facet', facet, fields=False))
            meshes.append(dict(kind=kind, n=n, points=x.tolist(), elements=elements, boundary=boundary, ids=ids))
    return rows, meshes


def spectrum(k, x):
    r = rigid(x); values = np.linalg.eigvalsh((k+k.T)*.5)
    threshold = 1e-10*float(values[-1])
    residual = [norm(matvec(k, r[:, i]))/(norm(k)*norm(r[:, i])) for i in range(6)]
    energy = [dot(r[:, i], matvec(k, r[:, i]))*.5 for i in range(6)]
    return dict(eigenvalues=values.tolist(), inherited_rank_threshold=threshold,
                modes_below_inherited_threshold=int(np.sum(values <= threshold)),
                negative_below_minus_threshold=int(np.sum(values < -threshold)),
                physical_six_rigid_normalized_residual=residual,
                physical_six_rigid_energy_kn_m=energy,
                classification='algebraic mode counts and actual rigid residual must be read together')


def assembled(meta, by_id):
    x = np.array(meta['points']); size = 6*len(x); k = np.zeros((size, size))
    facets = []
    for name, nodes in zip(meta['ids'], meta['elements']):
        dofs = [6*i+a for i in nodes for a in range(6)]
        ke = np.array(by_id[name]['stiffness']); k[np.ix_(dofs, dofs)] += ke
        facets.append(np.array_equal(ke, np.array(by_id[name+'-facet']['stiffness'])))
    result = dict(kind=meta['kind'], n=meta['n'], nodes=len(x), elements=len(meta['elements']),
                  dofs=size, spectrum=spectrum(k, x),
                  element_stiffness_exactly_equals_same_corner_faceted_control=facets)
    if meta['n'] == 2:
        boundary = [6*i+a for i in meta['boundary'] for a in range(6)]
        free = sorted(set(range(size))-set(boundary))
        solves = []
        # Fixed boundary prescriptions, no external nodal loads. Rigid motion is
        # a compatible exact zero-stress field independently of shell curvature.
        for axis in range(3, 6):
            target = rigid(x)[:, axis]*1e-3
            rhs = -k[np.ix_(free, boundary)]@target[boundary]
            u = target.copy(); t = time.monotonic()
            u[free] = np.linalg.solve(k[np.ix_(free, free)], rhs)
            force = matvec(k, u)
            solves.append(dict(name='rigid-boundary-'+str(axis-3),
                               elapsed_seconds=time.monotonic()-t,
                               interior_translation_force_l2_kn=norm(force[[i for i in free if i % 6 < 3]]),
                               interior_rotation_moment_l2_kn_m=norm(force[[i for i in free if i % 6 >= 3]]),
                               mixed_dof_algebraic_relative_residual=norm(force[free])/max(norm(force[boundary]), 1e-12),
                               displacement_error=metric(u.reshape(-1, 6)[:, :3], target.reshape(-1, 6)[:, :3]),
                               rotation_error=metric(u.reshape(-1, 6)[:, 3:], target.reshape(-1, 6)[:, 3:]),
                               boundary_translation_force_l2_kn=norm(force[[i for i in boundary if i % 6 < 3]]),
                               boundary_rotation_moment_l2_kn_m=norm(force[[i for i in boundary if i % 6 >= 3]]),
                               energy_kn_m=.5*dot(u, force),
                               solution=u.tolist(), prescribed_rigid=target.tolist(), residual=force.tolist()))
        result.update(boundary_dofs=boundary, free_dofs=free, global_equilibrium_solves=solves)
    return result


def main(output, probe):
    output.mkdir(parents=True, exist_ok=True)
    requests, meshes = cases()
    raw = json.dumps(requests, allow_nan=False).encode()
    (output/'requests.json').write_bytes(raw)
    started = time.monotonic()
    p = subprocess.Popen([str(probe)], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate(raw)
    (output/'probe.stdout.json').write_bytes(stdout); (output/'probe.stderr.log').write_bytes(stderr)
    if p.returncode:
        raise subprocess.CalledProcessError(p.returncode, [str(probe)])
    results = json.loads(stdout)
    by_id = {r['id']: r for r in results}
    if len(by_id) != len(requests):
        raise ValueError('probe inventory mismatch')
    observations = []
    for req in requests:
        r = by_id[req['id']]
        value = dict(id=req['id'], production_guard=r['production_guard'], geometry=geometry(req['points']))
        if 'stiffness' in r:
            value['spectrum'] = spectrum(np.array(r['stiffness']), req['points'])
            value['production_stiffness_exact_when_accepted'] = r['production_stiffness_exact_when_accepted']
            value['fields'] = r['fields']
        if 'original_node_indices' in req:
            dofs = [6*i+a for i in req['original_node_indices'] for a in range(6)]
            back = np.empty((36, 36)); back[np.ix_(dofs, dofs)] = r['stiffness']
            original = by_id[req['id'].split('-order-')[0]]['stiffness']
            value['stiffness_after_physical_id_permutation'] = metric(back, original)
        observations.append(value)
    report = dict(passed=False, b01_complete=False, diagnostic_only=True,
                  production_scope_expanded=False, curved_element_acceptance=None,
                  process=dict(pid=p.pid, exit_code=p.returncode, elapsed_seconds=time.monotonic()-started,
                               request_sha256=hashlib.sha256(raw).hexdigest(),
                               probe_sha256=hashlib.sha256(probe.read_bytes()).hexdigest()),
                  observations=observations, meshes=[assembled(m, by_id) for m in meshes])
    (output/'internal-observations.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(requests=len(requests), assembled_meshes=len(meshes),
                          global_solves=sum(len(m.get('global_equilibrium_solves', [])) for m in report['meshes']))))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--probe', type=Path, required=True)
    a = p.parse_args(); main(a.output.resolve(), a.probe.resolve())
