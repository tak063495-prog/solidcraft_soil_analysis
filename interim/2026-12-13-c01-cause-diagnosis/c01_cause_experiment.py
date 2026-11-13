"""Explicit verification-only contact quadrature/reference refinement experiment.

Does not edit production System, Water, Case, controls, or the retained studies.
"""
import argparse
import copy
from math import fsum
from pathlib import Path
import time
import c01_cause_diagnostics as diag


def contact_system(mesh, configuration, order):
    import numpy as np
    from joint_reference import System
    import fem_reference as fem
    if not 3 <= order <= 8:
        raise ValueError('contact order must be 3..8')
    system = System(mesh, configuration)
    if order == 4:
        return system  # exactly the original reference geometry/point ordering
    faces = {}
    for point in system.points:
        faces.setdefault((point['joint_id'], point['face_index']), point)
    xyz = np.asarray(mesh['nodes'])
    points = []
    for p in faces.values():
        na, nb = p['nodes_a'], p['nodes_b']
        if len(na) != 6:
            raise ValueError('experiment is restricted to quadratic matched faces')
        pa = xyz[na]
        for u, wu in zip(*fem.jacobi(order, 1)):
            for v, wv in zip(*fem.jacobi(order, 0)):
                l = np.array([1-u-(1-u)*v, u, (1-u)*v])
                dl = np.array([[-1., -1.], [1., 0.], [0., 1.]])
                shape = np.r_[l*(2*l-1), [4*l[a]*l[b] for a,b in ((0,1),(1,2),(2,0))]]
                gradient = np.vstack(((4*l-1)[:,None]*dl, [4*(l[a]*dl[b]+l[b]*dl[a]) for a,b in ((0,1),(1,2),(2,0))]))
                jac = (pa-pa[0]).T@gradient
                area = np.cross(jac[:,0], jac[:,1])*wu*wv
                weight = np.linalg.norm(area)
                q = copy.deepcopy(p)
                q.update(position_m=(pa[0]+shape@(pa-pa[0])).tolist(), normal=(area/weight).tolist(), weight_m2=float(weight), shape=shape.tolist())
                points.append(q)
    system.points = points
    system.normal = np.asarray([p['normal'] for p in points])
    system.weight = np.asarray([p['weight_m2'] for p in points])
    system.kn = np.asarray([p['material']['normal_stiffness_kpa_m'] for p in points])
    system.kt = np.asarray([p['material']['tangential_stiffness_kpa_m'] for p in points])
    system.mu = np.asarray([p['material']['friction_coefficient'] for p in points])
    system.c = np.zeros((len(points), 3, xyz.size))
    for i, p in enumerate(points):
        for a, b, w in zip(p['nodes_a'], p['nodes_b'], p['shape']):
            for axis in range(3):
                system.c[i, axis, 3*a+axis] -= w
                system.c[i, axis, 3*b+axis] += w
    return system


def solve(output, steps, order, control):
    import numpy as np
    old, reference = diag.dependencies()
    from joint_water_reference import Water
    mesh = old.read(diag.ROOT/'examples/joint-water-column-tet10-mesh.json')
    ground = old.fem.normalize_ground(old.read(diag.ROOT/'examples/ground-joint-water-mc.json'))
    case = old.fem.normalize_stage_case(old.read(diag.ROOT/'examples/stages-joint-water-cycle.json'))
    case = old.configuration(case, steps, 5e-5)
    if control == 'switch':
        case['joint_system']['global_step_control'] = dict(switch_bisections=4)
        case['settings']['minimum_increment'] = min(case['settings']['minimum_increment'], 1/(steps*32))
    if control == 'path':
        case['joint_system']['global_step_control'] = dict(path_error=dict(slip_scale_m=1e-4, dissipation_scale_kpa_m=.01))
        case['settings']['minimum_increment'] = min(case['settings']['minimum_increment'], 1/(steps*65536))
    output.mkdir(parents=True, exist_ok=False)
    for name, value in (('mesh.json', mesh), ('ground.json', ground), ('case.json', case)):
        diag.write(output/name, value)

    class DiagnosticWater(Water):
        def __init__(self, mesh, ground, case, layout=None):
            super().__init__(mesh, ground, case, layout)
            self.contact = contact_system(mesh, case['joint_system'], order)
            self.faces = list({(p['joint_id'], p['face_index']):p for p in self.contact.points}.values())

    reference.Water = DiagnosticWater  # only this newly loaded diagnostic module
    system = contact_system(mesh, case['joint_system'], order)
    writer = old.Writer(output/'independent')
    started = time.monotonic()
    sampled = []

    def sink(row, fields):
        writer(row, fields)
        # Fixed-order area integrals: compare different GP rules without inventing
        # a correspondence/interpolation of path-dependent GP states.
        integral = {key:fsum(float(x)*float(w) for x,w in zip(fields[key],system.weight)) for key in ('plastic_dissipation','elastic_energy','endpoint_work','numerical_dissipation')}
        integral['plastic_slip_integral'] = [fsum(float(x)*float(w) for x,w in zip(fields['plastic_slip'][:,i],system.weight)) for i in range(3)]
        integral['traction_integral'] = [fsum(float(x)*float(w) for x,w in zip(fields['traction'][:,i],system.weight)) for i in range(3)]
        sampled.append(dict(stage=row['stage'], progress=row['progress'], integrals=integral))

    kwargs = dict(maximum_m=5e-5)
    if control == 'switch':
        kwargs['switch_bisections'] = 4
    if control == 'path':
        kwargs['path_error'] = dict(slip_scale_m=1e-4, dissipation_scale_kpa_m=.01)
    try:
        value = reference._reference(mesh, ground, case, steps, 1, 7, False, True, sink, **kwargs)
        old.write(output/'result.json', value)
    except Exception as error:
        diag.write(output/'failure.json', dict(passed=False, complete=False, error_type=type(error).__name__, error=str(error), elapsed_seconds=time.monotonic()-started))
        raise
    finally:
        writer.close()
    residual = value['maximum_relative_force_residual']
    if residual > 1e-10 or len(value['frames']) != 7:
        raise ValueError('original reference force gate or seven-phase completion failed')
    diag.write(output/'integrals.json', dict(samples=sampled))
    diag.write(output/'geometry.json', diag.geometry(system, mesh))
    diag.write(output/'summary.json', dict(schema='c01-cause-experiment/1', starting_commit=diag.BASE, passed=False, c01_decision='未達', complete=True, contact_order=order, contact_gp_count=len(system.points), global_steps=steps, control=control, accepted_states=len(value['path']), maximum_force_residual=residual, original_independent_force_gate=1e-10, elapsed_seconds=time.monotonic()-started, uses_native_results=False, production_changed=False))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--steps', type=int, required=True)
    p.add_argument('--contact-order', type=int, choices=range(3,9), default=4)
    p.add_argument('--control', choices=('none','switch','path'), default='none')
    p.add_argument('--pc-series', action='store_true')
    args = p.parse_args()
    if args.steps not in (16,32,64,128,256,512) and not (args.pc_series and args.steps in (1024,2048)):
        p.error('large series needs explicit --pc-series on user PC')
    solve(args.output, args.steps, args.contact_order, args.control)


if __name__ == '__main__':
    main()
