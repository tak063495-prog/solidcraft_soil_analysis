"""Record projected/compatible raw MCC field differences on the same mesh.

Neither formulation generates the other's solution. Raw GP stress/history is
compared without projecting/replacing it. Plastic volume is a quadrature-based
cumulative-history indicator, not an exact plastic-zone boundary measurement.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import fem_reference as ref
from assess_projected_volume_diagnostic import fields, read, sha
import assess_projected_volume_diagnostic as field_module
import mcc_projected_cycle_diagnostic as cycle


def bindings(report, key):
    paths = report[key]
    for name, digest in paths.items():
        assert sha(Path(name)) == digest, name
    return {Path(name).name: Path(name) for name in paths}


def integration_rules(mesh):
    rules = {}
    for order in (5, 7):
        shapes, weights = [], []
        for e in range(len(mesh['tetrahedra'])):
            q = ref.volume(mesh, e, independent=True, independent_order=order)
            shapes.append(q['shape'])
            weights.append(q['weight'])
        rules[order] = np.stack(shapes), np.stack(weights)
    return rules


def displacement_integrals(mesh, rules, a, b, origin_a, origin_b):
    active = np.asarray(a['active_elements'])
    cells = ref.cells(mesh)[active]
    ua, ub = np.asarray(a['displacement_m']), np.asarray(b['displacement_m'])
    results = {}
    for order, (shapes, weights) in rules.items():
        shape, weight = shapes[active], weights[active]
        volume = float(weight.sum())
        values = {}
        for name, delta, reference in (
                ('from_seed', ua-ub, ub),
                ('from_construction_origin', (ua-origin_a)-(ub-origin_b), ub-origin_b)):
            difference = np.einsum('eqn,ena->eqa', shape, delta[cells])
            ref_value = np.einsum('eqn,ena->eqa', shape, reference[cells])
            error_squared = float(np.einsum('eq,eqa,eqa->', weight, difference, difference))
            reference_squared = float(np.einsum('eq,eqa,eqa->', weight, ref_value, ref_value))
            values[name] = dict(difference_squared_integral_m5=error_squared,
                                reference_squared_integral_m5=reference_squared,
                                rms_difference_m=float(np.sqrt(error_squared/volume)),
                                relative_l2=float(np.sqrt(error_squared/reference_squared)) if reference_squared > 0 else None)
        results[order] = dict(active_volume_m3=volume, fields=values)
    changes = {}
    for name in results[5]['fields']:
        left, right = results[5]['fields'][name], results[7]['fields'][name]
        changes[name] = abs(left['rms_difference_m']-right['rms_difference_m'])
        assert changes[name] < 1e-12, 'quadratic displacement integral failed quadrature refinement'
        for key in ('difference_squared_integral_m5', 'reference_squared_integral_m5'):
            assert abs(left[key]-right[key])/max(abs(right[key]), 1e-40) < 1e-10
    return dict(quadrature_orders=[5, 7], order5=results[5], order7=results[7],
                rms_change_under_quadrature_refinement_m=changes,
                definition='Physical volume integral of quadratic nodal displacement differences; raw GP constitutive fields are not interpolated')


def compare(mesh_path, projected_path, compatible_path, output):
    mesh, projected, compatible = map(read, (mesh_path, projected_path, compatible_path))
    assert projected['schema'] == 'soil-fem-projected-mcc-cycle-diagnostic/1'
    assert compatible['schema'] == 'soil-fem-independent-mcc-rate-fem/1'
    p_files = bindings(projected, 'input_and_source_sha256')
    c_files = bindings(compatible, 'input_sha256')
    assert sha(mesh_path) == sha(p_files['mesh.json']) == sha(c_files['mesh.json'])
    assert sha(p_files['ground.json']) == sha(c_files['ground.json'])
    assert sha(p_files['mesh-identity.json']) == sha(c_files['mesh-identity.json'])
    assert cycle.physical_case(read(p_files['case.json'])) == cycle.physical_case(read(c_files['case.json']))
    weights = np.stack([ref.volume(mesh, e)['weight'] for e in range(len(mesh['tetrahedra']))])
    assert np.all(weights > 0)
    rules = integration_rules(mesh)
    origin_a, origin_b = np.asarray(projected['reference_displacement_m']), np.asarray(compatible['reference_displacement_m'])
    rows = []
    assert len(projected['frames']) == len(compatible['frames'])
    for a, b in zip(projected['frames'], compatible['frames']):
        assert a['name'] == b['name'] and a['active_elements'] == b['active_elements']
        sa, sb = np.asarray(a['material_state']), np.asarray(b['material_state'])
        assert sa.shape == sb.shape == (*weights.shape, 14)
        active = np.asarray(a['active_elements'])
        w = weights[active]
        summaries = {}
        for label, state in (('projected', sa), ('compatible', sb)):
            local = state[active]
            plastic = np.linalg.norm(local[..., 8:], axis=-1) > 1e-10
            summaries[label] = dict(active_volume_m3=float(w.sum()),
                active_points=int(w.size), active_plastic_history_points=int(np.count_nonzero(plastic)),
                plastic_history_quadrature_volume_m3=float(w[plastic].sum()),
                volume_mean_stress_kpa=np.einsum('eq,eqa->a', w, local[..., :6]).tolist(),
                volume_mean_pc_kpa=float((w*np.exp(local[..., 6])).sum()/w.sum()),
                volume_mean_void_ratio=float((w*(np.exp(local[..., 7])-1)).sum()/w.sum()))
            summaries[label]['volume_mean_stress_kpa'] = [v/float(w.sum()) for v in summaries[label]['volume_mean_stress_kpa']]
        rows.append(dict(name=a['name'], raw_full_field_difference=fields(a, b),
                         material_summaries=summaries,
                         physical_volume_displacement_difference=displacement_integrals(mesh, rules, a, b, origin_a, origin_b),
                         projected_force_residual=a['relative_work_conjugate_force_residual'],
                         compatible_force_residual=b['relative_free_force_residual']))
    paths = [mesh_path, projected_path, compatible_path, Path(__file__), Path(ref.__file__), Path(field_module.__file__), Path(cycle.__file__)]
    report = dict(schema='soil-fem-mcc-discretization-difference/1', assessment_recorded=True,
                  elements=len(mesh['tetrahedra']), projected_reference_steps=projected['reference_steps'],
                  compatible_reference_steps=compatible['reference_steps'], frames=rows,
                  goal_complete=False, p6_complete=False, spatial_convergence_verified=False,
                  sha256={str(p): sha(p) for p in paths},
                  limitations='Different spatial formulations on one mesh; raw GP differences are not implementation accuracy failures or continuum error bounds; independent integration/spatial refinement remains required')
    with output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps([dict(name=r['name'], raw_displacement_relative=r['raw_full_field_difference']['displacement_m']['relative_l2'],
                          volume_integrated_displacement_relative=r['physical_volume_displacement_difference']['order7']['fields']['from_seed']['relative_l2'],
                          maximum_plastic_strain_difference=r['raw_full_field_difference']['plastic_strain']['maximum_absolute'],
                          plastic_history_volume_m3={k: v['plastic_history_quadrature_volume_m3'] for k, v in r['material_summaries'].items()}) for r in rows]), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ('mesh', 'projected', 'compatible', 'output'):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    compare(*(getattr(args, k).resolve() for k in ('mesh', 'projected', 'compatible', 'output')))
