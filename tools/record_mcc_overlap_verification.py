"""Seal analytic overlap checks and the comparable planar MCC mesh studies."""
import argparse
import json
import math
from pathlib import Path
from tet10_overlap import np
import fem_reference as ref
import freeze_sources
from record_mcc_spatial_step import archive_check
from mcc_projected_volume_diagnostic import read, sha


def benchmark(mesh, ground, case, layout):
    assert len(ground['materials']) == len(ground['bodies']) == 1
    material = ground['materials'][0]
    law = material['mcc']
    assert material['unit_weight_kn_m3'] == 20
    initial = ground['initialization']
    assert initial['method'] == 'equilibrate' and initial['water_level_m'] is None
    assert initial['surface_surcharge_kpa'] == 0 and initial['geostatic']['initial_boundary_factor'] == 0
    assert initial['geostatic']['seed']['kind'] == 'integration_points'
    assert case['settings']['kinematics'] == 'side_rollers'
    assert not case['boundary']['constraints'] and not case['boundary'].get('nodal_loads')
    assert case['boundary']['surface_loads'] == [dict(kind='pressure', face_set='free_surface', pressure_kpa=5.)]
    assert case['pressure_on_original_surface_only']
    stages = case['stages']
    assert len(stages) == 3
    assert stages[0]['kind'] == stages[1]['kind'] == 'boundary'
    assert stages[0]['target_factor'] == 1 and stages[1]['target_factor'] == 0
    assert stages[2]['kind'] == 'excavate' and stages[2]['bodies'] == ['Pit'] and stages[2]['remaining_surface_pressure_kpa'] == 0
    points = initial['geostatic']['seed']['points']
    assert len(points) == len(layout)
    maximum_position_error, maximum_stress_error, maximum_pc_error, maximum_e_error = 0., 0., 0., 0.
    minimum_p = float('inf')
    for point, actual in zip(points, layout):
        assert (point['element_index'], point['local_index']) == (actual['element_index'], actual['local_index'])
        assert point['material_id'] == mesh['material_ids'][point['element_index']]
        maximum_position_error = max(maximum_position_error, float(abs(np.asarray(point['position_m'])-actual['position_m']).max()))
        x, y, z = point['position_m']
        p = 20*(10-x/3-z)
        assert p > 0
        minimum_p = min(minimum_p, p)
        state = point['state']
        maximum_stress_error = max(maximum_stress_error, float(abs(np.asarray(state['effective_stress_compression_positive_kpa'])-[p, p, p, 0, 0, 0]).max()))
        maximum_pc_error = max(maximum_pc_error, abs(state['preconsolidation_pressure_kpa']-4*p))
        e = law['specific_volume_ncl_at_reference']-law['lambda']*math.log(4*p/law['reference_pressure_kpa'])+law['kappa']*math.log(4)-1
        maximum_e_error = max(maximum_e_error, abs(state['void_ratio']-e))
        assert not any(state['plastic_strain_compression_positive'])
    assert max(maximum_position_error, maximum_stress_error, maximum_pc_error) < 1e-9
    assert maximum_e_error < 1e-12
    removed = set(mesh['excavation_sets']['Pit'])
    volume = sum(p['weight_m3'] for p in layout)
    removed_volume = sum(p['weight_m3'] for p in layout if p['element_index'] in removed)
    assert abs(volume-768) < 1e-8 and abs(removed_volume-16) < 1e-8
    topology = ref.topology(mesh)
    area, area_vector = 0., np.zeros(3)
    for face in mesh['face_sets']['free_surface']:
        q = ref.face(mesh, topology[tuple(sorted(face))][0], face)
        assert abs(q['position'][:, 2]-(10-q['position'][:, 0]/3)).max() < 1e-10
        area += np.linalg.norm(q['area_vector'], axis=1).sum()
        area_vector += q['area_vector'].sum(axis=0)
    assert abs(area-96*math.sqrt(10/9)) < 1e-8
    assert abs(area_vector-[32., 0., 96.]).max() < 1e-8
    return dict(elements=len(mesh['tetrahedra']), integration_points=len(points),
                original_volume_m3=volume, removed_volume_m3=removed_volume,
                loaded_area_m2=float(area), loaded_oriented_area_vector_m2=area_vector.tolist(),
                minimum_input_mean_pressure_kpa=minimum_p,
                maximum_initial_position_error_m=maximum_position_error,
                maximum_initial_stress_recipe_error_kpa=maximum_stress_error,
                maximum_initial_pc_recipe_error_kpa=maximum_pc_error,
                maximum_initial_void_ratio_recipe_error=maximum_e_error,
                seed_recipe='p=20(10-x/3-z); isotropic K0=1; pc/p=4; explicit MCC reference void ratio; zero initial plastic history',
                boundary_recipe='side rollers; original free-surface pressure 0->5->0 kPa; Pit removal with zero remaining pressure')


def record(root, output):
    output.mkdir(parents=True, exist_ok=False)
    freeze_sources.freeze(output)
    geometry = root/'runs/tet10-overlap-verification-20261025-v1'
    studies = root/'runs/mcc-mesh-overlap-20261025-v2'
    verified = read(geometry/'verification.json')
    assert verified['passed'] and len(verified['analytic_cases']) == 4
    assert verified['changed_active_domain_refused'] and verified['curved_geometry_refused'] and verified['coordinate_change_refused']
    assert sha(root/'tools/tet10_overlap.py') == read(geometry/'source-snapshot.json')['tools/tet10_overlap.py']
    summaries, files = [], [geometry/'verification.json']
    for relative in ('runs/projected-mcc-development-20261025-v1/terrain-ultrafine-job',
                     'runs/projected-mcc-spatial-inputs-20261025-v1/mesh-3m/ultrafine-job',
                     'runs/projected-mcc-spatial-inputs-20261025-v1/mesh-2m/ultrafine-job'):
        job = root/relative
        result, status, manifest = [read(job/name) for name in ('result.json', 'status.json', 'job.json')]
        assert result['complete'] and status['state'] == 'succeeded' and status['artifact_sha256'] == sha(job/'result.json')
        assert manifest['binding']['identity'] == result['provenance']
        for name, digest in manifest['binding']['files_sha256'].items():
            assert sha(job/name) == digest, name
        mesh, ground, case = [read(job/'inputs'/name) for name in ('mesh.json', 'ground.json', 'case.json')]
        summaries.append(benchmark(mesh, ground, case, result['integration_points']))
        files += [job/name for name in ('result.json', 'job.json', 'status.json')]
    comparisons = []
    for name in ('rust-224-413.json', 'rust-413-723.json'):
        path = studies/name
        comparison = read(path)
        assert comparison['comparison_recorded'] and not comparison['spatial_convergence_verified']
        for item, digest in comparison['sha256'].items():
            assert sha(Path(item)) == digest, item
        for frame in comparison['frames']:
            for field in ('from_seed', 'from_construction_origin', 'stage_increment'):
                assert frame[field]['maximum_element_coverage_relative_error'] < 1e-8
                assert frame[field]['maximum_element_coverage_error_m3'] < 1e-8
        comparisons.append(dict(left_elements=comparison['left_elements'], right_elements=comparison['right_elements'],
                                 final_from_seed=comparison['frames'][-1]['from_seed'],
                                 final_stage_increment=comparison['frames'][-1]['stage_increment']))
        files.append(path)
    archives = {str(p): archive_check(p) for p in (geometry, studies)}
    assert read(output/'source-snapshot.json') == freeze_sources.manifest()
    report = dict(schema='soil-fem-mcc-overlap-verification-record/1', evidence_integrity_passed=True,
                  analytic_polynomial_cases_passed=4, benchmark_initial_and_boundary_recipes_verified=True,
                  benchmark_summaries=summaries, numerical_mesh_comparisons=comparisons,
                  source_archives_verified=archives, spatial_convergence_verified=False,
                  native_ui_connected=False, goal_complete=False, p6_complete=False,
                  evidence_sha256={str(p): sha(p) for p in files},
                  limitations='Affine meshes only; 723-cell independent field comparison, independent integration refinement, 1792-cell solution and all other P0-P6 work remain pending')
    with (output/'verification.json').open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps(dict(evidence_integrity_passed=True, analytic_cases=4,
                          verified_benchmark_elements=[s['elements'] for s in summaries], spatial_convergence_verified=False)), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    record(Path(__file__).resolve().parents[1], args.output.resolve())
