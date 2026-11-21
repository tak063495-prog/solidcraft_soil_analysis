"""Recheck extracted evidence with the recorded source checkout, without solving."""
import argparse
import json
import os
from pathlib import Path
import sys
for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools'))
import fem_reference as fem
import probe_low_preload_nc as low
import verify_wet_structures as wet


def measure(folder):
    nc = folder/'nc'
    inputs = nc/'input'
    low.check_inputs(inputs, excavation_increment=1/128, minimum_increment=1/4194304,
                     dense_lu=True, tangent_perturbation=1e-9)
    mesh, ground, case = [low.read(inputs/(n+'.json')) for n in low.NAMES]
    run = low.read(nc/'experiment/run.json')
    summary = low.inspect_run(mesh, ground, case, run)
    assert not run['complete'] and summary['stage_progress'] < 1
    assert summary == low.read(nc/'nc-experiment-measured.json')
    soil = folder/'soil'
    mesh = wet.read(soil/'mesh.json')
    ground = fem.normalize_ground(wet.read(soil/'ground.json'))
    case = fem.normalize_stage_case(wet.read(soil/'stages.json'))
    native = wet.read(soil/'stage-results.json')
    expected = wet.reference(mesh, ground, case, water_modes=('soil',))
    def serial(value):
        if isinstance(value, np.ndarray): return value.tolist()
        if isinstance(value, np.generic): return value.item()
        if isinstance(value, dict): return {str(k): serial(v) for k, v in value.items()}
        if isinstance(value, (list, tuple)): return [serial(v) for v in value]
        return value
    fields = [{k: v for k, v in row.items() if k != 'darcy'} for row in expected]
    assert serial(fields) == wet.read(soil/'independent-3d-fields.json')
    rows = wet.compare(mesh, ground, case, native, expected)
    saved = wet.read(soil/'verification.json')
    assert rows == saved['frames']
    fine = wet.reference(mesh, ground, case, order=7, water_modes=('soil',))
    difference = max(float(np.linalg.norm(a['u']-b['u'])/max(np.linalg.norm(b['u']), 1e-12))
                     for a, b in zip(expected, fine, strict=True))
    assert difference == saved['independent_gl5_to_gl7_relative_difference']
    assert native['last_committed'] == wet.read(soil/'clean-result.json')['last_committed']
    for row in saved['checkpoints']:
        index = saved['checkpoints'].index(row)
        resumed = wet.read(soil/f'resumed-{index}.json')
        assert resumed['complete']
        for key in ('last_committed', 'frames', 'flow_audits', 'excavation_audits', 'transition_audits', 'attempts'):
            assert wet.without_time(resumed.get(key)) == wet.without_time(native.get(key))
    import meshio
    water = wet.Water(mesh, ground, case)
    point_soil = np.asarray([water.soil[p['element_index']] for p in native['integration_points']])
    for index, ref in enumerate(expected):
        for point in native['frames'][index]['state']['material_states']:
            assert point['equivalent_plastic_strain'] == point['plastic_work_kpa'] == 0
            assert not np.asarray(point['plastic_strain_compression_positive']).any()
        gp = meshio.read(soil/'exports'/f'integration-points-{index:03}.vtu').point_data
        ids = gp['integration_point_id']
        np.testing.assert_array_equal(gp['material_stress_reference'], (~point_soil[ids]).astype(np.uint8))
        np.testing.assert_allclose(gp['pore_pressure_kpa'], ref['pore'][ids], rtol=1e-10, atol=1e-8)
        total = ref['stress'][ids].copy(); total[:, :3] += ref['pore'][ids, None]
        np.testing.assert_allclose(gp['total_stress_compression_positive_kpa'], total, rtol=1e-7, atol=1e-7)
    return dict(schema='soil-fem-nc-soil-remeasurement/1', nc_report_exact=True,
                soil_frame_rows_exact=True, independent_3d_fields_exact=True, gl5_gl7_exact=True, saved_restart_checks=8,
                exported_role_pore_total_stress_checked=True, all_soil_case_plastic_histories_zero=True, nc_passed=False,
                limited_soil_passed=True, p3b_p4_rows_complete=False, p6_complete=False)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--evidence', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    a = p.parse_args()
    report = measure(a.evidence)
    with a.output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2)
    print(json.dumps(report, indent=2))
