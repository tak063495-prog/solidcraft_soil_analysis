"""Version-2 display fields against real MCC VTUs, including two-material jumps.

This verifies frame extraction only. Native capture, authoring, packet loading
and GPU rendering of the new result family are separate, unfinished work.
"""
import argparse
import copy
import json
from pathlib import Path
import sys
import meshio
import numpy as np
import export_stages
from export_solidcraft_view import scene_frame
from freeze_sources import freeze
from native_workflow_snapshot import sha


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def verify(study, output):
    root = Path(__file__).resolve().parents[1]
    output.mkdir(parents=True, exist_ok=False)
    freeze(output)
    cases = [('terrain', study/'inputs/mesh.json', study/'terrain-job/result.json'),
        ('two-material', root/'examples/hardening-column-tet10-mesh.json', study/'two-material-result.json')]
    reports = []
    for name, mesh_path, result_path in cases:
        mesh, run = read(mesh_path), read(result_path)
        target = output/name
        target.mkdir()
        export_stages.export(mesh_path, result_path, target/'vtu')
        frames = []
        for index, frame in enumerate(run['frames']):
            view, points = scene_frame(mesh, run, index)
            assert view['schema'] == 'soil-fem-view-frame/2'
            assert points['schema'] == 'soil-fem-view-integration-points/2'
            with (target/f'frame-{index:03}.json').open('x', encoding='utf-8') as stream:
                json.dump(dict(frame=view, points=points), stream, allow_nan=False)
            vtk = meshio.read(target/f'vtu/stage-{index:03}.vtu')
            gp = meshio.read(target/f'vtu/integration-points-{index:03}.vtu')
            primary_vtu = meshio.read(target/f'vtu/mean-pressure-dofs-{index:03}.vtu')
            for field, values in view['cell_fields'].items():
                key = field if field in vtk.cell_data else 'volume_average_'+field
                assert key in vtk.cell_data, field
                assert np.allclose(np.asarray(values, float), vtk.cell_data[key][0], rtol=1e-12, atol=1e-12, equal_nan=True), field
            for field, values in points['fields'].items():
                assert field in gp.point_data, field
                assert np.allclose(np.asarray(values, float), gp.point_data[field], rtol=1e-13, atol=1e-13, equal_nan=True), field
            primary = view['mixed_primary']
            for field, vtk_field in [('original_node_ids', 'original_node_id'), ('material_ids', 'material_id'),
                ('mean_effective_stress_compression_positive_kpa', 'mean_effective_stress_compression_positive_kpa'),
                ('cumulative_material_volume_strain_compression_positive', 'cumulative_material_volume_strain_compression_positive'),
                ('frozen_compatibility_offset_m3', 'frozen_compatibility_offset_m3')]:
                assert np.array_equal(primary[field], primary_vtu.point_data[vtk_field]), field
            assert np.array_equal(primary['positions_m'], primary_vtu.points)
            assert len(set(zip(primary['original_node_ids'], primary['material_ids']))) == len(primary['original_node_ids'])
            p = np.asarray(points['fields']['primary_mean_effective_stress_compression_positive_kpa'])
            material = np.asarray(points['fields']['material_mean_effective_stress_compression_positive_kpa'])
            work = np.asarray(points['fields']['work_conjugate_effective_stress_compression_positive_kpa'])
            raw = np.asarray(points['fields']['effective_stress_compression_positive_kpa'])
            assert np.allclose(work[:,:3].mean(axis=1), p, rtol=1e-12, atol=1e-12)
            assert np.array_equal(raw[:,:3].mean(axis=1), material)
            assert np.array_equal(work[:,3:], raw[:,3:])
            jumps = {}
            for node, value in zip(primary['original_node_ids'], primary['mean_effective_stress_compression_positive_kpa']):
                jumps.setdefault(node, []).append(value)
            jump = max((max(v)-min(v) for v in jumps.values()), default=0.)
            frames.append(dict(index=index, material_primary_gap_kpa=float(abs(p-material).max()),
                primary_material_boundary_jump_kpa=float(jump), raw_points=len(p), primary_rows=len(primary['original_node_ids'])))
        if name == 'terrain':
            assert max(f['material_primary_gap_kpa'] for f in frames) > .1
        if name == 'two-material':
            assert max(f['primary_material_boundary_jump_kpa'] for f in frames) > .1
        refused = []
        for label, mutate in [
            ('relabel-as-legacy', lambda d: d.update(schema='soil-fem-stage-result/6')),
            ('missing-volume-row', lambda d: d['frames'][-1]['state']['mixed_mcc']['projected_volume']['volume_dofs'].pop()),
            ('changed-frozen-offset', lambda d: d['frames'][-1]['state']['mixed_mcc']['projected_volume']['compatibility_offsets'][0].update(integral_m3=10.)),
            ('changed-work-force', lambda d: d['frames'][-1].update(relative_force_residual=.1))]:
            damaged = copy.deepcopy(run)
            mutate(damaged)
            try:
                scene_frame(mesh, damaged, len(damaged['frames'])-1)
            except (ValueError, AssertionError):
                refused.append(label)
            else:
                raise AssertionError('corrupted view source accepted: '+label)
        reports.append(dict(name=name, frames=frames, refused=refused,
            mesh_sha256=sha(mesh_path), result_sha256=sha(result_path)))
        print(json.dumps(dict(case=name, frames=len(frames), refused=refused)), flush=True)
    report = dict(schema='soil-fem-projected-view-fields-verification/1', passed=True, cases=reports,
        scope='Version-2 frame extraction and VTU comparison only',
        native_authoring_connected=False, native_packet_connected=False, native_viewport_verified=False,
        independent_physics_acceptance=False, goal_complete=False, p6_complete=False,
        source_manifest_sha256=sha(output/'source-snapshot.json'))
    with (output/'verification.json').open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps(dict(passed=True, cases=len(reports))), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('study', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    verify(args.study.resolve(), args.output.resolve())
