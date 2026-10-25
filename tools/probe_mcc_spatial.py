"""Sample identical physical positions on different 3D MCC meshes.

Displacements use the actual Tet10 isoparametric inverse/shape functions.
Projected pressure/volume are explicitly primary P1 fields, not GP stress.
Sampling and finite probe errors are sensitivity diagnostics, not convergence
acceptance. Removed points are explicitly inactive rather than interpolated.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import fem_reference as ref
import mixed_mcc_io


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)


def make_points(output, counts):
    if min(counts) < 2 or max(counts) > 50:
        raise ValueError('probe rule counts must be 2..50')
    rules = [np.polynomial.legendre.leggauss(n) for n in counts]
    points, weights = [], []
    for xi, wx in zip(*rules[0]):
        x = 6*(xi+1)
        height = 10-x/3
        for yi, wy in zip(*rules[1]):
            y = 4*(yi+1)
            for ti, wt in zip(*rules[2]):
                points.append([float(x), float(y), float(height*(ti+1)/2)])
                weights.append(float(6*4*height/2*wx*wy*wt))
    volume = float(sum(weights))
    centroid = np.asarray(weights)@np.asarray(points)/volume
    assert abs(volume-768) < 1e-10 and abs(centroid-[5.5, 4., 49/12]).max() < 1e-11
    write(output, dict(schema='soil-fem-fixed-planar-slope-probes/1', points_m=points,
                       original_domain_weights_m3=weights, tensor_rule_counts=counts,
                       original_domain_volume_m3=volume, original_domain_centroid_m=centroid.tolist(),
                       geometry='0<x<12, 0<y<8, 0<z<10-x/3; Pit handled by active elements',
                       limitations='Finite common quadrature probes; filtered post-excavation rule is not exact cut-domain integration'))


def locate(mesh, points):
    xyz = np.asarray(mesh['nodes'])
    cells = np.asarray(ref.cells(mesh))
    geometry = xyz[cells]
    low, high = geometry.min(axis=1), geometry.max(axis=1)
    locations, maximum_patch = [], 0.
    for point in points:
        candidates = np.flatnonzero(np.all((point >= low-1e-9) & (point <= high+1e-9), axis=1))
        found = []
        for e in candidates:
            p = geometry[e]
            natural = np.linalg.solve((p[1:4]-p[0]).T, point-p[0])
            for _ in range(20):
                shape, gradient = ref.basis(natural[None, :], len(p) == 10)
                error = shape[0]@p-point
                if np.linalg.norm(error) < 1e-11:
                    break
                jacobian = (p-p[0]).T@gradient[0]
                natural -= np.linalg.solve(jacobian, error)
            else:
                continue
            if min(*natural, 1-natural.sum()) >= -1e-9:
                maximum_patch = max(maximum_patch, float(abs(shape[0]@p-point).max()))
                found.append((int(e), natural.copy(), shape[0].copy()))
        if not found:
            raise ValueError(f'probe outside supplied mesh: {point.tolist()}')
        locations.append(found)
    assert maximum_patch < 1e-10
    return locations, maximum_patch


def sample(points_file, entries, output):
    probes = read(points_file)
    assert probes['schema'] == 'soil-fem-fixed-planar-slope-probes/1'
    points = np.asarray(probes['points_m'])
    cases = []
    bindings = {str(points_file): sha(points_file)}
    for label, mesh_name, result_name in entries:
        mesh_path, result_path = Path(mesh_name).resolve(), Path(result_name).resolve()
        mesh, result = read(mesh_path), read(result_path)
        bindings.update({str(p): sha(p) for p in (mesh_path, result_path)})
        if result['schema'] == mixed_mcc_io.PROJECTED_SCHEMA:
            assert result['complete']
            mixed_mcc_io.validate(mesh, result)
            frames = [dict(name=f['name'], **f['state']) for f in result['frames']]
        elif result['schema'] in ('soil-fem-projected-mcc-cycle-diagnostic/1', 'soil-fem-independent-mcc-rate-fem/1'):
            source_bindings = result.get('input_and_source_sha256', result.get('input_sha256', {}))
            assert source_bindings
            for name, digest in source_bindings.items():
                assert sha(Path(name)) == digest, name
            meshes = [digest for name, digest in source_bindings.items() if Path(name).name == 'mesh.json']
            assert meshes == [sha(mesh_path)], 'probe mesh does not match reference mesh'
            frames = result['frames']
        else:
            raise ValueError('completed projected native or independent MCC result required')
        locations, patch = locate(mesh, points)
        cells = np.asarray(ref.cells(mesh))
        origin = np.asarray(result['reference_displacement_m'])
        sampled_frames = []
        for frame in frames:
            u = np.asarray(frame['displacement_m'])
            active = frame['active_elements']
            mixed = frame.get('mixed_mcc')
            primary = mixed if mixed is not None else frame
            p = {(v['node_index'], v['material_id']): v['compression_positive_kpa'] for v in primary.get('pressure_dofs', [])}
            volume_data = mixed['projected_volume'] if mixed is not None else frame
            theta = {(v['node_index'], v['material_id']): v['compression_positive'] for v in volume_data.get('volume_dofs', [])}
            samples, construction, pressures, volumes, valid = [], [], [], [], []
            for candidates in locations:
                retained = [item for item in candidates if active[item[0]]]
                if not retained:
                    valid.append(False)
                    for values in (samples, construction, pressures, volumes):
                        values.append(None)
                    continue
                e, natural, shape = retained[0]
                cell = cells[e]
                value = shape@u[cell]
                assert all(abs(other_shape@u[cells[other_e]]-value).max() < 1e-9 for other_e, _, other_shape in retained)
                samples.append(value.tolist())
                construction.append((shape@(u-origin)[cell]).tolist())
                valid.append(True)
                keys = [(int(n), mesh['material_ids'][e]) for n in mesh['tetrahedra'][e]]
                n = np.r_[1-natural.sum(), natural]
                pressures.append(float(n@np.asarray([p[k] for k in keys])) if p else None)
                volumes.append(float(n@np.asarray([theta[k] for k in keys])) if theta else None)
            sampled_frames.append(dict(name=frame['name'], active_probes=valid,
                                       displacement_m=samples, construction_displacement_m=construction,
                                       primary_mean_pressure_kpa=pressures, primary_compression_volume=volumes))
        cases.append(dict(label=label, elements=len(cells), nodes=len(mesh['nodes']),
                          coordinate_interpolation_patch_error_m=patch, frames=sampled_frames))
    assert bindings == {name: sha(Path(name)) for name in bindings}
    write(output, dict(schema='soil-fem-common-mcc-spatial-probes/1', cases=cases,
                       points_m=points.tolist(), original_domain_weights_m3=probes['original_domain_weights_m3'],
                       sha256={**bindings, str(Path(__file__).resolve()): sha(Path(__file__)), str(Path(ref.__file__).resolve()): sha(Path(ref.__file__))},
                       sampled=True, goal_complete=False, p6_complete=False, spatial_convergence_accepted=False,
                       field_definitions='Isoparametric Tet10 displacement; material-split P1 primary pressure and volume; no interpolation or replacement of raw GP stresses',
                       limitations='Finite common probes; pointwise primary pressure can be multi-valued at a material boundary; not a continuum convergence acceptance'))
    print(json.dumps([dict(label=c['label'], elements=c['elements'], active_probes=[sum(f['active_probes']) for f in c['frames']], patch_error_m=c['coordinate_interpolation_patch_error_m']) for c in cases]), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest='command', required=True)
    make = commands.add_parser('make')
    make.add_argument('output', type=Path)
    make.add_argument('--counts', nargs=3, type=int, default=[12, 8, 10])
    run = commands.add_parser('sample')
    run.add_argument('points', type=Path)
    run.add_argument('output', type=Path)
    run.add_argument('--case', nargs=3, action='append', required=True)
    args = parser.parse_args()
    if args.command == 'make':
        make_points(args.output.resolve(), args.counts)
    else:
        sample(args.points.resolve(), args.case, args.output.resolve())
