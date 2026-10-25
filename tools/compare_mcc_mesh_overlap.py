"""Volume-integrated displacement differences across nonmatching affine meshes.

Requires completed, bound results and matching configuration headers. This is
a numerical comparison, not a continuum-convergence acceptance. In particular,
equality of mesh-specific initial GP seeds needs separate provenance evidence.
"""
import copy
import hashlib
import json
from pathlib import Path
import argparse
from tet10_overlap import Overlay, np
import tet10_overlap
import mixed_mcc_io
import mcc_projected_cycle_diagnostic as cycle


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(mesh_path, result_path):
    mesh, result = read(mesh_path), read(result_path)
    files = [mesh_path, result_path]
    if result['schema'] == mixed_mcc_io.PROJECTED_SCHEMA:
        folder = result_path.parent
        manifest, status = read(folder/'job.json'), read(folder/'status.json')
        assert status['state'] == 'succeeded' and result['complete']
        assert status['artifact_sha256'] == sha(result_path)
        assert manifest['binding']['identity'] == result['provenance']
        for name, digest in manifest['binding']['files_sha256'].items():
            assert sha(folder/name) == digest, name
        assert sha(folder/'inputs/mesh.json') == sha(mesh_path)
        mixed_mcc_io.validate(mesh, result)
        frames = [dict(name=f['name'], **f['state']) for f in result['frames']]
        ground, case = read(folder/'inputs/ground.json'), read(folder/'inputs/case.json')
        files += [folder/name for name in ('job.json', 'status.json', 'inputs/ground.json', 'inputs/case.json')]
    elif result['schema'] in ('soil-fem-projected-mcc-cycle-diagnostic/1', 'soil-fem-independent-mcc-rate-fem/1'):
        bindings = result.get('input_and_source_sha256', result.get('input_sha256', {}))
        assert bindings
        for name, digest in bindings.items():
            assert sha(Path(name)) == digest, name
        paths = {Path(name).name: Path(name) for name in bindings}
        assert sha(paths['mesh.json']) == sha(mesh_path)
        ground, case = read(paths['ground.json']), read(paths['case.json'])
        frames = result['frames']
        files += [paths['ground.json'], paths['case.json']]
    else:
        raise ValueError('completed projected native job or independent MCC result required')
    configuration = copy.deepcopy(ground)
    geo = configuration['initialization'].get('geostatic')
    mesh_specific_seed = bool(geo and geo['seed']['kind'] == 'integration_points')
    if mesh_specific_seed:
        geo['seed'] = dict(kind='integration_points', note='Mesh-specific point values require separate physical-recipe provenance')
    return mesh, frames, np.asarray(result['reference_displacement_m']), case, configuration, mesh_specific_seed, files


def compare(mesh_a, result_a, mesh_b, result_b, output):
    left, right = load(mesh_a, result_a), load(mesh_b, result_b)
    ma, fa, oa, ca, ga, seed_a, files_a = left
    mb, fb, ob, cb, gb, seed_b, files_b = right
    assert cycle.physical_case(ca) == cycle.physical_case(cb), 'case configuration differs'
    assert ga == gb, 'material/initialization headers differ'
    assert sorted(ma['regions'], key=lambda r: r['id']) == sorted(mb['regions'], key=lambda r: r['id'])
    assert len(fa) == len(fb)
    inputs = files_a+files_b+[Path(__file__).resolve(), Path(tet10_overlap.__file__).resolve(),
                            Path(mixed_mcc_io.__file__).resolve(), Path(cycle.__file__).resolve()]
    bindings = {str(path): sha(path) for path in inputs}
    overlay = Overlay(ma, mb, progress=lambda e, n, pairs: print(f'Overlay {e}/{n}, intersections {pairs}', flush=True))
    frames = []
    previous_a, previous_b = np.zeros_like(oa), np.zeros_like(ob)
    for a, b in zip(fa, fb):
        assert a['name'] == b['name'], 'stage names differ'
        raw = overlay.compare(a['displacement_m'], b['displacement_m'], a['active_elements'], b['active_elements'])
        construction = overlay.compare(np.asarray(a['displacement_m'])-oa, np.asarray(b['displacement_m'])-ob,
                                       a['active_elements'], b['active_elements'])
        stage = overlay.compare(np.asarray(a['displacement_m'])-previous_a, np.asarray(b['displacement_m'])-previous_b,
                                a['active_elements'], b['active_elements'])
        frames.append(dict(name=a['name'], from_seed=raw, from_construction_origin=construction, stage_increment=stage))
        previous_a, previous_b = np.asarray(a['displacement_m']), np.asarray(b['displacement_m'])
    assert bindings == {str(path): sha(path) for path in inputs}
    report = dict(schema='soil-fem-mcc-nonmatching-displacement-integral/1', comparison_recorded=True,
                  left_elements=len(ma['tetrahedra']), right_elements=len(mb['tetrahedra']),
                  geometry=overlay.metadata, frames=frames, sha256=bindings,
                  case_and_material_headers_matched=True,
                  mesh_specific_initial_seed_identity_verified=not(seed_a or seed_b),
                  goal_complete=False, p6_complete=False, spatial_convergence_verified=False,
                  definition='Squared quadratic displacement difference integrated over convex tetrahedron intersections with per-element coverage checks',
                  limitations='Affine geometric cells only; no raw GP stress transfer; initial GP recipe and nodal boundary physical equivalence need separate provenance; nonlinear/time/spatial convergence and continuum error bounds are not certified')
    with output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps([dict(name=f['name'], from_seed_relative_l2=f['from_seed']['relative_l2'],
                          from_seed_rms_m=f['from_seed']['rms_difference_m'],
                          stage_increment_relative_l2=f['stage_increment']['relative_l2'],
                          construction_relative_l2=f['from_construction_origin']['relative_l2']) for f in frames]), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ('mesh_a', 'result_a', 'mesh_b', 'result_b', 'output'):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    compare(*(getattr(args, key).resolve() for key in ('mesh_a', 'result_a', 'mesh_b', 'result_b', 'output')))
