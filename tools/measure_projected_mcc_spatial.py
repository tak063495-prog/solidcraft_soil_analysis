"""Area-integrated displacement measures for the planar CAD slope/Pit study.

Uses the same physical floor/surrounding surfaces on every mesh, rather than
mesh-dependent nodal means. This measures sensitivity, not a convergence gate.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import fem_reference as ref
from mcc_projected_volume_diagnostic import read, sha


def operators(mesh):
    topology = ref.topology(mesh)
    removed = set(mesh['excavation_sets']['Pit'])
    xyz = np.asarray(mesh['nodes'])
    rows, audits = {}, {}
    for name, area_expected, centroid_expected in (
            ('Pit:floor', 16., [6., 4., 7.]),
            ('Pit:surrounding_surface', 80*np.sqrt(10/9), [6., 4., 8.])):
        weights = np.zeros(len(xyz))
        area, moment = 0., np.zeros(3)
        for face in mesh['face_sets'][name]:
            owner = next(i for i in topology[tuple(sorted(face))] if i not in removed)
            q = ref.face(mesh, owner, face)
            da = np.linalg.norm(q['area_vector'], axis=1)
            area += da.sum()
            moment += da@q['position']
            np.add.at(weights, q['nodes'], da@q['shape'])
        centroid = moment/area
        errors = dict(area_m2=abs(float(area)-area_expected),
                      centroid_m=float(abs(centroid-centroid_expected).max()),
                      constant_patch=abs(float(weights.sum()/area)-1),
                      affine_patch_m=float(abs(weights@xyz/area-centroid_expected).max()))
        assert max(errors.values()) < 1e-10, (name, errors)
        rows[name] = weights/area
        audits[name] = dict(area_m2=float(area), centroid_m=centroid.tolist(), analytic_geometry_errors=errors)
    return rows, audits


def measure(entries, output):
    reports, bindings = [], {}
    for label, mesh_file, result_file in entries:
        mesh_file, result_file = Path(mesh_file).resolve(), Path(result_file).resolve()
        mesh, result = read(mesh_file), read(result_file)
        bindings.update({str(p): sha(p) for p in (mesh_file, result_file)})
        rows, audit = operators(mesh)
        if result['schema'] in ('soil-fem-projected-mcc-cycle-diagnostic/1', 'soil-fem-independent-mcc-rate-fem/1'):
            for name, digest in result.get('input_and_source_sha256', result.get('input_sha256', {})).items():
                assert sha(Path(name)) == digest, name
            frames = result['frames']
        elif result['schema'] == 'soil-fem-projected-mcc-stage-result/1':
            assert result['complete'], 'incomplete solve cannot be measured as complete'
            frames = [dict(name=f['name'], displacement_m=f['state']['displacement_m']) for f in result['frames']]
        else:
            raise ValueError('explicit projected reference or completed native result required')
        origin = np.asarray(result['reference_displacement_m'])
        assert np.array_equal(origin, frames[0]['displacement_m'])
        previous = origin
        measures = []
        for frame in frames:
            u = np.asarray(frame['displacement_m'])
            measures.append(dict(name=frame['name'], surfaces={
                name: dict(area_mean_from_seed_m=(row@u).tolist(),
                           area_mean_from_construction_origin_m=(row@(u-origin)).tolist(),
                           area_mean_stage_increment_m=(row@(u-previous)).tolist()) for name, row in rows.items()}))
            previous = u
        reports.append(dict(label=label, elements=len(mesh['tetrahedra']), nodes=len(mesh['nodes']),
                            surface_quadrature_audit=audit, frames=measures))
    assert bindings == {name: sha(Path(name)) for name in bindings}
    report = dict(schema='soil-fem-projected-mcc-spatial-measures/1', measured=True,
                  goal_complete=False, p6_complete=False, spatial_convergence_accepted=False,
                  cases=reports, sha256={**bindings, str(Path(__file__).resolve()): sha(Path(__file__))},
                  limitations='Area-weighted sensitivity measures only; field convergence, integrator refinement and independent comparisons still required')
    with output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps([dict(label=r['label'], elements=r['elements'],
                          excavation_floor_mean_increment_mm=1000*r['frames'][-1]['surfaces']['Pit:floor']['area_mean_stage_increment_m'][2]) for r in reports]), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path)
    parser.add_argument('--case', nargs=3, action='append', required=True, metavar=('LABEL', 'MESH', 'RESULT'))
    args = parser.parse_args()
    measure(args.case, args.output.resolve())
