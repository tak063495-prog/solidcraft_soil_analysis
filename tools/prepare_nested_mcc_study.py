"""Prepare fresh original-recipe MCC inputs on a nested planar slope mesh.

No solved displacement or material history is transferred to the refined mesh.
The old explicit GP seed recipe is checked before generating the new inventory.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import numpy as np
import freeze_sources
import make_planar_overburden_seed as seed


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(base, mesh, solver, case, output):
    audit_path = mesh.with_suffix('.refinement.json')
    audit, ground, original_mesh = map(read, (audit_path, base/'inputs/ground.json', base/'inputs/mesh.json'))
    assert audit['refined'] and audit['child_elements'] == 8*audit['parent_elements']
    assert sha(base/'inputs/mesh.json') in audit['sha256'].values(), 'refinement parent differs from base input'
    for name, digest in audit['sha256'].items():
        assert sha(Path(name)) == digest, name
    if len(ground['materials']) != 1 or len(ground['bodies']) != 1:
        raise ValueError('planar single-body/material benchmark required')
    assert ground['initialization']['water_level_m'] is None and ground['initialization']['surface_surcharge_kpa'] == 0
    material = ground['materials'][0]
    assert material['unit_weight_kn_m3'] == 20.
    for point in ground['initialization']['geostatic']['seed']['points']:
        x, y, z = point['position_m']
        p = 20*(10-x/3-z)
        state = point['state']
        assert np.allclose(state['effective_stress_compression_positive_kpa'], [p, p, p, 0, 0, 0], atol=1e-10, rtol=1e-12)
        assert abs(state['preconsolidation_pressure_kpa']-4*p) < 1e-9
        assert not any(state['plastic_strain_compression_positive'])
        law = material['mcc']
        e = law['specific_volume_ncl_at_reference']-law['lambda']*np.log(4*p/law['reference_pressure_kpa'])+law['kappa']*np.log(4)-1
        assert abs(e-state['void_ratio']) < 1e-12
    refined_mesh = read(mesh)
    assert refined_mesh['regions'] == original_mesh['regions']
    files = [base/'inputs/mesh.json', base/'inputs/ground.json', mesh, audit_path, solver, case]
    bindings = {str(p): sha(p) for p in files}
    output.mkdir(parents=True, exist_ok=False)
    freeze_sources.freeze(output)
    inputs, binary_folder = output/'inputs', output/'solver'
    inputs.mkdir()
    binary_folder.mkdir()
    binary = binary_folder/'soil-fem-3d.exe'
    shutil.copyfile(solver, binary)
    shutil.copyfile(mesh, inputs/'mesh.json')
    shutil.copyfile(audit_path, output/'parent-refinement.json')
    shutil.copyfile(base/'inputs/case.json', inputs/'case.json')
    shutil.copyfile(case, inputs/'case-ultrafine.json')
    identity = inputs/'mesh-identity.json'
    subprocess.run([str(binary), 'mesh-identity', str(inputs/'mesh.json'), str(identity)], check=True)
    seed.create(inputs/'mesh.json', identity, base/'inputs/ground.json', inputs/'ground.json',
                10., -1/3, 0., 1., 1., 4.)
    layout, generated = read(identity), read(inputs/'ground.json')
    assert generated['initialization']['geostatic']['seed']['mesh_sha256'] == layout['mesh_sha256']
    assert len(layout['integration_points']) == 8*len(read(base/'inputs/mesh-identity.json')['integration_points'])
    assert bindings == {str(p): sha(p) for p in files}
    report = dict(schema='soil-fem-nested-mcc-input-preparation/1', input_preparation_passed=True,
                  parent_elements=audit['parent_elements'], elements=len(refined_mesh['tetrahedra']),
                  nodes=len(refined_mesh['nodes']), integration_points=len(layout['integration_points']),
                  excavation_volume_m3=audit['excavation_volumes_m3'],
                  original_input_sha256=bindings,
                  prepared_sha256={str(p): sha(p) for p in (inputs/'mesh.json', inputs/'ground.json', inputs/'case.json',
                                                          inputs/'case-ultrafine.json', identity, binary)},
                  goal_complete=False, p6_complete=False, spatial_convergence_verified=False,
                  limitations='Inputs only; K0 seed is unbalanced; equilibrium and nonlinear comparisons remain required')
    with (output/'preparation.json').open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps({k: report[k] for k in ('elements', 'nodes', 'integration_points', 'excavation_volume_m3')}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ('base', 'mesh', 'solver', 'case', 'output'):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    prepare(*(getattr(args, k).resolve() for k in ('base', 'mesh', 'solver', 'case', 'output')))
