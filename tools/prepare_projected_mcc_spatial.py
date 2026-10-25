"""Re-mesh the same planar SolidCraft slope/Pit and regenerate point seeds.

This benchmark preserves the original overburden recipe, not the old mesh's
GP inventory. Its K0 seed is explicitly unbalanced and must be equilibrated.
No solved Rust or reference state is transferred between meshes.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import numpy as np
import freeze_sources
import mesh_layers
import make_planar_overburden_seed as seed


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(base, cad, solver, output, sizes):
    if len(set(sizes)) != len(sizes) or min(sizes) <= 0:
        raise ValueError('distinct positive mesh sizes required')
    ground = read(base/'inputs/ground.json')
    mesh = read(base/'inputs/mesh.json')
    manifest = read(cad/'bodies.json')
    # Reject silently changing the physical case or fitting a new seed law.
    if len(ground['materials']) != 1 or len(ground['bodies']) != 1:
        raise ValueError('this spatial benchmark is the single-body planar slope')
    mat = ground['materials'][0]
    assert ground['initialization']['water_level_m'] is None
    assert ground['initialization']['surface_surcharge_kpa'] == 0
    assert mat['unit_weight_kn_m3'] == 20.
    points = ground['initialization']['geostatic']['seed']['points']
    for point in points:
        x, y, z = point['position_m']
        p = 20*(10-x/3-z)
        state = point['state']
        assert np.allclose(state['effective_stress_compression_positive_kpa'], [p, p, p, 0, 0, 0], atol=1e-10, rtol=1e-12)
        assert abs(state['preconsolidation_pressure_kpa']-4*p) < 1e-9
        assert not any(state['plastic_strain_compression_positive'])
        a = mat['mcc']
        e = a['specific_volume_ncl_at_reference']-a['lambda']*np.log(4*p/a['reference_pressure_kpa'])+a['kappa']*np.log(4)-1
        assert abs(e-state['void_ratio']) < 1e-12
    assert [(r['body'], r['source_sha256'], r['volume_m3']) for r in mesh['regions']] == [
        (r['body'], r['source_sha256'], r['volume_m3']) for r in manifest['regions']]
    files = [base/'inputs'/n for n in ('mesh.json', 'ground.json', 'case.json')]
    files += [solver, cad/'bodies.json']+[cad/r['step'] for r in manifest['regions']+manifest['excavation_tools']]
    bindings = {str(p): sha(p) for p in files}
    output.mkdir(parents=True, exist_ok=False)
    freeze_sources.freeze(output)
    geometry = output/'geometry'
    geometry.mkdir()
    shutil.copyfile(cad/'bodies.json', geometry/'bodies.json')
    for r in manifest['regions']+manifest['excavation_tools']:
        shutil.copyfile(cad/r['step'], geometry/r['step'])
    binary = output/'solver/soil-fem-3d.exe'
    binary.parent.mkdir()
    shutil.copyfile(solver, binary)
    reports = []
    for size in sizes:
        inputs = output/f'mesh-{size:g}m'/'inputs'
        inputs.mkdir(parents=True)
        m, identity, g = [inputs/name for name in ('mesh.json', 'mesh-identity.json', 'ground.json')]
        mesh_layers.generate(geometry/'bodies.json', m, size, 2, 3)
        subprocess.run([str(binary), 'mesh-identity', str(m), str(identity)], check=True)
        seed.create(m, identity, base/'inputs/ground.json', g, 10., -1/3, 0., 1., 1., 4.)
        shutil.copyfile(base/'inputs/case.json', inputs/'case.json')
        generated, layout = read(m), read(identity)
        assert read(g)['initialization']['geostatic']['seed']['mesh_sha256'] == layout['mesh_sha256']
        assert generated['excavation_sets']['Pit']
        reports.append(dict(mesh_size_m=size, elements=len(generated['tetrahedra']),
                            nodes=len(generated['nodes']), integration_points=len(layout['integration_points']),
                            removed_elements=len(generated['excavation_sets']['Pit']),
                            inputs_sha256={p.name: sha(p) for p in (m, identity, g, inputs/'case.json')},
                            seed_provenance=read(g.with_suffix('.seed-provenance.json')),
                            mesh_provenance=read(m.with_suffix('.provenance.json'))))
    assert bindings == {str(p): sha(p) for p in files}, 'original input or binary changed'
    report = dict(schema='soil-fem-projected-mcc-spatial-inputs/1', input_preparation_passed=True,
                  goal_complete=False, p6_complete=False, spatial_convergence_accepted=False,
                  original_input_sha256=bindings, meshes=reports,
                  limitations='Inputs only; neither equilibrium nor nonlinear/spatial-convergence acceptance')
    with (output/'preparation.json').open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps([dict(mesh_size_m=r['mesh_size_m'], elements=r['elements'], nodes=r['nodes'], integration_points=r['integration_points']) for r in reports]), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ('base', 'cad', 'solver', 'output'):
        parser.add_argument(name, type=Path)
    parser.add_argument('--sizes', type=float, nargs='+', default=[3., 2.])
    args = parser.parse_args()
    prepare(*(getattr(args, k).resolve() for k in ('base', 'cad', 'solver', 'output')), args.sizes)
