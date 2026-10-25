"""Refine only the numerical controls/mesh of the common MC column."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import numpy as np
import fem_reference as ref
from verify_incompressibility import mesh as cube


def prepare(base):
    original = json.loads((base/'mesh-3d.json').read_text(encoding='utf-8'))
    stages = json.loads((base/'stages-3d.json').read_text(encoding='utf-8'))
    records = []
    for n in (2, 4):
        output = base/f'3d-n{n}-fine'
        output.mkdir(exist_ok=False)
        mesh = cube(n, True)
        xyz = np.asarray(mesh['nodes'])*[.2, 1., .2]
        mesh['nodes'] = xyz.tolist()
        mesh['node_sets'] = dict(all=list(range(len(xyz))),
            bottom=np.flatnonzero(xyz[:, 2] == 0.).tolist(),
            top=np.flatnonzero(xyz[:, 2] == .2).tolist())
        mesh['face_sets'] = dict(top=[list(f) for f, owners in ref.topology(mesh).items()
            if len(owners) == 1 and np.all(xyz[list(f), 2] == .2)])
        mesh['regions'] = copy.deepcopy(original['regions'])
        if n == 2:
            assert mesh == original
        case = copy.deepcopy(stages)
        case['settings'].update(initial_increment=1/32, material_maximum_substep=1e-5)
        for name, data in (('mesh.json', mesh), ('stages.json', case)):
            with (output/name).open('x', encoding='utf-8') as stream:
                json.dump(data, stream, indent=2)
        records.append(dict(subdivisions=n, nodes=len(xyz), elements=len(mesh['tetrahedra10']),
            files={name: hashlib.sha256((output/name).read_bytes()).hexdigest()
                   for name in ('mesh.json', 'stages.json')}))
    with (base/'refinement-inputs.json').open('x', encoding='utf-8') as stream:
        json.dump(dict(cases=records, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()), stream, indent=2)
    print(json.dumps(records), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('base', type=Path)
    prepare(parser.parse_args().base.resolve())
