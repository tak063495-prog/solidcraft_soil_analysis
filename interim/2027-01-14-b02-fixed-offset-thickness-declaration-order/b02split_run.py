"""One-factor B-02 controls on the two existing meshes."""
import argparse
import copy
from pathlib import Path
import sys

FOLLOW = Path(__file__).resolve().parents[1] / '2027-01-02-b02-offset-composite-cloud-followup'
sys.path.insert(0, str(FOLLOW))
from b02follow_common import command, fingerprint, legacy, sha, write

OFFSET = 0.025
THIN = 0.05
CONTROLS = ((THIN, 'default'), (0.1, 'default'), (0.2, 'default'),
            (THIN, 'elements'), (THIN, 'faces'), (THIN, 'plates'))
ELEMENT_ARRAYS = ('tetrahedra', 'tetrahedra10', 'material_ids', 'region_ids')


def inputs(level, thickness, order):
    if type(level) is not int or level not in (0, 1) or (thickness, order) not in CONTROLS:
        raise ValueError('outside the fixed 12/96 one-factor design')
    mesh, ground, stages, _ = legacy['cases'].inputs(level, 4, OFFSET, True)
    system = stages['beam_system']
    for section in system['plate_sections']:
        section['thickness_m'] = thickness
    original = fingerprint([mesh, ground, stages])
    changed_mesh, changed_stages = copy.deepcopy((mesh, stages))
    element_ids = list(range(len(mesh['tetrahedra'])))
    if order == 'elements':
        element_ids.reverse()
        for key in ELEMENT_ARRAYS:
            changed_mesh[key].reverse()
    if order == 'faces':
        for faces in changed_mesh['face_sets'].values():
            faces.reverse()
    if order == 'plates':
        changed_stages['beam_system']['plates'].reverse()
    # Invert every possible declaration map; any unrelated mutation is an error.
    restored_mesh, restored_stages = copy.deepcopy((changed_mesh, changed_stages))
    for key in ELEMENT_ARRAYS:
        restored_mesh[key] = [changed_mesh[key][i] for i in element_ids]
    for key in mesh['face_sets']:
        if order == 'faces':
            restored_mesh['face_sets'][key].reverse()
    if order == 'plates':
        restored_stages['beam_system']['plates'].reverse()
    if (restored_mesh, restored_stages) != (mesh, stages):
        raise ValueError('physical inputs differ after undoing declaration permutation')
    mapping = dict(tetrahedra_declared_to_original=element_ids, order=order,
                   canonical_input_sha256=original, thickness_m=thickness,
                   offset_m=OFFSET, offset_over_thickness=OFFSET/thickness)
    return changed_mesh, ground, changed_stages, mapping


def name(level, thickness, order):
    return f'tet{12 * 8**level}-t{thickness:g}-{order}'


def run(output, observer):
    output.mkdir(parents=True, exist_ok=False)
    plan = []
    for level in (0, 1):
        for thickness, order in CONTROLS:
            mesh, ground, stages, mapping = inputs(level, thickness, order)
            folder = output / name(level, thickness, order)
            for key, obj in [('mesh', mesh), ('ground', ground), ('stages', stages), ('order-map', mapping)]:
                write(folder / (key + '.json'), obj)
            plan.append(dict(name=folder.name, level=level, **mapping))
    write(output / 'plan.json', dict(cases=plan, plate_gl=5, observer_sha256=sha(observer),
                                     driver_sha256=sha(__file__)))
    for case in plan:
        folder = output / case['name']
        base = [folder / (k + '.json') for k in ('mesh', 'ground', 'stages')]
        command([observer, 'solve', *base, folder / 'native'], folder / 'native.log')
        command([sys.executable, FOLLOW / 'b02follow_reference.py', folder,
                 folder / 'independent.jsonl', '--order', '5'], folder / 'independent.log')
        gates = legacy['compare'].compare(folder)
        write(folder / 'existing-gates.json', gates)
        print(folder.name, gates['maximum_native_reassembled_residual'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--observer', type=Path, required=True)
    args = parser.parse_args()
    run(args.output.resolve(), args.observer.resolve())
