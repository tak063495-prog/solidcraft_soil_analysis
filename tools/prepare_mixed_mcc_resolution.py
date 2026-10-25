"""Refine the existing independent-input MCC column without changing physics.

The copied solver and local reference sources stay immutable during a study.
No previously computed displacement, stress or history is an input.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import shutil
from prepare_mixed_mcc_reference import column


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def write(path, value):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)


def prepare(baseline, solver, output, resolutions):
    root = Path(__file__).resolve().parents[1]
    if not resolutions or any(n < 2 or n % 2 for n in resolutions):
        raise ValueError('positive even vertical subdivisions are required')
    ground = read(baseline/'ground.json')
    case = read(baseline/'n4/stages-finest.json')
    before = copy.deepcopy(case)
    output.mkdir(exist_ok=False)
    shutil.copyfile(baseline/'ground.json', output/'ground.json')
    (output/'solver').mkdir()
    shutil.copyfile(solver, output/'solver/soil-fem-3d.exe')
    reference = output/'reference-sources'
    reference.mkdir()
    # Include the transitive point_state helper, even though the independent
    # numerical kernel only imports this function from the export module.
    for name in ('mcc_rate_fem.py', 'mcc_rate_midpoint.py', 'fem_reference.py',
                 'verify_high_order.py', 'high_order_io.py'):
        shutil.copyfile(root/'tools'/name, reference/name)
    for n in resolutions:
        where = output/f'n{n}'
        where.mkdir()
        write(where/'mesh.json', column(n))
        write(where/'stages-finest.json', case)
        refined = copy.deepcopy(case)
        refined['settings'].update(initial_increment=1/128, material_maximum_substep=1e-7)
        for name in ('initial_increment', 'material_maximum_substep'):
            refined['settings'].pop(name)
        original = copy.deepcopy(before)
        for name in ('initial_increment', 'material_maximum_substep'):
            original['settings'].pop(name)
        assert refined == original, 'refinement changed physical inputs'
        refined = copy.deepcopy(case)
        refined['settings'].update(initial_increment=1/128, material_maximum_substep=1e-7)
        write(where/'stages-ultrafine.json', refined)
    assert ground == read(output/'ground.json')
    write(output/'preparation.json', dict(
        schema='soil-fem-mixed-mcc-resolution-preparation/1',
        development_only=True, goal_complete=False, p6_complete=False,
        physical_origin=str(baseline),
        unchanged_physics='Dry K0=1 NC two-layer column, pressure 100->130->100 kPa; upper half removed, remaining pressure 100 kPa',
        state_origin='Only original K0/OCR/NCL input; no solved Rust or reference field',
        vertical_subdivisions=resolutions,
        baseline_inputs_sha256={str(p):sha(p) for p in (baseline/'ground.json', baseline/'n4/stages-finest.json')},
        sources_sha256={str(p):sha(p) for p in (Path(__file__), root/'tools/prepare_mixed_mcc_reference.py')},
        created_inputs_sha256={p.relative_to(output).as_posix():sha(p) for p in output.rglob('*') if p.is_file()}))
    print('Prepared mixed MCC spatial refinements:', resolutions, flush=True)
    print('Immutable solver SHA256:', sha(output/'solver/soil-fem-3d.exe'), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('baseline', type=Path)
    parser.add_argument('solver', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--resolutions', type=int, nargs='+', default=[4, 8, 16, 32])
    args = parser.parse_args()
    prepare(args.baseline.resolve(), args.solver.resolve(), args.output.resolve(), args.resolutions)
