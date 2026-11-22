"""Bounded independent MC/contact comparison, with no native solution as input."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent / '2026-12-18-c03-plastic-bulk-cause-diagnosis'))
from c03diag_study import Model, fem, read, run


def measure(tet10, cutback):
    import numpy as np
    n, opening, phase, until = (32, .12, 'reload', 1.) if tet10 else (4, .3, 'open', .25)
    examples = ROOT / 'examples'
    mesh = read(examples / ('joint-water-column-tet10-mesh.json' if tet10 else 'joint-water-column-mesh.json'))
    ground = fem.normalize_ground(read(examples / 'ground-joint-water-mc.json'))
    for material in ground['materials']:
        material['mohr_coulomb']['cohesion_kpa'] = 30.
    case = fem.normalize_stage_case(read(examples / 'stages-joint-water-cycle.json'))
    for stage in case['stages']:
        stage['initial_increment'] = 1/n
        if stage['name'] == 'open':
            stage['target_factor'] = opening
    case['settings'].update(initial_increment=1/n, material_maximum_substep=1e-4)
    case['joint_system']['maximum_increment_m'] = 1e-4
    model = Model(mesh, ground, case)
    rows, rejected = [], []
    prefix = hashlib.sha256()
    started = time.monotonic()

    class Reached(Exception):
        pass

    def save(row, fields):
        residual = float(np.linalg.norm(model.t.T @ fields['reaction'].ravel()))
        progress = row['progress']
        load = model.base + ((1-progress)*(model.seed_force-model.base)
                             if row['stage'] == 'geostatic' else row['factor']*model.boundary)
        residual /= max(1., np.linalg.norm(load))
        assert residual < 1e-10
        assert all(np.isfinite(value).all() for value in fields.values())
        rows.append(dict(stage=row['stage'], progress=progress, reassembled_force=residual))
        if row['stage'] in ('geostatic', 'shear'):
            payload = {key: value.tolist() for key, value in fields.items()}
            prefix.update(json.dumps([row, payload], allow_nan=False, sort_keys=True).encode())
        if row['stage'] == phase and progress >= until:
            raise Reached

    def record(attempt):
        if not attempt['accepted']:
            rejected.append(attempt)

    try:
        run(model, n, 1e-4, 1e-4, save, cutback=cutback, attempt_sink=record)
    except Reached:
        stop = 'observation point'
    except ArithmeticError as error:
        stop = str(error)
    return dict(steps=n, opening=opening, cutback=cutback, phase=phase, until=until, stop=stop,
                prefix_sha256=prefix.hexdigest(), accepted=len(rows), rejected=rejected,
                last=rows[-1], observed_phase=[r for r in rows if r['stage']==phase],
                maximum_reassembled_force=max(r['reassembled_force'] for r in rows),
                seconds=time.monotonic()-started, numerical_acceptance=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tet10', action='store_true')
    parser.add_argument('--cutback', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(measure(args.tet10, args.cutback), stream, indent=2, allow_nan=False)
        stream.write('\n')
