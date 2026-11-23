"""Prepare two fresh runs of the same 96 Tet10 stage problem."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

NAMES = ('mesh', 'ground', 'stages')


def encoded(value):
    return (json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n').encode()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('output', type=Path)
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument('--source', type=Path, help='original PC mesh.json, ground.json, stages.json')
    group.add_argument('--cloud-original', action='store_true')
    group.add_argument('--cloud-face-sets', action='store_true')
    args = p.parse_args()
    if args.source:
        raw = [(args.source / (name + '.json')).read_bytes() for name in NAMES]
    else:
        root = Path(__file__).resolve().parents[2]
        sys.path.insert(0, str(root / 'tools'))
        from refine_tet10 import refine
        if args.cloud_face_sets:
            from nc_collapse_sweep import profile
            mesh, ground, case = profile(1.4, 100, True)
        else:
            from probe_nc_midpoint import inputs
            mesh, ground, case = inputs('midpoint')
        mesh, _ = refine(mesh)
        case['settings']['maximum_free_dofs'] = 3 * len(mesh['nodes'])
        raw = [encoded(value) for value in (mesh, ground, case)]
    mesh, _, case = [json.loads(data) for data in raw]
    if len(mesh['tetrahedra10']) != 96 or len(mesh['tetrahedra']) != 96:
        p.error('expected 96 quadratic tetrahedra')
    settings = case['settings']
    if (settings['relative_linear_tolerance'] != 1e-10
            or settings.get('linear_solver', {'method': 'dense_lu'}) != {'method': 'dense_lu'}
            or settings.get('dense_lu_refinement', False) is not False
            or 'excavate_upper' not in [s['name'] for s in case['stages']]):
        p.error('expected the baseline NC excavation case with dense LU and tolerance 1e-10')
    args.output.mkdir(parents=True, exist_ok=False)
    for name, data in zip(NAMES, raw):
        (args.output / (name + '.json')).write_bytes(data)
        print(name, hashlib.sha256(data).hexdigest())
    settings['dense_lu_refinement'] = True
    (args.output / 'stages-refined.json').write_bytes(encoded(case))


if __name__ == '__main__':
    main()
