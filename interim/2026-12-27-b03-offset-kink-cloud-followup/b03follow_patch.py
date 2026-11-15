"""Constant membrane/curvature patch residuals; no extra equilibrium solves."""
from collections import Counter
import json
from pathlib import Path
import numpy as np
from b03follow_math import mv, norm


def audit(folder):
    rows = []
    eps = np.array([[2e-5, 7e-6, 0.], [7e-6, -1e-5, 0.], [0., 0., 0.]])
    curve = np.array([[3e-5, 1e-5, 0.], [1e-5, -2e-5, 0.], [0., 0., 0.]])
    for path in sorted((Path(folder)/'assemblies').glob('*.json')):
        record = json.loads(path.read_text())
        elements = record['elements']
        edges = Counter(tuple(sorted((nodes[i], nodes[j]))) for nodes in elements for i, j in ((0, 1), (1, 2), (2, 0)))
        boundary = set()
        for nodes in elements:
            for i, j, mid in ((0, 1, 3), (1, 2, 4), (2, 0, 5)):
                if edges[tuple(sorted((nodes[i], nodes[j])))] == 1:
                    boundary.update((nodes[i], nodes[j], nodes[mid]))
        interior = [6*n+a for n in range(len(record['positions'])) if n not in boundary for a in range(6)]
        rotation = np.asarray(record['rotation'])
        for side in ('native', 'independent'):
            fields = record[side+'_fields']['polynomial']
            strains, curves, shear, drilling = [], [], [], []
            for field in fields:
                axes = np.asarray(field['local_axes'])
                for sample in field['samples']:
                    for key, target, dest in (('membrane_strain', eps, strains), ('curvature_per_m', curve, curves)):
                        xx, yy, xy = sample[key]
                        local = np.array([[xx, xy/2, 0.], [xy/2, yy, 0.], [0., 0., 0.]])
                        dest.append(rotation.T@axes.T@local@axes@rotation-target)
                    shear.extend(sample['transverse_shear_strain'])
                    drilling.append(sample['drilling_mismatch_rad'])
            internal = mv(record[side+'_stiffness'], record['polynomial'])
            rows.append(dict(id=record['id'], side=side, interior_dofs=len(interior),
                constant_membrane_strain_tensor_max_abs=float(np.max(np.abs(strains))),
                constant_curvature_tensor_max_abs_per_m=float(np.max(np.abs(curves))),
                shear_should_be_zero_max_abs=max(abs(v) for v in shear),
                drilling_should_be_zero_max_abs_rad=max(abs(v) for v in drilling),
                interior_force_residual_l2_kn=norm(internal[[i for i in interior if i%6<3]]),
                interior_moment_residual_l2_kn_m=norm(internal[[i for i in interior if i%6>=3]])))
    return dict(diagnostic_only=True, new_solves=0, prescribed_fields=90,
        literature_patch_test_reproduced=False,
        reference='Analytic constant membrane strain and curvature of the existing flat polynomial field; not a new external benchmark or tolerance.',
        cases=rows)
