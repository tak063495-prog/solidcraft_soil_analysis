"""Additional fixed-increment global audit of the sealed GP replay data."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from bounded_mean_release_control import load


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(tangent, controls, output):
    if output.exists():
        raise ValueError('new audit output required')
    current = json.loads((HERE/'archive-sha256.json').read_text())['entries']
    prior = json.loads((ROOT/'interim/2026-11-10-nc-bounded-refinement/archive-sha256.json').read_text())['entries']
    if sha(tangent) != current['diagnosis/tangent/measurements.json']['sha256']:
        raise ValueError('sealed material replay differs')
    if sha(controls) != prior['native-complement-refined128/controls.json']['sha256']:
        raise ValueError('sealed control differs')
    source = ROOT/'tools/bounded_mean_release_control.py'
    if sha(source) != current['sources/tools/bounded_mean_release_control.py']['sha256']:
        raise ValueError('existing holding-load definition source differs')
    m, c = json.loads(tangent.read_text()), json.loads(controls.read_text())
    B, w, R = (np.asarray(c[k]) for k in ('B', 'weights', 'transform'))
    records = []
    for witness in m['material_replays']['native-complement-refined128']:
        base = witness['points']
        D0 = np.asarray([x['implicit_tangent'] for x in base])
        K0 = np.einsum('qai,qab,qbj,q->ij', B, D0, B, w, optimize=True)
        v0 = np.linalg.svd(K0)[2][-1]
        target = load(*[np.asarray(c[k]) for k in ('final_load_kn', 'holding_kn', 'horizontal_holding_kn', 'vertical_holding_kn')], witness['attempt']['target'])
        def metrics(D):
            K = np.einsum('qai,qab,qbj,q->ij', B, D, B, w, optimize=True)
            _, sv, vr = np.linalg.svd(K)
            eig = np.linalg.eigvals(K)
            j = np.argmin(np.abs(eig))
            return dict(smallest_singular=float(sv[-1]), condition_2=float(sv[0]/sv[-1]),
                        symmetric_part_min=float(np.linalg.eigvalsh((K+K.T)/2)[0]),
                        smallest_magnitude_eigenvalue=dict(real=float(eig[j].real), imaginary=float(eig[j].imag)),
                        relative_global_operator_difference=float(np.linalg.norm(K-K0)/np.linalg.norm(K0)),
                        right_weak_mode_absolute_alignment=float(abs(vr[-1]@v0)))
        row = dict(label=witness['label'], base=metrics(D0), partitions={}, finite_differences=[])
        for count, points in witness['partition_sensitivity'].items():
            if any('error' in x for x in points) or [x['point'] for x in points] != [x['point'] for x in base]:
                raise ValueError('partition GP replay error or order differs')
            D = np.asarray([x['implicit_tangent'] for x in points])
            stress = np.asarray([x['updated_state']['effective_stress_compression_positive_kpa'] for x in points])
            force = -np.einsum('qai,qa,q->i', B, stress, w)
            row['partitions'][count] = dict(**metrics(D),
                fixed_increment_relative_force_residual=float(np.linalg.norm(R.T@target-force)/max(np.linalg.norm(target), 1.)),
                new_global_equilibrium_not_solved=True)
        for index in range(3):
            entries = [x['finite_differences'][index] for x in base]
            row['finite_differences'].append(dict(width=entries[0]['width'],
                changed_branch_points=sum(any(a != point['plastic_substeps'] or b != point['plastic_substeps']
                  for a, b in entry['plastic_substeps_plus_minus']) for point, entry in zip(base, entries)),
                operators={kind: metrics(np.asarray([entry[kind] for entry in entries]))
                           for kind in ('central', 'forward', 'backward')}))
        records.append(row)
    result = dict(passed=False, input_tangent_sha256=sha(tangent), input_controls_sha256=sha(controls),
                  source_sha256=sha(Path(__file__)), records=records,
                  material_partition_acceptance=False, full_release_accepted=False, static_limit_confirmed=False,
                  scope='sealed same-history fixed-increment GP data; not a new equilibrium path')
    output.write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('tangent', 'controls', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    a = p.parse_args()
    run(a.tangent.resolve(), a.controls.resolve(), a.output.resolve())
