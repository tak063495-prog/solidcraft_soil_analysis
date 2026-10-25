"""Compare fixed physical probes without treating sample norms as acceptance."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare(source, left, right, output):
    data = json.loads(source.read_text(encoding='utf-8-sig'))
    assert data['schema'] == 'soil-fem-common-mcc-spatial-probes/1'
    for name, digest in data['sha256'].items():
        assert sha(Path(name)) == digest, name
    cases = {c['label']: c for c in data['cases']}
    assert len(cases) == len(data['cases']) and left != right
    a, b = cases[left], cases[right]
    assert len(a['frames']) == len(b['frames'])
    weights = np.asarray(data['original_domain_weights_m3'])
    reports = []
    for af, bf in zip(a['frames'], b['frames']):
        assert af['name'] == bf['name'] and af['active_probes'] == bf['active_probes']
        active = np.asarray(af['active_probes'], dtype=bool)
        w = weights[active]
        row = dict(name=af['name'], active_probes=int(active.sum()),
                   filtered_probe_weight_sum_m3=float(w.sum()), fields={})
        for name in ('displacement_m', 'construction_displacement_m', 'primary_mean_pressure_kpa', 'primary_compression_volume'):
            av, bv = [v for v, keep in zip(af[name], active) if keep], [v for v, keep in zip(bf[name], active) if keep]
            if any(v is None for v in av+bv):
                row['fields'][name] = None
                continue
            av, bv = np.asarray(av), np.asarray(bv)
            assert np.isfinite(av).all() and np.isfinite(bv).all()
            difference = av-bv
            if av.ndim == 1:
                difference, bv = difference[:, None], bv[:, None]
            squared = np.square(difference).sum(axis=1)
            reference = np.square(bv).sum(axis=1)
            row['fields'][name] = dict(
                weighted_probe_rms=float(np.sqrt(w@squared/w.sum())),
                weighted_probe_relative_l2=float(np.sqrt(w@squared/max(float(w@reference), 1e-60))),
                maximum_probe_component_difference=float(abs(difference).max()))
        reports.append(row)
    report = dict(schema='soil-fem-mcc-common-probe-comparison/1', left=left, right=right,
                  left_elements=a['elements'], right_elements=b['elements'], frames=reports,
                  goal_complete=False, p6_complete=False, spatial_convergence_accepted=False,
                  sha256={str(source): sha(source), str(Path(__file__).resolve()): sha(Path(__file__))},
                  limitations='Finite common probe norms; need spatial and sample-rule refinement; primary pressure is distinct from GP material stress')
    with output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(json.dumps([dict(name=r['name'], displacement=r['fields']['displacement_m']) for r in reports]), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('source', type=Path)
    parser.add_argument('left')
    parser.add_argument('right')
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    compare(args.source.resolve(), args.left, args.right, args.output.resolve())
