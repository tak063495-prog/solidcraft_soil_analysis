"""Measure shared physical probes, endpoint work and original reassembly gates."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

from b02split_run import CONTROLS, THIN, name
from b02follow_common import legacy, metric, norm, write
from b02follow_report import compare, observe
import b02cause_report as prior

METHODS = ('native', 'independent')
SPATIAL = ('q_kn_m', 'instantaneous_work_kn_m', 'cumulative_work_kn_m')


def original_forces(row, free):
    if 'state' in row:
        internal = np.asarray(row['soil_internal_kn']).copy()
        internal[np.asarray(row['interface_dofs'], int)] += row['interface_internal_kn']
    else:
        internal = np.asarray(row['total_internal_kn'])
    applied = np.asarray(row.get('state', row)['applied_nodal_load_kn'])
    residual = (internal - applied)[free]
    denom = max(norm(applied), 1.)
    return dict(free_residual_components_kn=residual.tolist(), l2_kn=norm(residual),
                denominator_kn=denom, relative_fixed_order=norm(residual)/denom)


def spatial(coarse, fine):
    result = {}
    for key in SPATIAL:
        a, b = [], []
        for left, right in zip(coarse, fine, strict=True):
            if (left['stage_index'], left['progress']) != (right['stage_index'], right['progress']):
                raise ValueError('unmatched accepted state')
            lval, rval = np.asarray(left[key]).ravel(), np.asarray(right[key]).ravel()
            if lval.shape != rval.shape:
                raise ValueError('unmatched physical probe')
            a.extend(lval)
            b.extend(rval)
        result[key] = metric(a, b)  # Same norm and 96-element denominator as 71655b7.
    return result


def report(run, output):
    if output.exists():
        raise FileExistsError(output)
    cases, records, pairs, contrasts = {}, {}, {}, {}
    accepted = 0
    for case in json.loads((run/'plan.json').read_bytes())['cases']:
        key = case['name']
        folder = run/key
        coarse = run/name(0, case['thickness_m'], 'default')
        mesh = json.loads((folder/'mesh.json').read_bytes())
        fixed = set(mesh['node_sets']['bottom'])
        free = [3*n+i for n in range(len(mesh['nodes'])) if n not in fixed for i in range(3)]
        gates = json.loads((folder/'comparison.json').read_bytes())
        summaries = {}
        for method in METHODS:
            observations, details = observe(folder, coarse, method, case)
            records[key, method] = observations
            rows = prior.rows(folder, method)
            if len(rows) != 20 or len(observations) != 20:
                raise ValueError(f'incomplete history: {key}/{method}')
            accepted += len(rows)
            raw = []
            for row, gate, detail, obs in zip(rows, gates['all_point_component_comparisons'],
                                               details, observations, strict=True):
                point = (row['stage_index'], row['progress'])
                if point != (obs['stage_index'], obs['progress']) or row['accepted_index'] != gate['accepted_index']:
                    raise ValueError('unmatched accepted index')
                raw.append(dict(accepted_index=row['accepted_index'], stage_index=point[0], progress=point[1],
                                original_relative=gate[method+'_reassembled_relative_residual'],
                                structure_force_l2_kn=detail['explicit_free_structure_force_l2_kn'],
                                structure_moment_l2_kn_m=detail['explicit_free_structure_moment_l2_kn_m'],
                                instantaneous_work_kn_m=obs['instantaneous_work_kn_m'],
                                cumulative_work_kn_m=obs['cumulative_work_kn_m'],
                                **original_forces(row, free)))
            write(folder/(method+'-reassembly.json'), raw)
            peak = max(raw, key=lambda v: v['original_relative'])
            summaries[method] = dict(max_relative=peak['original_relative'], peak_accepted_index=peak['accepted_index'],
                                     peak_l2_kn=peak['l2_kn'], peak_denominator_kn=peak['denominator_kn'],
                                     max_structure_force_kn=max(v['structure_force_l2_kn'] for v in raw),
                                     max_structure_moment_kn_m=max(v['structure_moment_l2_kn_m'] for v in raw),
                                     above_original_gate=[v['accepted_index'] for v in raw if v['original_relative'] > 1e-10])
        cases[key] = dict(existing_field_gates_passed=gates['existing_field_gates_passed'], residuals=summaries,
                          input_only_vs_native=compare(records[key,'native'], records[key,'independent']))
        print(key, flush=True)
    for thickness, order in CONTROLS:
        tag = name(1, thickness, order)
        pairs[tag] = {m: spatial(records[name(0, thickness, order),m], records[tag,m]) for m in METHODS}
        for level in (0,1):
            key = name(level, thickness, order)
            baseline = name(level, THIN, 'default')
            contrasts[key] = {m: compare(records[key,m],records[baseline,m]) for m in METHODS}
    write(output, dict(passed=False, b02_complete=False, offset_m=.025, quadrature_gl=5,
                       native_histories=len(cases), independent_histories=len(cases), accepted_states=accepted,
                       cases=cases, spatial_12_to_96=pairs, within_mesh_vs_thin_default=contrasts,
                       original_relative_reassembly_gate=1e-10))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report(args.run.resolve(), args.output.resolve())
