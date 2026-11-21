"""Audit actual failed Newton iterates; never certify removal or destruction.

Input trace uses the existing solver and is not an independent 3D reference.
"""
import argparse
import json
from pathlib import Path
import subprocess

import numpy as np
import probe_low_preload_nc as low
from verify_high_order import prepare


def groups(path):
    current = None
    with path.open(encoding='utf-8') as stream:
        for line in stream:
            event = json.loads(line)
            if event['kind'] == 'progress':
                p = event['event']
                if p['kind'] == 'attempt_started':
                    if current is not None:
                        yield current
                    current = (dict(start=p, events=[], outcome=None)
                               if p['stage_index'] == 3 and .9512958 <= p['progress'] <= .951297 else None)
                elif current is not None and p['kind'] in ('accepted', 'rejected'):
                    current['outcome'] = p
            elif current is not None:
                current['events'].append(event)
    if current is not None:
        yield current


def tangent_report(event):
    k = np.asarray(event['matrix_columns']).T
    rhs = np.asarray(event['reduced_residual'])
    singular = np.linalg.svd(k, compute_uv=False)
    x = np.linalg.solve(k, rhs)
    residual = lambda v: float(np.linalg.norm(k@v-rhs)/max(np.linalg.norm(rhs), 1e-30))
    correction = x.tolist()
    refinement = [residual(x)]
    for _ in range(4):
        x += np.linalg.solve(k, rhs-k@x)
        refinement.append(residual(x))
    return dict(iteration=event['iteration'], projected_load=event['projected_load'],
                dimension=len(k), numpy_correction=correction, singular_values=singular.tolist(),
                condition_number=float(singular[0]/singular[-1]),
                symmetric_eigenvalues=np.linalg.eigvalsh((k+k.T)/2).tolist(),
                numpy_linear_relative_residual=refinement[0],
                numpy_refinement_residuals=refinement,
                interpretation='actual discrete Newton matrix; numpy replay only, no native acceptance')


def measure(folder, inputs, executable, output):
    output.mkdir(parents=True, exist_ok=False)
    low.check_inputs(inputs, excavation_increment=1/128, minimum_increment=1/4194304,
                     dense_lu=True, tangent_perturbation=1e-9)
    run = low.read(folder/'run.json')
    mesh = low.read(inputs/'mesh.json'); ground = low.read(inputs/'ground.json')
    case = low.read(inputs/'stages.json')
    trials = list(groups(folder/'trace.jsonl'))
    observed = [g for g in trials if any(e['kind'] == 'iterate' for e in g['events'])]
    if not observed or not all(g['outcome'] is not None for g in observed):
        raise ValueError('missing complete attempt trace')
    last = observed[-1]
    with (folder/'trace.jsonl').open(encoding='utf-8') as stream:
        final_line = None
        for line in stream:
            if line.strip():
                final_line = line
        terminal = json.loads(final_line) if final_line else {}
    if terminal.get('kind') != 'progress' or last['outcome'] != terminal['event']:
        raise ValueError('terminal failed attempt is outside captured scope')
    if last['outcome']['kind'] != 'rejected':
        raise ValueError('expected final rejected attempt')
    q, offsets = prepare(mesh)
    parameters = {m['id']: m['mcc'] for m in ground['materials']}
    requests = []; mapping = []; line_search = []; matrices = []; partitions = []
    for e in last['events']:
        if e['kind'] == 'iterate':
            old = e['committed']
            begin = len(requests)
            for ei, quadrature in enumerate(q):
                if not e['active'][ei]:
                    continue
                for qi in range(len(quadrature['b'])):
                    pi = int(offsets[ei]+qi)
                    requests.append(dict(parameters=parameters[mesh['material_ids'][ei]],
                                         state=old['material_states'][pi], strain=e['strain_increment'][pi],
                                         partitions=e['partitions'][pi]))
                    mapping.append(dict(iteration=e['iteration'], projected_load=e['projected_load'],
                                        point=pi, state=e['states'][pi]))
            partitions.append(dict(iteration=e['iteration'], projected_load=e['projected_load'],
                                   partitions=e['partitions'], material_begin=begin,
                                   material_end=len(requests), relative_force=e['relative_force_residual']))
        elif e['kind'] == 'tangent':
            matrices.append(tangent_report(e))
        elif e['kind'] in ('line_search', 'repartition'):
            line_search.append({k:v for k,v in e.items() if k not in ('states','displacement_m')})
    low.write(output/'material-requests.json', requests)
    r = subprocess.run([str(executable.resolve()), str(output/'material-requests.json'),
                        str(output/'material-results.json')], capture_output=True, timeout=120, check=False)
    (output/'material.log').write_bytes(r.stdout+r.stderr)
    if r.returncode:
        raise ValueError('material probe error')
    results = low.read(output/'material-results.json')['points']
    point_reports = []
    if len(results) != len(mapping):
        raise ValueError('material count mismatch')
    for native, p in zip(mapping, results, strict=True):
        if 'error' in p:
            raise ValueError('successful native material replay failed')
        if p['updated_state'] != native['state'] or p['local_residuals']['replayed_state'] != native['state']:
            raise ValueError('fixed partition failed exact state replay')
        steps = p['local_residuals']['substeps']
        point_reports.append(dict(iteration=native['iteration'], projected_load=native['projected_load'],
            point=native['point'], plastic_substeps=p['plastic_substeps'],
            normalized_yield=p['normalized_yield'],
            maximum_reconstructed_local_residual=max(s['reconstructed_local_residual_norm'] for s in steps),
            exact_state_replay=True,
            finite_difference_relative_errors=[d['relative_to_implicit'] for d in p['finite_differences']],
            branch_crossing_columns=[sum(a != b for a,b in d['plus_minus_plastic_substeps']) for d in p['finite_differences']]))
    report = dict(schema='soil-fem-nc-failed-newton-diagnostic/1', passed=False,
                  independent_reference=False, complete=run['complete'],
                  highest_progress=run['last_committed']['stage_progress'],
                  accepted_count=low.accepted_count(run), failure_reason=run['failure_reason'],
                  force_tolerance=case['settings']['relative_force_tolerance'],
                  linear_tolerance=case['settings']['relative_linear_tolerance'],
                  traced_attempts=[dict(start=g['start'], outcome=g['outcome'],
                      iterate_count=sum(e['kind']=='iterate' for e in g['events'])) for g in observed],
                  last_attempt=dict(start=last['start'], outcome=last['outcome'],
                      partitions=partitions, line_search=line_search, matrices=matrices,
                      points=point_reports),
                  history_and_force=low.inspect_run(mesh, ground, case, run),
                  limitations=['native-iterate diagnostic, not independent 3D solution',
                               'local residual multiplier reconstructed from exact per-substep states',
                               'no full removal / no destruction or stability acceptance'])
    low.write(output/'measurements.json', report)
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--trace', required=True, type=Path)
    p.add_argument('--inputs', required=True, type=Path)
    p.add_argument('--material-probe', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    a = p.parse_args()
    report = measure(a.trace, a.inputs, a.material_probe, a.output)
    print(json.dumps({k:report[k] for k in ('passed','complete','highest_progress','failure_reason')}, indent=2))
