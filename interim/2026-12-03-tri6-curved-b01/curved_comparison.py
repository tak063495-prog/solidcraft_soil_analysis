"""Fresh draft Q4 solves and same-physical-vertex diagnostics; no new gates."""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import importlib.util
import json
import math
from pathlib import Path
import numpy as np
from b01_cases import HERE, read, write


def run(study, out):
    out.mkdir(parents=True, exist_ok=False)
    spec = importlib.util.spec_from_file_location('b01_draft_curved', HERE/'draft-reference/reference.py')
    curved = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(curved)
    case = next(c for c in curved.catalog()['benchmarks'] if c['id'] == 'scordelis_lo')
    rows, attached = [], []
    for record in read(study/'measurements.json')['cases']:
        folder = study/record['case']
        if 'comparison' not in record:
            continue
        c = read(folder/'stages.json')
        states = [json.loads(s) for s in (folder/'native/accepted.jsonl').read_text().splitlines()]
        full = next(s for s in states if s['stage_index'] == 2 and s['progress'] == 1.)
        state = full['state']
        n = record['n']
        native = np.asarray(state['beams']['node_displacement_m'])
        rotation = np.asarray(state['beams']['node_rotation_rad'])
        if record['kind'] == 'attached':
            q = float(native[n//2, 2])
            attached.append(dict(n=n, crown_y0_uz_m=q,
                                 soil_displacement_m=state['displacement_m'],
                                 soil_stress_kpa=[v['effective_stress_compression_positive_kpa'] for v in state['material_states']],
                                 native_load_work_kn_m=record['comparison']['maximum_by_field']['load_work_kn_m'],
                                 accepted_index=full['accepted_index']))
            continue
        result, x, u = curved.solve(case, n)
        curved.write_nodes(out/f'curved-roof-n{n}.csv', x, u)
        xyz = np.asarray([v['position_m'] for v in c['beam_system']['nodes']])
        # Both use precisely the same analytic vertices; facet mids are excluded
        # because they are different physical points and have no Q4 node ID.
        assert xyz[:len(x)].shape == x.shape
        assert np.max(np.abs(xyz[:len(x)]-x)) < case['probe']['coordinate_tolerance']
        delta = native[:len(x)]-u
        q = float(native[n, 2])
        rows.append(dict(n=n, native_signed_probe=q, draft_signed_probe=result['displacement'],
                         literature_signed_probe=case['reference']['signed_value'],
                         native_minus_draft=q-result['displacement'],
                         native_relative_literature_error=abs((q+.3024)/.3024),
                         draft_result=result, same_vertex_displacement_difference=delta.tolist(),
                         same_vertex_native_rotations=rotation[:len(x)].tolist(),
                         same_vertex_displacement_relative_l2=float(np.linalg.norm(delta)/np.linalg.norm(u)),
                         same_vertex_maximum_absolute=float(np.max(abs(delta))),
                         curve_acceptance=None, b01_complete=False))
    for collection, field in ((rows, 'native_signed_probe'), (attached, 'crown_y0_uz_m')):
        for i, row in enumerate(collection):
            row['previous_level_signed_difference'] = None
            row['three_level_difference_order'] = None
            if i:
                row['previous_level_signed_difference'] = row[field]-collection[i-1][field]
            if i >= 2:
                d1 = collection[i-1][field]-collection[i-2][field]
                d2 = row[field]-collection[i-1][field]
                if d1*d2 > 0 and abs(d2) < abs(d1):
                    row['three_level_difference_order'] = math.log(abs(d1/d2))/math.log(row['n']/collection[i-1]['n'])
            if collection is attached and i:
                for key in ('soil_displacement_m', 'soil_stress_kpa'):
                    a = np.asarray(collection[i-1][key])
                    b = np.asarray(row[key])
                    row[key+'_previous_level_relative_l2'] = float(np.linalg.norm(a-b)/max(np.linalg.norm(b), 1e-12))
                    row[key+'_previous_level_maximum_absolute'] = float(np.max(abs(a-b)))
    write(out/'curved-comparison.json', dict(passed=False, b01_complete=False, acceptance=None,
          reference_kind='draft PR20 Q4 analytic directors; literature is not an exact solution',
          native_reference_pair='full-load same physical vertices only; all-field Tri6 comparison is separate',
          roof=rows, attached=attached))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('study', type=Path)
    p.add_argument('output', type=Path)
    a = p.parse_args()
    run(a.study.resolve(), a.output.resolve())
