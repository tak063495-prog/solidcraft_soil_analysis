"""Combine disjoint small studies and audit their unchanged raw output files."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from audit import audit
from reference import dump


def report(directories, out):
    combined = None; paths = {}; directories = list(map(Path, directories))
    for directory in directories:
        study = json.loads((directory/'study.json').read_text(encoding='utf-8'))
        if not study['complete']:
            raise ValueError('incomplete study is never successful')
        if combined is None:
            combined = study
        else:
            if set(combined['cases']) & set(study['cases']):
                raise ValueError('overlapping studies; no silently replaced results')
            combined['cases'].update(study['cases']); combined['raw_sha256'].update(study['raw_sha256'])
            combined['elapsed_seconds'] += study['elapsed_seconds']
        for name, expected in study['raw_sha256'].items():
            path = directory/name
            if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                raise ValueError('raw history changed: '+name)
            paths[name] = path
    result = dict(schema='dry-multi-joint-observations/1', complete=True, accepted=None,
                  A_01='未達', NC_01='未達・対象外', P6='未達・対象外',
                  raw_run_count=len(paths), raw_frame_count=0, cases={}, audits={})
    for name, entry in combined['cases'].items():
        base_name = f'{name}-steps32-gp2-identity.json'
        # The report can also handle a different last step from the PC study.
        steps = sorted(map(int, entry['load_steps']))
        base_name = f'{name}-steps{steps[-1]}-gp2-identity.json'
        doc = json.loads(paths[base_name].read_text(encoding='utf-8'))
        result['audits'][name] = audit(doc)
        order = {}
        for variant, row in entry['order'].items():
            order[variant] = dict(discrete_history_equal=row['discrete_history_equal'],
                                  open_state_differences=row['open_state_differences'],
                                  maximum_absolute_by_field={key: v['maximum_absolute'] for key, v in row['differences'].items()},
                                  accepted=None, tolerance='未定義')
        e = dict(reference=entry['reference'], order=order,
                 work_step_series={key: v['final_work_defect_kn_m'] for key, v in entry['load_steps'].items()},
                 work_maximum_series={key: v['maximum_work_defect_kn_m'] for key, v in entry['load_steps'].items()},
                 step_pair_relative_l2=[{key: v['l2_relative'] for key, v in pair['differences'].items()}
                                        for pair in entry['step_pairs']],
                 gp_pair_relative_l2=[{key: v['l2_relative'] for key, v in pair['differences'].items()}
                                      for pair in entry['gp_pairs']], accepted=None)
        if len(entry['step_pairs']) >= 2:
            first, second = entry['step_pairs'][:2]
            # Signed/negative orders are retained. This descriptive statistic
            # uses sampled history relative L2; it is not an extrapolated field.
            ratio = second['coarse']/first['coarse']
            e['observed_history_orders'] = {}
            for field in first['differences']:
                a = first['differences'][field]['l2_relative']; b = second['differences'][field]['l2_relative']
                e['observed_history_orders'][field] = float(np.log(a/b)/np.log(ratio)) if a and b and ratio > 1 else None
        if doc['member_ids']:
            nodes = [n for n in doc['frames'][-1]['nodes'] if n['id'].startswith('B2:')]
            e['right_body_load_sharing_kn'] = {key: np.sum([n[key] for n in nodes], axis=0).tolist()
                                              for key in ('bulk_internal_kn', 'contact_internal_kn', 'structure_internal_kn', 'reaction_kn')}
        result['cases'][name] = e
    for path in paths.values():
        result['raw_frame_count'] += len(json.loads(path.read_text(encoding='utf-8'))['frames'])
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    dump(out/'study.json', combined); dump(out/'observations.json', result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--studies', type=Path, nargs='+', required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); report(a.studies, a.output)


if __name__ == '__main__':
    main()
