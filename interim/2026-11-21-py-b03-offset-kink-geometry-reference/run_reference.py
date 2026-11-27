"""Rebuild all small, independent B-03 diagnostics; exit 2 means acceptance undefined."""
import argparse
from pathlib import Path

from geometry_reference import (compare_documents, dump, read_document, series,
                                synthetic_document, write_document)
from legacy_refusals import audit_repository


def run(output, repository):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    diagnostics = series()
    dump(diagnostics, output/'series.json')
    legacy = audit_repository(repository)
    dump(legacy, output/'legacy-refusals.json')
    fixture = synthetic_document()
    write_document(fixture, output/'synthetic-reference.json')
    write_document(fixture, output/'synthetic-reference.csv')
    roundtrip = compare_documents(read_document(output/'synthetic-reference.json'),
                                  read_document(output/'synthetic-reference.csv'))
    dump(roundtrip, output/'synthetic-roundtrip.json')
    cases = diagnostics['cases']
    representative = next(c for c in cases if c['kind'] == 'kink' and c['t_m'] == .1
                          and c['theta_deg'] == 90 and c['e_m'] == [.05, .05])
    summary = {'accepted': None, 'tolerance': '未定義', 'b03_complete': False,
               'case_count': len(cases),
               'kink_count': sum(c['kind'] == 'kink' for c in cases),
               'step_count': sum(c['kind'] == 'step' for c in cases),
               'range_rejected_case_count': sum(c['range_rejected'] for c in cases),
               'max_closed_vector_error_m': max(c['closed_vector_error_m'] for c in cases),
               'max_closed_kink_norm_error_m': max(c['closed_norm_error_m'] for c in cases if c['kind'] == 'kink'),
               'max_closed_step_norm_error_m': max(c['closed_norm_error_m'] for c in cases if c['kind'] == 'step'),
               'representative_90deg': representative,
               'legacy_rejections_confirmed': all(c['same_rejection_confirmed'] for c in legacy['cases']),
               'synthetic_roundtrip_consistent': roundtrip['contract_consistent'],
               'synthetic_is_body_output': False,
               'all_candidate_choices_undecided': all(d['accepted'] is None for c in cases for d in c['candidates'].values())}
    for name in representative['candidates']:
        defined = [c['candidates'][name] for c in cases if c['candidates'][name]['defined']]
        summary.setdefault('candidate_totals', {})[name] = {
            'defined_count': len(defined), 'undefined_count': len(cases)-len(defined),
            'max_geometry_error_m': max(d['max_geometry_error_m'] for d in defined),
            'max_represented_rigid_residual': max(d['represented_six_mode_max_residual'] for d in defined),
            'max_represented_work_cross_product_residual': max(d['represented_work_cross_product_residual'] for d in defined),
            'original_rigid_subspace_dimensions': sorted(set(d['original_rigid_subspace_dimension'] for d in defined)),
            'accepted': None,
        }
    dump(summary, output/'summary.json')
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_directory')
    parser.add_argument('--repository-root', default=str(Path(__file__).resolve().parents[2]))
    args = parser.parse_args()
    result = run(args.output_directory, args.repository_root)
    print({k: v for k, v in result.items() if k not in ('representative_90deg', 'candidate_totals')})
    raise SystemExit(2)
