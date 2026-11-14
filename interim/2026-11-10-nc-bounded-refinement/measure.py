"""Audit all saved 3D fields and explicitly limited actual-matrix evidence."""
import argparse
from collections import Counter, deque
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from linear_evidence_codec import read_rows

spec = importlib.util.spec_from_file_location('previous_nc_measure', ROOT/'interim/2026-11-10-nc-mean-stress-release/measure.py')
previous = importlib.util.module_from_spec(spec); spec.loader.exec_module(previous)


def read(path):
    return json.loads(path.read_text())


def raw_rows(path):
    with gzip.open(path, 'rt') as stream:
        for line in stream:
            yield json.loads(line)


def audit(folder):
    original = previous.rows
    def fields_only(path):
        # The previous field/force/history/domain gates remain unchanged.
        # Its matrix loop only accumulates statistics; audit summaries and
        # the retained actual matrices separately, without inventing matrices.
        return iter(()) if path.name == 'linear-systems.jsonl' else original(path)
    previous.rows = fields_only
    try:
        result = previous.audit(folder)
    finally:
        previous.rows = original
    result.pop('maximum_remeasured_linear_residual_including_rejections')
    controls = read(folder/'controls.json')
    if controls['binding']['source'] != read(folder/'restart.json')['binding']['source']:
        raise ValueError('source/checkpoint binding differs')
    if not result['force_decomposition_bounded']:
        raise ValueError('auxiliary load cannot enter this experiment')
    manifest = read(folder/'streaming-digests.json')
    for name, record in manifest.items():
        if name == 'linear-systems.jsonl':
            continue
        checksum = hashlib.sha256(); size = 0; count = 0
        with gzip.open(folder/(name+'.gz'), 'rb') as stream:
            for line in stream:
                checksum.update(line); size += len(line); count += 1
        if checksum.hexdigest() != record['raw_sha256'] or size != record['raw_size'] or count != record['rows']:
            raise ValueError('complete saved stream changed: '+name)
    tail = deque(maxlen=128); count = 0; maximum = 0.; canonical_maximum = 0.; layout_difference = 0.; layout_difference_count = 0
    for row in raw_rows(folder/'linear-systems.jsonl.gz'):
        if not np.isfinite(row['rhs']).all() or not np.isfinite(row['solution']).all():
            raise ValueError('nonfinite linear summary')
        # The actual solve may use a strided/F-order array. Rebuilding its
        # values as contiguous arrays can change BLAS summation roundoff.
        # Keep both evaluations and their difference; never replace the gate
        # used by the solver or pretend these floats were bit-identical.
        canonical = row['independently_recomputed_true_relative_residual']
        if not np.isfinite(canonical):raise ValueError('nonfinite canonical residual')
        difference = abs(row['true_relative_residual']-canonical)
        layout_difference = max(layout_difference, difference)
        layout_difference_count += difference != 0.
        canonical_maximum = max(canonical_maximum, canonical)
        maximum = max(maximum, row['true_relative_residual']); count += 1; tail.append(row)
    if count != manifest['linear-systems.jsonl']['rows']:
        raise ValueError('linear summary count differs')
    actual = list(read_rows(folder/'terminal-linear-systems.jsonl.f64.gz'))
    if len(actual) != len(tail):
        raise ValueError('terminal actual-system count differs')
    maximum_actual = 0.; worst = None
    for row, summary in zip(actual, tail):
        raw = (json.dumps(row, allow_nan=False)+'\n').encode()
        A = np.array(row['matrix']); rhs = np.array(row['rhs']); solution = np.array(row['solution'])
        if hashlib.sha256(raw).hexdigest() != summary['raw_record_sha256'] or hashlib.sha256(np.asarray(A, dtype='<f8').tobytes()).hexdigest() != summary['matrix_sha256']:
            raise ValueError('terminal actual system differs from recorded stream')
        residual = float(np.linalg.norm(A@solution-rhs)/max(np.linalg.norm(rhs), 1e-30))
        if residual != summary['independently_recomputed_true_relative_residual']:
            raise ValueError('terminal actual linear residual differs')
        maximum_actual = max(maximum_actual, residual)
        left, singular, right = np.linalg.svd(A, full_matrices=False)
        condition = float(singular[0]/singular[-1]) if singular[-1] else None
        witness = dict(label=row['label'], iteration=row['iteration'], matrix_sha256=summary['matrix_sha256'],
            smallest_singular_value=float(singular[-1]), largest_singular_value=float(singular[0]), condition_number_2=condition,
            rhs_projection_on_smallest_left_mode=float(abs(left[:, -1]@rhs)/max(np.linalg.norm(rhs), 1e-30)),
            correction_norm_m=float(np.linalg.norm(solution)), true_linear_residual=residual,
            nonsymmetric_matrix_not_symmetrized=True, physical_static_limit_certified=False)
        if worst is None or (condition is None or (worst['condition_number_2'] is not None and condition > worst['condition_number_2'])):
            worst = witness
    refinement = list(raw_rows(folder/'linear-refinement.jsonl.gz'))
    improved = [x for x in refinement if x['returned_relative_residual'] < x['original_relative_residual']]
    rescued = [x for x in refinement if x['original_relative_residual'] > 1e-10 and x['returned_relative_residual'] <= 1e-10]
    attempts = list(original(folder/'attempts.jsonl'))
    result.update(linear_systems=count, all_linear_matrices_archived=False, actual_terminal_matrices=len(actual),
        maximum_recorded_true_linear_residual_including_rejections=maximum,
        maximum_runtime_canonical_linear_residual_including_rejections=canonical_maximum,
        original_vs_canonical_residual_maximum_absolute_difference=layout_difference,
        original_vs_canonical_residual_difference_count=layout_difference_count,
        maximum_remeasured_terminal_true_linear_residual_including_rejections=maximum_actual,
        complete_raw_linear_stream_sha256_not_reconstructible=manifest['linear-systems.jsonl']['raw_sha256'],
        worst_retained_actual_matrix=worst, refinement_calls=len(refinement), improved_solutions=len(improved),
        linear_gate_rescued_without_relaxation=len(rescued),
        rejection_reasons=dict(Counter(x['stop_reason'] for x in attempts if not x['accepted'])),
        final_rejected_attempt=next((x for x in reversed(attempts) if not x['accepted']), None))
    return result


def measure(evidence):
    paths = sorted(p for p in evidence.iterdir() if (p/'controls.json').exists() and (p/'report.json').exists())
    result = dict(passed=False, full_removal_accepted=False, history_converged=False, static_limit_confirmed=False,
        production_adopted=False, native_committed_increments=0, paths={}, step_sensitivity={}, linear_policy_sensitivity={},
        material_sensitivity={}, prehistory_sensitivity={}, cross_map={}, restart={}, full_release_freeze_verified=False,
        auxiliary_load_paths_accepted=False, starting_commit='5bee3b47f61878e614cae1e7268a57449d2f12d1')
    for path in paths:
        result['paths'][path.name] = audit(path)
        print('audited '+path.name, file=sys.stderr, flush=True)
    for kind in ['native', 'independent']:
        base = f'{kind}-complement-refined128'
        result['step_sensitivity'][kind] = previous.compare(evidence/base, evidence/f'{kind}-complement-refined256')
        result['linear_policy_sensitivity'][kind] = previous.compare(evidence/f'{kind}-complement-plain128', evidence/base)
    base = evidence/'independent-complement-refined128'
    result['material_sensitivity'] = previous.compare(base, evidence/'independent-complement-refined128-m8')
    result['prehistory_sensitivity'] = previous.compare(base, evidence/'independent-complement-refined128-pre800')
    for order in ['complement', 'mean']:
        result['cross_map'][order] = previous.compare(evidence/f'native-{order}-refined128', evidence/f'independent-{order}-refined128')
    for path in paths:
        if path.name.endswith('-resumed'):
            direct = path.name.removesuffix('-resumed')
            record = previous.restart_compare(evidence, direct, path.name)
            # SHA of every full original system plus rhs/solution is compared.
            # Only terminal actual matrices are included in this new archive.
            record['actual_linear_systems_exact'] = False
            record['all_linear_system_record_hashes_rhs_solution_exact'] = True
            record['compared_system_record_hashes_rhs_solutions'] = record.pop('compared_actual_linear_systems')
            record['actual_terminal_matrices_archived'] = result['paths'][path.name]['actual_terminal_matrices']
            terminal_direct = (evidence/direct/'terminal-linear-systems.jsonl.f64.gz').read_bytes()
            terminal_resumed = (path/'terminal-linear-systems.jsonl.f64.gz').read_bytes()
            if terminal_direct != terminal_resumed:
                raise ValueError('restart retained actual terminal matrices differ')
            record['retained_actual_terminal_systems_byte_identical'] = True
            result['restart'][direct] = record
    result['independent_prehistories'] = previous.prehistory(evidence)
    prefix = read(evidence/'fresh-native-prefix/verification.json')
    fresh = read(evidence/'fresh-native-prefix/result.json')
    original = read(evidence/'prepared/origin/prior/fresh-head-midpoint/job/recovery-checkpoint-0.json')['run']
    if next(x['state'] for x in fresh['frames'] if x['name']=='unload') != next(x['state'] for x in original['frames'] if x['name']=='unload'):
        raise ValueError('current production prefix differs from sealed native history')
    result['fresh_native_prefix'] = prefix
    result['processes'] = read(evidence/'processes.json')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, required=True); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(measure(args.evidence.resolve()), indent=2, allow_nan=False)+'\n')
