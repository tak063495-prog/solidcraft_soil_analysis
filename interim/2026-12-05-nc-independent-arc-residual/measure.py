"""Measure actual original bordered systems; no solve of a physical path here."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys

import numpy as np
import scipy
from scipy.linalg import lu_factor

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(ROOT/'tools'), str(ROOT/'interim/2026-11-09-nc-pseudo-arclength')]
from independent_arc import solve_bordered
from diagnose_nc_dense_lu import decimal_residual


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def matrix_metrics(row):
    a, b, x = (np.array(row[key]) for key in ('matrix', 'rhs', 'solution'))
    residual = a@x-b
    rows = np.max(abs(a), axis=1); columns = np.max(abs(a), axis=0)
    scaled = a / rows[:, None]; scaled /= np.max(abs(scaled), axis=0)
    factor, pivots = lu_factor(a)
    y, diagnostics = solve_bordered(a, b)
    return dict(dimension=len(b), condition_number_2=float(np.linalg.cond(a)),
        stiffness_condition_number_2=float(np.linalg.cond(a[:-1, :-1])),
        equilibrated_condition_number_2=float(np.linalg.cond(scaled)),
        row_max_abs_range=[float(rows.min()), float(rows.max())],
        column_max_abs_range=[float(columns.min()), float(columns.max())],
        displacement_column_max_abs_range=[float(columns[:-1].min()), float(columns[:-1].max())],
        lambda_column_max_abs=float(columns[-1]), arc_row_max_abs=float(rows[-1]),
        direct_pivot_swaps=int(np.sum(pivots != np.arange(len(pivots)))),
        direct_abs_u_pivot_range=[float(abs(np.diag(factor)).min()), float(abs(np.diag(factor)).max())],
        direct_u_growth=float(np.max(abs(np.triu(factor)))/np.max(abs(a))),
        rhs_norm_2=float(np.linalg.norm(b)), solution_norm_2=float(np.linalg.norm(x)),
        original_float64_relative_residual=float(np.linalg.norm(residual)/max(np.linalg.norm(b), 1e-30)),
        original_decimal80_relative_residual=decimal_residual(a.tolist(), b.tolist(), x.tolist()),
        original_force_block_absolute_residual_2=float(np.linalg.norm(residual[:-1])),
        original_arc_row_absolute_residual=float(abs(residual[-1])),
        relative_backward_error=float(np.linalg.norm(residual)/(np.linalg.norm(a)*np.linalg.norm(x)+np.linalg.norm(b))),
        improved_float64_relative_residual=float(np.linalg.norm(a@y-b)/max(np.linalg.norm(b), 1e-30)),
        improved_decimal80_relative_residual=decimal_residual(a.tolist(), b.tolist(), y.tolist()),
        solver=diagnostics)


def run(work):
    result = dict(source='Linux cloud new computation, not PC transcription', pc_data=False,
        base_commit='e278dbf17c771e87098b6180e9c65d9e4a6ad2e7',
        implementation_commit='10b81658985ba6bbd8ddecb647e991859ed99d52',
        platform=platform.platform(), python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
        mesh=dict(nodes=54, tet10=12, note='128/400 below are independent history integration step counts, not mesh sizes'),
        acceptance_norm='float64 ||A x-b||_2 / max(||b||_2,1e-30), original unscaled A,b and actual pre-normalized x',
        tolerance=1e-10, passed=False, nc01_complete=False, references={}, traces={})
    for n in (128, 400):
        source=read(work/f'reference{n}/reference.json')
        result['references'][str(n)]={key:source[key] for key in ('steps','integration_method','elapsed_seconds','failure_reason','maximum_relative_free_force_residual','last_system_maximum_true_linear_residual')}
        result['references'][str(n)].update(last_accepted_stage=source['last_accepted']['name'],last_accepted_progress=source['last_accepted']['progress'],new_cloud_history=True)
    for name in ('before128','before400','after400'):
        report=read(work/name/'report.json'); sample_path=work/name/'samples.json'
        samples=read(sample_path) if sample_path.is_file() else []
        if len(samples)!=report['equilibrium_samples']:
            raise ValueError('reported sample count differs from saved trace: '+name)
        result['traces'][name]=dict(report=report, samples=len(samples), peak_sampled_lambda=max((x['progress'] for x in samples),default=None),
            step_size_m=2e-5, requested_count=(16 if name=='before128' else 8),
            fold_confirmed=False, limitations='short bounded rate trace, no equilibrium corrector or full history acceptance')
    for name in ('before128','before400'):
        summary=read(work/f'{name}-systems-summary.json')
        result['traces'][name].update(actual_bordered_solves=summary['actual_bordered_solves'],worst_actual_system=matrix_metrics(summary['worst']))
    rows=read(work/'after400/linear-systems.json')
    measured=[];max_difference=0.
    for row in rows:
        a,b,x=(np.array(row[k]) for k in ('matrix','rhs','solution'))
        value=float(np.linalg.norm(a@x-b)/max(np.linalg.norm(b),1e-30));measured.append(value)
        max_difference=max(max_difference,abs(value-row['true_relative_residual']))
    result['traces']['after400'].update(actual_bordered_solves=len(rows),maximum_remeasured_original_residual=max(measured,default=None),
        saved_reported_residual_maximum_difference=max_difference,
        all_saved_actual_systems_within_unchanged_gate=bool(result['traces']['after400']['samples'] and rows) and all(v<=1e-10 for v in measured),
        maximum_unscaled_direct_residual_on_new_path=max((r['solver']['unscaled_direct_relative_residual'] for r in rows),default=None),
        worst_actual_system=(matrix_metrics(rows[int(np.argmax(measured))]) if rows else None))
    result['pc_failure_reproduced_on_linux']=False
    result['limitations']=['PC-2 actual matrix is unavailable; its exact pivot/error cannot be identified',
        'Linux same steps400 short interval did not exceed 1e-10 even before correction',
        'history400/800 full arc128, material4/8 and checkpoint comparisons remain PC work',
        'reference histories stop on the old stiffness linear gate, which is not modified']
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--work',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.write_text(json.dumps(run(args.work),ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
