"""Reproduce observations from sealed new reassemblies, using Decimal only."""
import argparse
from decimal import Decimal, localcontext
import hashlib
import json
from pathlib import Path
from a01hp_common import read, write


def norm(values):
    with localcontext() as c:
        c.prec=90
        return sum((Decimal(str(x))**2 for x in values),Decimal(0)).sqrt()


def flatten(matrix):
    return [v for row in matrix for v in row]


def difference(a,b):
    return norm([Decimal(x)-Decimal(y) for x,y in zip(flatten(a),flatten(b))])


def collect(raw,output):
    low,high,binary,frozen=(read(raw/(n+'.json')) for n in ('rebuilt50','rebuilt80','binary-rebuilt80','frozen80'))
    binary_low=read(raw/'binary-rebuilt50.json')
    double=read(raw/'double-with-both-svd.json')
    names=('analytic','fd_h_div_1','fd_h_div_2','fd_h_div_4')
    with localcontext() as c:
        c.prec=90
        matrices={name:dict(
            spectrum_50=low['records'][name]['spectrum_descending'],
            spectrum_80=high['records'][name]['spectrum_descending'],
            binary_input_spectrum_50=binary_low['records'][name]['spectrum_descending'],
            binary_input_spectrum_80=binary['records'][name]['spectrum_descending'],
            rank_50=low['records'][name]['numerical_rank'],rank_80=high['records'][name]['numerical_rank'],
            cutoff_50=low['records'][name]['observation_cutoff'],cutoff_80=high['records'][name]['observation_cutoff'],
            high_precision_matrix_change_50_to_80=str(difference(low['raw'][name]['total_matrix'],high['raw'][name]['total_matrix'])),
            input_binary64_matrix_change_at_80=str(difference(high['raw'][name]['total_matrix'],binary['raw'][name]['total_matrix'])),
            binary_input_matrix_change_50_to_80=str(difference(binary_low['raw'][name]['total_matrix'],binary['raw'][name]['total_matrix'])),
            frozen_layout_target_matrix_change_at_80=str(difference(high['raw'][name]['total_matrix'],frozen['raw'][name]['total_matrix'])),
            differences_80=high['records'][name]['differences'],
            right_rigid_overlap_squared=[r['middle_rigid_overlap_squared'] for r in high['records'][name]['two_smallest']],
            contact_rigid_spectrum=high['records'][name]['contact_rigid_spectrum'],
            # A sign convention fixes each individual vector. For the two
            # degenerate analytic directions only the projector is comparable.
            left_projector_50_to_80=str(difference(low['records'][name]['left_projector_two'],high['records'][name]['left_projector_two'])),
            right_projector_50_to_80=str(difference(low['records'][name]['right_projector_two'],high['records'][name]['right_projector_two'])),
            signed_projections=[r['projections'] for r in high['records'][name]['two_smallest']],
            projection_norms=high['records'][name]['projections_two'],
            true_null_projection_norms=high['records'][name]['projections_true_numerical_null'],
            svd_reconstruction_error=high['records'][name]['svd_reconstruction_error']) for name in names}
        fd=[Decimal(high['records']['fd_h_div_'+str(d)]['spectrum_descending'][-2]) for d in (1,2,4)]
        fd_ratios=[str(fd[i]/fd[i+1]) for i in range(2)]
        dps_projection_change={name:str(abs(Decimal(low['records'][name]['projections_two']['saved_residual']['left_norm'])-Decimal(high['records'][name]['projections_two']['saved_residual']['left_norm']))) for name in names}
        # Signed residual balance in the same left singular vectors.
        signed_balance=[]
        for mode in high['records']['analytic']['two_smallest']:
            p=mode['projections']
            signed_balance.append(str(Decimal(p['reassembled_residual']['left_component'])-(Decimal(p['external_load']['left_component'])-Decimal(p['internal_bulk']['left_component'])-Decimal(p['internal_contact']['left_component']))))
        probe_branch_changes=[]
        for divisor,gps in high['fd_branch_paths'].items():
            for q,probes in enumerate(gps):
                for probe in probes:
                    for sign in ('plus','minus'):
                        if probe[sign]['path']!=high['branches'][q]['path']:
                            probe_branch_changes.append(dict(divisor=divisor,gp=q,axis=probe['axis'],sign=sign))
        source_matrix=dict(
            double_spectrum_with_vectors=double['records']['saved_pre_predictor']['spectrum_with_vectors_descending'],
            double_spectrum_values_only=double['records']['saved_pre_predictor']['spectrum_descending'],
            promoted_exact_decimal_spectrum_80=high['records']['saved_pre_predictor']['spectrum_descending'],
            promoted_exact_binary64_spectrum_80=binary['records']['saved_pre_predictor']['spectrum_descending'],
            promoted_exact_binary64_spectrum_50=binary_low['records']['saved_pre_predictor']['spectrum_descending'],
            saved_solver_spectrum_binary64_promoted_80=binary['records']['saved_solver_matrix']['spectrum_descending'],
            saved_predictor_change=high['saved_predictor_matrix_change'])
        result=dict(passed=False,a01_complete=False,diagnostic_only=True,
            diagnostic_completed=True,new_limited_acceptance=False,element_count=18,reduced_dofs=48,
            source='12-28 sealed stopped native open stage, progress 0.125, iteration 8',
            exact_analytic_null_directions=2,exact_fd_null_directions=1,
            structural_reason='Middle body has six rigid DOFs. Only two distinct final slip GPs constrain four independent directions. Their final Coulomb tangents each have rank two. FD radial stiffness adds one independent constraint; rotation about the line of the two GP positions remains free.',
            matrices=matrices,fd_small_value_ratio_h_to_half=fd_ratios,
            fd_probe_branch_changes=probe_branch_changes,
            projection_50_to_80_changes=dps_projection_change,
            residual_reassembly_difference=high['residual_reassembly_difference'],
            residual_sign=high['residual_sign'],signed_analytic_projection_balance_errors=signed_balance,
            birth_holding_active_at_stop=False,source_matrix=source_matrix,
            frozen_linearized_rhs_is_in_range=False,
            roundoff_is_sole_stopping_mechanism=False,
            fd_artificial_stiffness_h_squared_observed=True,
            physical_instability_proved=False,physical_instability_classification_d='判定不能',
            new_global_solves=0,cloud_specific_remainder_completed=True,
            known_immediate_cloud_diagnostics_remaining=[],pc_required_due_to_resources=False,
            user_decisions_pending=['湿潤先端の水圧/作用面と初期力・保持力','複数先端・部材との拘束/投影/所有・世代・仕事会計','数値方針の採否と(d)を評価する平衡枝追跡・安定性の定義'],
            cutoff_is_acceptance_tolerance=False,physical_percentages_summing_to_100=False,
            raw_files_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(raw.glob('*.json'))})
    write(output,result)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--raw',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();collect(a.raw,a.output)
