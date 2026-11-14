"""Aggregate new A-04 observations, preserving every scope limit and failure."""
import argparse
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
from study import read,write,sha


def compact(comparison):
    return {k:v for k,v in comparison.items() if k!='all_points'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--work',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();work=a.work.resolve();out=a.output.resolve()
    soil=read(work/'formal-soil-03/measurements.json');audits=read(work/'formal-soil-03/supplemental-audit.json');convergence=read(work/'formal-soil-03/common-grid-convergence.json')
    dry=read(work/'formal-dry-01/measurements.json');wet=read(work/'formal-hydrostatic-01/measurements.json');refusals=read(work/'formal-refusals-01/refusals.json')
    cases={k:dict(comparison=compact(v['comparison']),observation=v['observation'],existing_bounded_gates_passed=v['existing_bounded_gates_passed'],restarts=v.get('restarts',[]),tampering=v.get('tampering',[])) for k,v in soil['cases'].items()}
    stages=soil['cases'];checks=[r for v in stages.values() for r in v.get('restarts',[])]+dry['restarts']+wet['restarts']
    mutations=[r for v in stages.values() for r in v.get('tampering',[])]+dry['tampering']+wet['tampering']+[r for d in (dry,wet) for r in d['state_spoofing'] if r['outer_hash_recomputed']]
    undefined=[r for d in (dry,wet) for r in d['state_spoofing'] if not r['outer_hash_recomputed']]
    output=dict(passed=False,a04_complete=False,base=soil['base'],label_date='2026-11-29',execution_date_jst='2026-11-12',
        production_physics_or_controller_changed=False,bidirectional_implemented=False,
        scope_judgments={
            'soil_only_frozen_porosity_ONE_WAY_CAD_DEATH_MC_elastic':'限定受入の維持・再確認（既存の固定格子、小型2ケースと定義済みゲートだけ）',
            'dry_single_tip_B_death_rebirth':'限定受入の維持・再確認（A-02の既存範囲）',
            'fully_submerged_hydrostatic_fresh_MC_CONTINUOUS_B':'限定受入の維持・再確認（A-03の既存範囲）',
            'Richards_ONE_WAY_CONTINUOUS':'未達／未対応（現行拒否を保持）',
            'Richards_ONE_WAY_SOIL':'未達／未対応（現行拒否を保持）',
            'Richards_ONE_WAY_PRESCRIBED':'未達／未対応（現行拒否を保持）',
            'Richards_bidirectional':'未達／未対応（controller/schema/CLIなし、実装しない）',
            'Richards_interface_head_mapping':'未達／未対応（現行拒否を保持。節点投影や静水位で補完しない）',
            'Richards_joint_open_close_generations_in_physical_seconds':'未達／未対応（準静的JOINT世代と混同しない）',
            'Richards_joint_transient_mass_budget':'未達／未対応（土の台帳を転用しない）',
            'Richards_joint_hydraulic_mechanical_load_transfer':'未達／未対応（土GP伝達と静水圧面荷重を連成としない）',
            'Richards_joint_coupled_time_refinement':'未達／未対応（土だけの物理秒細分化を転用しない）'},
        fresh_native_paths=8,fresh_input_only_independent_paths=8,
        accepted_pairs=sum(v['comparison']['accepted_pairs'] for v in cases.values())+dry['comparison']['accepted_pairs']+wet['comparison']['accepted_pairs'],
        soil_cases=cases,soil_audits=audits,
        dry_control={**dry,'comparison':compact(dry['comparison'])},hydrostatic_control={**wet,'comparison':compact(wet['comparison'])},
        convergence={k:{**v,'pairs':[compact(c) for c in v['pairs']]} for k,v in convergence.items()},
        checkpoint_restart=dict(real_files=len(checks),suffix_states_exact=sum(r['suffix_states_exact'] for r in checks),whole_run_exact_except_elapsed=True,
            same_producer_binary=True,all_input_and_binary_hashes_compared=True),
        tamper_refusal=dict(authenticated_rehashed_checkpoint_cases=len(mutations),unknown_physical_time_field_cases=len(undefined),same_producer_binary=True,
            typed_unknown_field_refusal_is_not_replay_refusal=True),
        actual_input_cli_refusals=refusals,
        boundary_consistency_gap=read(work/'boundary-consistency-gap.json'),
        retained_failures=['initial dense solver dimension budget rejected','saturated preparation without declared cut reservoir cancelled/incomplete',
            'original derivative width first-cut nodal force comparison exceeded existing gate','supplemental resultant diagnostic initially misused nodal criterion; retained and corrected',
            'supplemental support reconstruction initially included free residual; retained and corrected to constrained projection',
            'earlier saturated series used mismatched exterior reservoir and hydraulic Dirichlet heads; excluded from limited acceptance, production refusal gap retained'],
        source_count= len(read(work/'formal-soil-03/source-fingerprints.json')),fresh_solve_sources_unchanged=soil['sources_unchanged'],
        numerical_acceptance=False,note='No all-field/global-resultant/JOINT transient threshold is invented. Saved-state audits are not fresh solves. P6/NC and all outstanding rows remain incomplete.')
    write(out/'measurements.json',output)
    write(out/'soil-pointwise-audit.json',audits);write(out/'common-grid-convergence.json',convergence);write(out/'scope-refusals.json',refusals)
    print('accepted pairs',output['accepted_pairs'],'checkpoints',len(checks),'suffix',output['checkpoint_restart']['suffix_states_exact'],'authenticated mutations',len(mutations),flush=True)
if __name__=='__main__':main()
