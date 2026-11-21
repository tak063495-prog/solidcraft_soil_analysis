import json,hashlib,shutil,sys
from pathlib import Path
sys.path.insert(0,'interim/2026-12-19-c04-transient-water-cause-diagnosis')
from c04cause_driver import write
study=Path('interim/2026-12-19-c04-transient-water-cause-diagnosis');ev=Path('work/c04-evidence')
obs=Path('work/c04-observation-final');full=json.loads((obs/'measurements.json').read_text())
shutil.copytree(obs,ev/'data/observation-final')
for name in ('diagnostic-fixture.json','limitations-and-failure-record.json'):
 shutil.copy2(study/name,ev/'reports'/name)
for path in Path('work/c04-validation-preparation').iterdir():
 if path.is_file():shutil.copy2(path,ev/'retained/preparation'/path.name)
main=dict(schema='c04-transient-water-cause-summary/1',passed=False,c04_complete=False,a04_complete=False,
 c01_complete=False,c02_complete=False,c03_complete=False,nc01_complete=False,p6_complete=False,ab_items_completed=False,
 pc_executed=False,original_acceptance_all_met=False,cloud_diagnostics_remaining=True,unique_cause_fully_identified=False,
 status='未達; limited observations are not limited acceptance',
 original_l2_definition=full['reduction_order']+'; ||a-b||2/max(||b||2,1e-12)',
 full_all_fields_and_stages_in_sealed_archive='data/observation-final/measurements.json',levels={},global_controls={},
 global_contributions=json.loads(Path('work/c04-global-contributions.json').read_text()),
 formal_solve_paths=12,formal_accepted_states=7008,preliminary_solve_paths_retained=10,preliminary_accepted_states=4832,
 total_scientific_completed_solve_calls=22,total_scientific_accepted_states=11840,
 original_saved_values_rewritten=False,richards_wet_joint_connected=False,current_richards_refusal_verified=True,
 old_hydraulic_rt0_mass_balance_retained=True,hydraulic_to_joint_map=None,
 old_sealed_replay=json.loads(Path('work/c04-old-sealed-replay/verification.json').read_text()),
 old_baseline_binding=json.loads(Path('work/c04-global-final/old-baseline-binding.json').read_text()))
for k,v in full['levels'].items():
 audit=v['native_global_audit'];g=audit['gp_correspondence']
 main['levels'][k]=dict(original_all_fields=v['comparison']['all_fields'],original_endpoints=v['comparison']['endpoints'],
  mismatch_count=v['mismatch_count'],margin_bins=v['margin_bins'],zero_saved_plastic_increment_both=v['zero_saved_plastic_increment_both'],
  native_parts=v['native_parts'],independent_parts=v['independent_parts'],dissipation_identity=v['dissipation_identity'],
  stage_dissipation_increment_difference_l2=v['stage_dissipation_increment_difference_l2'],common_saved_jump_replay=v['common_saved_jump_replay'],
  native_global_audit={x:y for x,y in audit.items() if x!='gp_correspondence'},independent_global_audit=v['independent_global_audit'],
  gp_correspondence={x:y for x,y in g.items() if not isinstance(y,list)})
quantities=('u','stress','strain','traction','plastic_slip','accumulated','plastic_dissipation','numerical_dissipation','endpoint_work','contact_status','open')
for k,v in full['global_controls'].items():
 main['global_controls'][k]=dict(accepted_pairs=v['comparison']['accepted_pairs'],all_fields={q:v['comparison']['all_fields'][q] for q in quantities},
  endpoints={q:v['comparison']['endpoints'][q] for q in quantities if q in v['comparison']['endpoints']},
  global_audit={q:w for q,w in v['global_audit'].items() if q!='gp_correspondence'})
write(study/'measurements.json',main)
shutil.copy2('work/c04-reference-assessment.json',study/'reference-assessment.json')
shutil.copy2('work/c04-global-final/old-baseline-binding.json',study/'old-baseline-binding.json')
# Check every formal path, not only the representative pair.
actual=[]
for path in sorted(Path('work/c04-global-final').iterdir()):
 if not (path/'fields.jsonl').is_file():continue
 count=sum(1 for _ in (path/'fields.jsonl').open())
 if path.name not in full['global_controls']:raise ValueError('formal audit missing '+path.name)
 actual.append(dict(name=path.name,accepted_states=count,complete=True))
if len(actual)!=12 or sum(x['accepted_states'] for x in actual)!=7008:raise ValueError(actual)
write(study/'solve-call-ledger.json',dict(formal_paths=actual,formal_path_count=12,formal_states=7008,
 preliminary_paths=10,preliminary_states=4832,total_completed_scientific_calls=22,total_accepted_states=11840,
 refused_cli_invocations=3,software_test_solves_counted=False,archive_replays_counted_as_solves=False))
shutil.copy2(study/'solve-call-ledger.json',ev/'reports/solve-call-ledger.json')
shutil.copy2(study/'measurements.json',ev/'reports/summary.json')
print('science finalized')
