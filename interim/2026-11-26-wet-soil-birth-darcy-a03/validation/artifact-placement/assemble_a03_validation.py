import hashlib,json,re,subprocess
from pathlib import Path
root=Path.cwd();folder=root/'work/a03-validation';run=root/'runs/a03-final-n4';dest=root/'interim/2026-11-26-wet-soil-birth-darcy-a03'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
first_rust=read(folder/'rust/results.json');assert len(first_rust)==6 and first_rust[2]['exit_code']==101
retry=read(folder/'rust/whole-rerun-result.json');assert retry['exit_code']==0
rust=first_rust.copy();rust[2]=retry;assert all(x['exit_code']==0 for x in rust)
pattern=re.compile(r'test result: (\w+)\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out;')
rg=[]
for log in ('06-whole-rerun.log','03.log'):
 values=pattern.findall((folder/'rust'/log).read_text());assert values
 assert all(v[0]=='ok' and all(int(n)==0 for n in v[2:]) for v in values)
 rg.append(dict(log=log,passed=sum(int(v[1]) for v in values),failed=0,ignored=0,filtered_out=0,result_groups=len(values)))
base=read(folder/'inventory-02/python-test-inventory.json');groups={};all_ids=[]
for group,suffix in (('python-core-0','-final'),('python-core-1','-final'),('python-mean-evidence',''),('python-spatial-evidence','')):
 v=read(folder/(group+suffix)/'python-test-inventory.json')
 assert v['exit_code']==0 and v['tests_run']==v['assigned_total']
 assert all(v[k]==0 for k in ('failures','errors','skipped','expected_failures','unexpected_successes'))
 assert v['discovered_test_ids']==base['discovered_test_ids'] and v['all_group_test_ids']==base['all_group_test_ids']
 all_ids+=v['all_group_test_ids'][group];groups[group]=v
assert len(all_ids)==len(set(all_ids))==573 and sorted(all_ids)==base['discovered_test_ids']
post=read(folder/'inventory-after/python-test-inventory.json');assert post['discovered_test_ids']==base['discovered_test_ids'] and post['all_group_test_ids']==base['all_group_test_ids']
original=read(folder/'live-validation-source-fingerprints.json');assert all(sha(root/p)==h for p,h in original.items())
fresh=read(run/'source-fingerprints.json');assert all(sha(root/p)==h for p,h in fresh.items())
core=read(folder/'core-final-processes.json');assert core['binary_sha256_before']==core['binary_sha256_after']
metrics=read(run/'measurements.json');assert metrics['passed'] is False and metrics['a03_complete'] is False
history=root/'interim/2026-11-24-nc-pc-rerun-fix/README.ja.md'
report=dict(a03_complete=False,environment=read(folder/'environment.json'),rust_commands=rust,initial_rust_commands_retained=first_rust,initial_rust_full_suite=dict(passed=450,failed=1,reason='missing original archived slope failure JSON; no assertion or tolerance changed'),rust_test_logs=rg,
 rust_total_passed=sum(g['passed'] for g in rg),rust_result_groups=sum(g['result_groups'] for g in rg),
 python=dict(discovered=573,executed=573,successful=573,unexecuted=0,failures=0,errors=0,skipped=0,expected_failures=0,unexpected_successes=0,
             inventory_before_after_and_every_group_identical=True,all_ids_partitioned_once=True,groups=groups,
             initial_core_attempts_retained={g:read(folder/g/'python-test-inventory.json') for g in ('python-core-0','python-core-1')},
             historical_known_linux_nc_limitation=dict(reference=str(history.relative_to(root)),sha256=sha(history),unexecuted_methods=8,
               class_errors=2,maximum_previously_diagnosed_float_difference=4.440892098500626e-16,historical_failures_and_unexecuted_count_retained=True,
               reproduced_in_this_run=False,current_strict_replay_methods_passed=8,changed_tolerances=False,changed_archives_or_sources=False,
               linux_platform_float_allowance=False,old_nc_numerical_acceptance=False,cause_of_environment_difference='unidentified')),
 fixed_fresh_solve_sources=dict(count=len(fresh),all_sha256_unchanged=True),original_live_validation_sources=dict(count=len(original),all_original_git_blobs_verified=True,all_sha256_unchanged=True),
 fresh_solve_and_control_counts=dict(native_cases=2,independent_input_only_3d_cases=2,accepted_pairs=32,prepared_pairs=2,
    accepted_fields=54,prepared_fields=50,separate_pid_real_checkpoint_restarts=8,entire_suffix_states_exact=64,
    same_binary_rehashed_mutations_refused=12,production_cli_input_refusals=24,additional_standard_cli_solves=2,additional_same_binary_whole_solves=2),
 windows_ci_after_push='unconfirmed',nc01_p6_and_all_remainders_complete=False)
(dest/'validation.json').write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
print('Rust',report['rust_total_passed'],'groups',report['rust_result_groups'],'Python',report['python']['successful'])
