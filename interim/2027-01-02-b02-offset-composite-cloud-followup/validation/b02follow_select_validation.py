"""Keep the failed first attempt; select one complete successful outcome per ID."""
from pathlib import Path
import json,shutil,sys
sys.path.insert(0,'interim/2027-01-02-b02-offset-composite-cloud-followup')
from b02follow_common import write
from b02follow_validation import GROUPS,source_binding,collect
root=Path('work/b02follow-full');p=root/'python';rerun=Path('work/b02follow-rerun-core-0')
D=Path('interim/2027-01-02-b02-offset-composite-cloud-followup')
initial=json.loads((p/'summary.json').read_text());assert initial['errors']==1 and initial['failures']==0 and initial['executed']==920
rerun_proof=json.loads((rerun/'source-proof.json').read_text());assert rerun_proof['sources_unchanged'] and rerun_proof['process']['exit_code']==0
bindings=[json.loads((folder/name).read_text()) for folder in (p,rerun) for name in ('source-binding-before.json','source-binding-after.json')]
assert all(v==bindings[0] for v in bindings) and source_binding()==bindings[0]
retained=D/'validation/retained-first-python';retained.mkdir(parents=True,exist_ok=False)
for name in ('summary.json','source-binding-before.json','source-binding-after.json','python-core-0-process.json','python-core-0.log'):
 shutil.copy2(p/name,retained/name)
shutil.copytree(p/'python-core-0',retained/'python-core-0')
# Only replace the new working selection, after preserving the failed attempt.
shutil.rmtree(p/'python-core-0');shutil.copytree(rerun/'python-core-0',p/'python-core-0')
for name in ('python-core-0-process.json','python-core-0.log'):
 shutil.copy2(rerun/name,p/name)
shutil.copy2(rerun/'source-proof.json',p/'core-0-rerun-source-proof.json')
reports=[json.loads((p/g/'python-test-inventory.json').read_text()) for g in GROUPS]
processes=[json.loads((p/(g+'-process.json')).read_text()) for g in GROUPS]
ids=[t for r in reports for t in r['test_elapsed_seconds']]
discovered=sorted(set(t for r in reports for t in r['discovered_test_ids']))
assigned=[t for r in reports for t in r['all_group_test_ids'][r['group']]]
assert sorted(ids)==sorted(assigned)==discovered and len(ids)==len(set(ids))==920
keys=('failures','errors','skipped','expected_failures','unexpected_successes')
assert all(not any(r[k] for k in keys) for r in reports) and all(v['exit_code']==0 for v in processes)
report=dict(processes=processes,discovered=len(discovered),executed=len(ids),executed_ids=sorted(ids),discovered_ids=discovered,assigned_ids=sorted(assigned),
 passed=sum(r['tests_run'] for r in reports),**{k:sum(r[k] for r in reports) for k in keys},
 unexecuted=len(set(discovered)-set(ids)),duplicate_ids=len(ids)-len(set(ids)),all_discovered_assigned_executed_exactly_once=True,
 source_count=len(bindings[0]),sources_unchanged=True,changed_sources=[],all_processes_exit_zero=True)
write(p/'summary.json',report)
collect(root.resolve())
v=json.loads((root/'validation.json').read_text());assert v['passed'] and v['python']['passed']==920
rerun_inventory=json.loads((rerun/'python-core-0/python-test-inventory.json').read_text())
v['initial_python_full_attempt']=dict(executed=initial['executed'],passed=initial['passed'],failures=initial['failures'],errors=initial['errors'],skipped=initial['skipped'],unexecuted=initial['unexecuted'],retained_at='validation/retained-first-python')
v['python_full_partition_execution_attempts']=initial['executed']+rerun_inventory['tests_run']
v['python_full_partition_attempt_outcomes']=dict(passed=initial['passed']+rerun_inventory['tests_run'],errors=initial['errors'],failures=initial['failures'])
v['selected_results']='One successful result for each of all 920 discovered IDs. The retained initial failed core-0 is replaced only in the selection by the full unchanged-source core-0 rerun; no duplicate counting.'
v['recovery']='P6-03 fixed BASE commit/tree absent under depth=1; explicit fixed filtered retrieval, 30 paths/29 unique frozen blobs verified; no code/reference/tolerance change.'
v['standalone_recovery_checks_not_added_to_full_counts']=dict(p603_restored=19,new_contracts_before_full=7)
v['archive_reproduction']='26 parts / 4067 entries verified; measurements and 127 derived files byte/SHA identical; no extra equilibrium solve'
v['source_commit_at_full_python']='84e87e03cf094e0e48468dad8271bdd84e4b7bec'
v['scientific_commit']='c4d63c56cd9caa0a0b78b20689c522afaa6041e8'
v['late_saved_orchestration_separate_from_tested_sources']=True
write(root/'validation.json',v)
write(D/'validation.json',v)
print(json.dumps(dict(passed=v['passed'],rust=v['rust'],python=v['python'],attempts=v['python_full_partition_execution_attempts']),indent=2))
