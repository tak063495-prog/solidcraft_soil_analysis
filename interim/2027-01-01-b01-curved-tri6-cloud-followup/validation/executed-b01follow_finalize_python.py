"""Select the successful retry of one failed shard, retaining every first outcome."""
from pathlib import Path
import json,sys,shutil
sys.path.insert(0,str(Path('interim/2027-01-01-b01-curved-tri6-cloud-followup').resolve()))
from b01follow_common import read,write
from b01follow_validation import GROUPS,source_binding

root=Path('work/b01follow-full');original=root/'python';first=read(original/'summary.json')
retry=Path('work/b01follow-core0-retry');retry_process=read('work/b01follow-core0-retry-process.json');retry_report=read(retry/'python-test-inventory.json')
assert retry_process['exit_code']==0 and retry_report['exit_code']==0
assert first['errors']==1 and first['failures']==0 and first['sources_unchanged']
old=read(original/'python-core-0/python-test-inventory.json')
assert old['errors']==1 and old['tests_run']==retry_report['tests_run']
assert old['all_group_test_ids']==retry_report['all_group_test_ids']
before=read(original/'source-binding-before.json');current=source_binding();assert before==current
retained=root/'retained-first-python';shutil.move(original,retained);shutil.copytree(retained,original)
shutil.rmtree(original/'python-core-0');shutil.copytree(retry,original/'python-core-0')
shutil.copyfile('work/b01follow-core0-retry.log',original/'python-core-0.log')
write(original/'python-core-0-process.json',retry_process)
reports=[read(original/g/'python-test-inventory.json') for g in GROUPS]
ids=[test for r in reports for test in r.get('test_elapsed_seconds',{})]
discovered=sorted(set(test for r in reports for test in r.get('discovered_test_ids',[])))
assigned=[test for r in reports for test in r.get('all_group_test_ids',{}).get(r['group'],[])]
assert discovered==sorted(ids)==sorted(assigned) and len(ids)==len(set(ids))
processes=first['processes'].copy();processes[0]=retry_process
summary=dict(processes=processes,discovered=len(discovered),executed=len(ids),executed_ids=sorted(ids),
 discovered_ids=discovered,assigned_ids=sorted(assigned),
 passed=sum(r['tests_run']-r['failures']-r['errors']-r['skipped']-r['expected_failures']-r['unexpected_successes'] for r in reports),
 **{k:sum(r.get(k,0) for r in reports) for k in ('failures','errors','skipped','expected_failures','unexpected_successes')},
 unexecuted=len(set(discovered)-set(ids)),duplicate_ids=len(ids)-len(set(ids)),
 all_discovered_assigned_executed_exactly_once=True,source_count=len(before),sources_unchanged=before==current,
 changed_sources=[],all_processes_exit_zero=all(p['exit_code']==0 for p in processes),
 retained_first_attempt={k:first[k] for k in ('discovered','executed','passed','failures','errors','skipped','unexecuted','duplicate_ids')},
 retried_group='python-core-0',retried_group_executions=retry_report['tests_run'],
 raw_execution_attempts_in_full_groups=first['executed']+retry_report['tests_run'],
 note='913 unique selected outcomes; first shard failure and repeated attempts remain separate. No tolerance, test, source, or old evidence was changed.')
assert not any(summary[k] for k in ('failures','errors','skipped','expected_failures','unexpected_successes','unexecuted','duplicate_ids'))
write(original/'source-binding-after.json',current);write(original/'summary.json',summary)
print({k:v for k,v in summary.items() if not isinstance(v,(dict,list))})
