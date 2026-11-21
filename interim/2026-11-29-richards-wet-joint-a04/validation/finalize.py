"""Assemble actual full-suite counts; keep failed class setup as unexecuted IDs."""
import hashlib,json,re,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent.parent
V=HERE/'validation/full-01'

def read(p):return json.loads(p.read_text())
def write(p,v):
 with p.open('x') as f:json.dump(v,f,ensure_ascii=False,indent=2);f.write('\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
processes=read(V/'processes.json');assert set(processes)=={'fmt','clippy','rust-test','rust-examples','python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence'}
rust={};pattern=r'test result: (\w+)\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out'
for label in ('rust-test','rust-examples'):
 groups=[dict(zip(('result','passed','failed','ignored','measured','filtered_out'),[m[0],*map(int,m[1:])])) for m in re.findall(pattern,(V/(label+'.log')).read_text())]
 assert groups and processes[label]['exit_code']==0
 rust[label]=dict(groups=groups,**{k:sum(g[k] for g in groups) for k in ('passed','failed','ignored','measured','filtered_out')})
python={};assigned=[];executed=[];not_executed=[];inventory=None
for group in ('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence'):
 d=read(V/group/'python-test-inventory.json');python[group]=d
 if inventory is None:inventory=d['discovered_test_ids']
 assert d['discovered_test_ids']==inventory
 ids=d['all_group_test_ids'][group];assigned+=ids
 if d['tests_run']==len(ids):executed+=ids
 else:
  assert d['tests_run']==0 and d['errors']==1 and group in ('python-mean-evidence','python-spatial-evidence')
  not_executed+=ids
assert len(assigned)==len(set(assigned)) and sorted(assigned)==inventory
assert len(executed)==len(set(executed)) and set(executed).isdisjoint(not_executed) and sorted(executed+not_executed)==inventory
binding=read(HERE/'validation/source-binding-during-tests.json');changed=[f for f,v in binding['files'].items() if sha(ROOT/f)!=v];assert not changed,changed
fp=read(ROOT/'work/a04/formal-soil-03/source-fingerprints.json');changed=[f for f,v in fp.items() if sha(ROOT/f)!=v];assert not changed,changed
for d in read(HERE/'document-prefix-preservation.json'):
 original=subprocess.check_output(['git','show','04dd721ab82dd2c044569edf9e0b775b605cf2ef:'+d['path']],cwd=ROOT)
 assert (ROOT/d['path']).read_bytes().startswith(original) and hashlib.sha256(original).hexdigest()==d['original_sha256']
summary=dict(passed=False,all_requested_commands_attempted_to_completion=True,all_checks_successful=False,
 fmt_success=processes['fmt']['exit_code']==0,clippy_success=processes['clippy']['exit_code']==0,
 rust=rust,rust_passed=sum(v['passed'] for v in rust.values()),rust_failed=sum(v['failed'] for v in rust.values()),
 rust_ignored=sum(v['ignored'] for v in rust.values()),rust_filtered_out=sum(v['filtered_out'] for v in rust.values()),
 python_discovered=len(inventory),python_executed=len(executed),python_passed=len(executed),python_failures=sum(d['failures'] for d in python.values()),
 python_setup_errors=sum(d['errors'] for d in python.values()),python_not_executed=len(not_executed),python_skipped=sum(d['skipped'] for d in python.values()),
 python_expected_failures=sum(d['expected_failures'] for d in python.values()),python_unexpected_successes=sum(d['unexpected_successes'] for d in python.values()),
 python_executed_ids=sorted(executed),python_not_executed_ids=sorted(not_executed),all_assigned_once=True,processes=processes,
 historical_NC=dict(mean='saved mean-stress/complement partition differs from own history',spatial='remeasured JSON differs',allow_platform_float_differences=False,
 unchanged_original_code_and_archives=True,past_unexecuted_records_preserved=True,NC_numerical_acceptance=False),
 source_bindings=dict(fresh_solve_pre_run_files=len(fp),fresh_solve_pre_run_files_unchanged=True,full_suite_during_run_files=len(binding['files']),
 full_suite_during_run_files_unchanged=True,note='Full-suite binding was recorded during tests, not before start. Production/test Git blobs remain those of the specified start. Late supplemental scripts are not fresh-solve inputs.'),
 windows_ci='new pushed head not yet confirmed; no old success/timeout/cancellation transferred')
assert summary['rust_failed']==summary['rust_ignored']==summary['rust_filtered_out']==0
assert summary['python_executed']==565 and summary['python_not_executed']==8 and summary['python_setup_errors']==2 and summary['python_skipped']==0
write(HERE/'validation.json',summary)
write(HERE/'validation/final-source-binding.json',dict(files=binding['files'],all_unchanged=True,fresh_pre_run_files=len(fp),all_fresh_source_files_unchanged=True))
print(json.dumps({k:summary[k] for k in ('rust_passed','python_discovered','python_executed','python_not_executed','python_setup_errors','passed')}))
