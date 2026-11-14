"""Count actual IDs and outcomes; software success never changes B-02 status."""
import collections,hashlib,json,pathlib,re,shutil
root=pathlib.Path(__file__).resolve().parents[1];work=root/'work/b02cause-validation';dest=root/'interim/2026-12-24-b02-offset-composite-cause-diagnosis'
def read(p):return json.loads(pathlib.Path(p).read_text(encoding='utf-8'))
rustdis=read(work/'rust-discovery/inventory.json');summaries=[];rustids=[]
executable_targets={str(pathlib.Path(x['executable']).resolve()):x['target'] for x in rustdis['artifacts']}
for name in ('rust-main','rust-examples'):
 target=None
 for line in (work/(name+'.log')).read_text().splitlines():
  m=re.search(r'Running .* \(([^)]+)\)',line)
  if m:target=executable_targets.get(str((root/m[1]).resolve()))
  t=re.match(r'test (.+) \.\.\. (ok|FAILED|ignored)$',line)
  if t and target:rustids.append(target+'::'+t[1])
  m=re.search(r'test result: (\w+)\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out',line)
  if m:summaries.append(dict(command=name,**dict(zip(('result','passed','failed','ignored','measured','filtered_out'),[m[1],*map(int,m.groups()[1:])]))))
rt={k:sum(r[k] for r in summaries) for k in ('passed','failed','ignored','measured','filtered_out')}
assert len(rustids)==len(set(rustids))
rust_missing=sorted(set(rustdis['test_ids'])-set(rustids));assert set(rustids)<=set(rustdis['test_ids'])
pyroot=work/'python-full';full=read(pyroot/'summary.json');groups={};executed=[];assigned=[];discovery=None;counts=collections.Counter();setup_errors=0
for row in full['groups']:
 name=row['group'];g=read(pyroot/name/'python-test-inventory.json');groups[name]=g
 ids=g['discovered_test_ids'];discovery=ids if discovery is None else discovery;assert discovery==ids
 executed+=list(g['test_elapsed_seconds']);assigned+=g['all_group_test_ids'][name]
 for k in ('tests_run','failures','errors','skipped','expected_failures','unexpected_successes'):counts[k]+=g[k]
 log=(pyroot/name/'python-tests.log').read_text()
 setup_errors+=len(re.findall(r'^ERROR: (?:setUpClass|tearDownClass|setUpModule|tearDownModule)\b',log,re.M))
assert len(executed)==len(set(executed)) and len(assigned)==len(set(assigned)) and sorted(assigned)==discovery
notexecuted=sorted(set(discovery)-set(executed));body_errors=counts['errors']-setup_errors
checks={n:read(work/(n+'-process.json')) for n in ('fmt','clippy','build','rust-discovery','rust-main','rust-examples')}
software_ok=all(c['exit_code']==0 for c in checks.values()) and all(r['exit_code']==0 for r in full['groups']) and not notexecuted and not rust_missing and all(counts[k]==0 for k in ('failures','errors','skipped','expected_failures','unexpected_successes')) and rt['failed']==rt['ignored']==rt['filtered_out']==0 and full['only_expected_post_solve_changes'] and full['core_sources_unchanged']
baseline=read(root/'work/b02cause-preparation/baseline-preservation.json');restoration=read(root/'work/b02cause-preparation/archive-restoration.json')
result=dict(schema='b02cause-validation/1',scope='software checks only; B-02 numerical acceptance remains false',diagnostic_only=True,
 passed=software_ok,software_checks_passed=software_ok,b02_complete=False,fixed_base_commit=baseline['base'],
 formal_code_commit='c792a4ec7c93800fd00c458ab2c178a4bef03dd0',tested_core_code_commit=(root/'work/b02cause-registered-code-commit.txt').read_text().strip(),postprocessing_source_sha256=read(dest/'measurements.json')['postprocessing_source_sha256'],checks=checks,
 rust=dict(discovered=rustdis['total'],executed=len(rustids),success=rt['passed'],failed=rt['failed'],not_executed=len(rust_missing),skip=rt['ignored'],ignored=rt['ignored'],filtered_out=rt['filtered_out'],measured=rt['measured'],doc_tests=0,missing_ids=rust_missing,test_ids=rustdis['test_ids'],summaries=summaries),
 python=dict(discovered=len(discovery),executed=counts['tests_run'],success=counts['tests_run']-counts['failures']-body_errors-counts['skipped']-counts['expected_failures']-counts['unexpected_successes'],failed=counts['failures'],errors=counts['errors'],setup_errors=setup_errors,body_errors=body_errors,not_executed=len(notexecuted),skip=counts['skipped'],expected_failures=counts['expected_failures'],unexpected_successes=counts['unexpected_successes'],missing_ids=notexecuted,duplicate_ids=0,test_ids=discovery,groups=groups),
 source_binding=full,formal_diagnostic_sources_unchanged=read(root/'work/b02cause-formal/source-binding-before.json')==read(root/'work/b02cause-formal/source-binding-after.json'),
 separate_smokes_not_double_counted=dict(initial_selector_executions=3,initial_selector_success=3,initial_group_executions=6,initial_group_success=6,final_new_contract_executions=8,final_new_contract_success=8,failed=0,skip=0,input_only_reference_paths=1,reference_states=20),
 baseline_restoration=dict(files=baseline['files'],total_bytes=sum(v['size'] for v in baseline['rows']),archive_files=restoration['archives'],archive_errors=len(restoration['errors']),missing=len(baseline['missing']),git_blob_mismatches=len(baseline['differences'])),
 archive_reproduction=read(dest/'validation/archive-reproduction.json'),preparations_and_failures=read(root/'work/b02cause-preparation/preparation-events.json'),
 production_rust_changed=False,nc_code_or_settings_changed=False,old_values_rewritten=False,numerical_tolerances_relaxed=False,skipped_to_bypass_failures=False,
 residual_768_diagnosis_repeated=False,new_maximum_tet10=96,cloud_diagnostics_remaining=True,pc_executed=False,windows_ci_after_push='未確定')
(dest/'validation.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
shutil.copytree(work,dest/'validation/full',dirs_exist_ok=False)
for name in ('b02cause_checks.py','b02cause_python_full.py','b02cause_collect_validation.py','b02cause_register.py','b02cause_commit.py','b02cause_restore.py','b02cause_preserve.py'):
 shutil.copy2(root/'work'/name,dest/'validation'/('executed-'+name))
print(dict(passed=software_ok,rust=rt,python=dict(counts),rust_ids=len(rustids),python_ids=len(discovery)))
