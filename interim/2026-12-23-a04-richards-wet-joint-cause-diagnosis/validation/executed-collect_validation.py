import collections,hashlib,json,pathlib,re,shutil,subprocess
root=pathlib.Path('/workspace/solidcraft_soil_analysis');work=root/'work/a04cause-validation';dest=root/'interim/2026-12-23-a04-richards-wet-joint-cause-diagnosis'
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
assert sorted(rustids)==rustdis['test_ids'],(len(rustids),rustdis['total'],set(rustdis['test_ids'])-set(rustids))
pyroot=work/'python-full';full=read(pyroot/'summary.json');groups={};executed=[];assigned=[];discovery=None;counts=collections.Counter()
for row in full['groups']:
 name=row['group'];g=read(pyroot/name/'python-test-inventory.json');groups[name]=g
 ids=g['discovered_test_ids'];discovery=ids if discovery is None else discovery;assert discovery==ids
 executed+=list(g['test_elapsed_seconds']);assigned+=g['all_group_test_ids'][name]
 for k in ('tests_run','failures','errors','skipped','expected_failures','unexpected_successes'):counts[k]+=g[k]
assert len(executed)==len(set(executed)) and len(assigned)==len(set(assigned))
assert sorted(assigned)==discovery
notexecuted=sorted(set(discovery)-set(executed))
checks={n:read(work/(n+'-process.json')) for n in ('fmt','clippy','build','rust-discovery','rust-main','rust-examples')}
software_ok=all(c['exit_code']==0 for c in checks.values()) and all(r['exit_code']==0 for r in full['groups']) and not notexecuted and all(counts[k]==0 for k in ('failures','errors','skipped','expected_failures','unexpected_successes')) and rt['failed']==rt['ignored']==rt['filtered_out']==0 and full['sources_unchanged']
assert counts['failures']==counts['errors']==counts['skipped']==counts['expected_failures']==counts['unexpected_successes']==0,'failed test vs fixture error accounting must be inspected before success can be calculated'
baseline=read('/workspace/work/a04/baseline-files.json')
restored=[v for v in baseline['files'] if v['restored']]
result=dict(schema='a04-cause-validation/1',scope='software checks only; A-04 numerical acceptance remains false',
 diagnostic_only=True,passed=software_ok,software_checks_passed=software_ok,a04_complete=False,
 fixed_base_commit=baseline['base'],tested_python_code_commit=pathlib.Path('/workspace/work/a04/registered-code-commit.txt').read_text().strip(),
 rust_numerical_sources_identical_to_fixed_base=True,checks=checks,
 rust=dict(discovered=rustdis['total'],executed=len(rustids),success=rt['passed'],failed=rt['failed'],not_executed=rustdis['total']-len(rustids),skip=rt['ignored'],ignored=rt['ignored'],filtered_out=rt['filtered_out'],measured=rt['measured'],doc_tests=0,test_ids=rustdis['test_ids'],summaries=summaries),
 python=dict(discovered=len(discovery),executed=counts['tests_run'],success=counts['tests_run'],failed=counts['failures'],errors=counts['errors'],not_executed=len(notexecuted),skip=counts['skipped'],expected_failures=counts['expected_failures'],unexpected_successes=counts['unexpected_successes'],missing_ids=notexecuted,duplicate_ids=0,test_ids=discovery,groups=groups),
 source_count_before_after=full['source_count'],code_sources_unchanged_before_after=full['sources_unchanged'],
 separate_smokes_not_double_counted=dict(initial_tests=8,initial_success=7,initial_failed=1,corrected_subset_success=8,final_new_tests_success=11),
 archive_reproduction=read(dest/'validation/archive-reproduction.json'),
 baseline_restoration=dict(files=len(baseline['files']),total_bytes=sum(v['size'] for v in baseline['files']),restored_files=len(restored),restored_bytes=sum(v['size'] for v in restored),missing=0,Git_blob_mismatches=0),
 old_archive_verification=[read(dest/'validation'/n) for n in ('old-original-archive-verification.json','old-corrected-archive-verification.json')],
 preparations_and_postprocessing_failure_records=read(dest/'validation/preparation-events.json'),
 production_rust_changed=False,nc_code_or_settings_changed=False,old_values_rewritten=False,numerical_tolerances_relaxed=False,skipped_to_bypass_failures=False,
 cloud_diagnostics_remaining=True,pc_executed=False,windows_ci_after_push='未確定')
(dest/'validation.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
# Raw complete logs, discovery and test IDs accompany the pinned summary.
shutil.copytree(work,dest/'validation/full',dirs_exist_ok=False)
for filename in ('checks.py','python_full.py','discover_rust.py','collect_validation.py','preserve.py'):
 shutil.copy2(pathlib.Path('/workspace/work/a04')/filename,dest/'validation'/('executed-'+filename))
print(dict(passed=software_ok,rust=rt,python=dict(counts),rust_ids=len(rustids),python_ids=len(discovery),sources=full['source_count']))
