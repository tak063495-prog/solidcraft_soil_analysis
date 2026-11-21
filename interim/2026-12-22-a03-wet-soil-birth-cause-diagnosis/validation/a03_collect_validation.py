import collections,hashlib,json,pathlib,re,subprocess
root=pathlib.Path.cwd();work=root/'work/a03-validation';dest=root/'interim/2026-12-22-a03-wet-soil-birth-cause-diagnosis'
def read(p):return json.loads(pathlib.Path(p).read_text(encoding='utf-8'))
rustdis=read(work/'rust-discovery/inventory.json');summaries=[];rustids=[]
executable_targets={str(pathlib.Path(x['executable']).resolve()):x['target'] for x in rustdis['artifacts']}
for name in ('rust-main','rust-examples'):
 text=(work/(name+'.log')).read_text();target=None
 for line in text.splitlines():
  match=re.search(r'Running .* \(([^)]+)\)',line)
  if match:
   target=executable_targets.get(str((root/match[1]).resolve()))
  test=re.match(r'test (.+) \.\.\. (ok|FAILED|ignored)$',line)
  if test and target:rustids.append(target+'::'+test[1])
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
assert read(work/'discovery-before/python-test-inventory.json')['discovered_test_ids']==discovery
notexecuted=sorted(set(discovery)-set(executed))
checks={n:read(work/(n+'-process.json')) for n in ('fmt','clippy','build','rust-main','rust-examples')}
software_ok=all(c['exit_code']==0 for c in checks.values()) and all(r['exit_code']==0 for r in full['groups']) and not notexecuted and all(counts[k]==0 for k in ('failures','errors','skipped','expected_failures','unexpected_successes')) and rt['failed']==rt['ignored']==rt['filtered_out']==0 and full['sources_unchanged']
# Successful-test counts are exact for the all-pass final run. Keep fixtures
# separate from started tests if a future reproduction has setup errors.
assert counts['failures']==counts['errors']==counts['skipped']==counts['expected_failures']==counts['unexpected_successes']==0,'inspect raw inventories for failed-test versus fixture-error accounting'
report=dict(schema='a03-cause-validation/1',scope='software regression checks only; not A-03 numerical acceptance',diagnostic_only=True,passed=software_ok,software_checks_passed=software_ok,a03_complete=False,
            fixed_base_commit='9a5e9d182299b6fce64cc18ea61a65be4702c4e1',tested_code_commit='681427c73a1998129b7b5f1625084b8ae3742a16',
            checks=checks,rust=dict(discovered=rustdis['total'],executed=len(rustids),success=rt['passed'],failed=rt['failed'],not_executed=rustdis['total']-len(rustids),skip=rt['ignored'],ignored=rt['ignored'],filtered_out=rt['filtered_out'],measured=rt['measured'],doc_tests=0,test_ids=rustdis['test_ids'],summaries=summaries),
            python=dict(discovered=len(discovery),executed=counts['tests_run'],success=counts['tests_run'],failed=counts['failures'],errors=counts['errors'],not_executed=len(notexecuted),skip=counts['skipped'],expected_failures=counts['expected_failures'],unexpected_successes=counts['unexpected_successes'],missing_ids=notexecuted,duplicate_ids=0,test_ids=discovery,groups=groups),
            source_count_before_after=full['source_count'],code_sources_unchanged_before_after=full['sources_unchanged'],
            separate_smoke_tests_not_double_counted=dict(preparation=6,a03_light=10),
            original_and_additional_refusals_are_not_software_test_success_counts=True,
            archive_reproduction=read(work/'archive-reproduction.json'),
            old_archive_restoration=dict(files=1590,bytes=28021593924,missing=0,Git_blob_mismatches=0),
            preparations_and_infrastructure_failures=read(root/'work/preparation-events.json'),
            old_values_rewritten=False,numerical_tolerances_relaxed=False,nc_code_or_settings_changed=False,skipped_to_bypass_failures=False,
            cloud_diagnostics_remaining=True,pc_executed=False,windows_ci_after_push='未確定')
(dest/'validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(dict(passed=software_ok,rust=rt,python=dict(counts),rust_ids=len(rustids),python_ids=len(discovery),sources=full['source_count']))
