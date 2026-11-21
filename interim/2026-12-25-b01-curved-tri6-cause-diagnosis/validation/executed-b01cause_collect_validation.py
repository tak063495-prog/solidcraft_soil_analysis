import collections,hashlib,json,pathlib,re,shutil,subprocess
root=pathlib.Path.cwd();work=root/'work/b01cause-validation';dest=root/'interim/2026-12-25-b01-curved-tri6-cause-diagnosis'
def read(p):return json.loads(pathlib.Path(p).read_text(encoding='utf-8'))
rustdis=read(work/'rust-discovery/inventory.json');mapping={str(pathlib.Path(a['executable']).resolve()):a['target'] for a in rustdis['artifacts']};rustids=[];summaries=[]
for name in ('rust-main','rust-examples'):
 target=None
 for line in (work/(name+'.log')).read_text().splitlines():
  m=re.search(r'Running .* \(([^)]+)\)',line)
  if m:target=mapping.get(str((root/m[1]).resolve()))
  m=re.match(r'test (.+) \.\.\. (ok|FAILED|ignored)$',line)
  if m and target:rustids.append(target+'::'+m[1])
  m=re.search(r'test result: (\w+)\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out',line)
  if m:summaries.append(dict(command=name,**dict(zip(('result','passed','failed','ignored','measured','filtered_out'),[m[1],*map(int,m.groups()[1:])]))))
rt={k:sum(r[k] for r in summaries) for k in ('passed','failed','ignored','measured','filtered_out')}
assert len(rustids)==len(set(rustids));assert set(rustids)<=set(rustdis['test_ids']);rust_missing=sorted(set(rustdis['test_ids'])-set(rustids))
pyroot=work/'python-full';full=read(pyroot/'summary.json');groups={};executed=[];assigned=[];discovery=None;counts=collections.Counter();setup_errors=0
for row in full['groups']:
 name=row['group'];g=read(pyroot/name/'python-test-inventory.json');groups[name]=g;ids=g['discovered_test_ids'];discovery=ids if discovery is None else discovery;assert discovery==ids
 executed+=list(g.get('test_elapsed_seconds',{}));assigned+=g['all_group_test_ids'][name]
 for k in ('tests_run','failures','errors','skipped','expected_failures','unexpected_successes'):counts[k]+=g.get(k,0)
 setup_errors+=len(re.findall(r'^ERROR: (?:setUpClass|tearDownClass|setUpModule|tearDownModule)\b',(pyroot/name/'python-tests.log').read_text(),re.M))
assert len(executed)==len(set(executed)) and len(assigned)==len(set(assigned)) and sorted(assigned)==discovery
missing=sorted(set(discovery)-set(executed));body_errors=counts['errors']-setup_errors
checks={n:read(work/(n+'-process.json')) for n in ('fmt','clippy','build','rust-discovery','rust-main','rust-examples')}
ok=all(c['exit_code']==0 for c in checks.values()) and all(r['exit_code']==0 for r in full['groups']) and not missing and not rust_missing and all(counts[k]==0 for k in ('failures','errors','skipped','expected_failures','unexpected_successes')) and rt['failed']==rt['ignored']==rt['filtered_out']==0 and full['sources_unchanged']
before=read(root/'work/b01cause-probe/source-binding-before.json');after=read(root/'work/b01cause-probe/source-binding-after.json');baseline=read(root/'work/b01cause-baseline-preservation.json')
result=dict(schema='b01cause-validation/1',scope='software checks; B-01 remains unmet',diagnostic_only=True,passed=bool(ok),software_checks_passed=bool(ok),b01_complete=False,
 fixed_base_commit=baseline['base'],scientific_commit=(root/'work/b01cause-science-commit.txt').read_text().strip(),tested_code_commit=(root/'work/b01cause-registered-code-commit.txt').read_text().strip(),checks=checks,
 rust=dict(discovered=rustdis['total'],executed=len(rustids),success=rt['passed'],failed=rt['failed'],not_executed=len(rust_missing),skip=rt['ignored'],ignored=rt['ignored'],filtered_out=rt['filtered_out'],measured=rt['measured'],doc_tests=0,missing_ids=rust_missing,test_ids=rustdis['test_ids'],summaries=summaries),
 python=dict(discovered=len(discovery),executed=counts['tests_run'],success=counts['tests_run']-counts['failures']-body_errors-counts['skipped']-counts['expected_failures']-counts['unexpected_successes'],failed=counts['failures'],errors=counts['errors'],setup_errors=setup_errors,body_errors=body_errors,not_executed=len(missing),skip=counts['skipped'],expected_failures=counts['expected_failures'],unexpected_successes=counts['unexpected_successes'],missing_ids=missing,duplicate_ids=0,test_ids=discovery,groups=groups),
 full_python_source_binding=full,formal_source_binding=dict(source_count=len(before),unchanged=before==after,scope='production Rust plus the then-existing SPEC/geometry/build/probe/driver; later reporter and CLI comparison are separately committed postprocessing'),
 baseline_preservation=dict(files=baseline['files'],bytes=baseline['bytes'],missing=baseline['missing'],differences=baseline['differences'],mode_differences=baseline['mode_differences']),
 archive_reproduction=read(dest/'validation/archive-reproduction.json'),separate_smokes_not_double_counted=dict(new_python_contracts_success=11,failed=0,skip=0),
 preparations_and_retained_failures=dict(diagnostic_compile_failure=1,initial_saved_cli_comparison_false_controls=4,initial_comparison_cause='nested elapsed_seconds excluded only at top level; preserved initial false and separate corrected saved-field audit; no new solves',partial_git_transfers_and_initial_main_checkout_differences='retained; fixed target restored; all 8550 files verified'),
 production_rust_changed=False,nc_code_or_settings_changed=False,old_values_rewritten=False,numerical_tolerances_relaxed=False,skipped_to_bypass_failures=False,
 new_maximum_tri6=8,new_maximum_structure_dofs=150,cloud_diagnostics_remaining=True,specification_decisions_needed=True,pc_executed=False,windows_ci_after_push='未確定')
(dest/'validation.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n');shutil.copytree(work,dest/'validation/full',dirs_exist_ok=False)
for name in ('b01cause_rust_full.py','b01cause_python_full.py','b01cause_collect_validation.py','b01cause_register.py','b01cause_commit.py','b01cause_preserve.py'):
 shutil.copyfile(root/'work'/name,dest/'validation'/('executed-'+name))
print(dict(passed=ok,rust=rt,python=dict(counts),rust_ids=len(rustids),python_ids=len(discovery)))
