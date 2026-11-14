import pathlib,json,re,hashlib,shutil,subprocess,sys,platform,importlib.metadata
root=pathlib.Path('/workspace/solidcraft_soil_analysis');src=root/'work/a02-validation';dest=root/'interim/2026-12-21-a02-joint-birth-cause-diagnosis'
commands=json.loads((src/'rust-commands.json').read_text());discovery=json.loads((src/'rust-discovery.json').read_text());rust=[]
for label in ('rust','rust-examples'):
 for row in re.findall(r'test result: (\w+)\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out',(src/(label+'.log')).read_text()):
  rust.append(dict(zip(('result','passed','failed','ignored','measured','filtered_out'),[row[0],*map(int,row[1:])])))
r={k:sum(x[k] for x in rust) for k in ('passed','failed','ignored','measured','filtered_out')};r.update(discovered=sum(v['discovered_count'] for v in discovery.values()),skipped=0,result_groups=len(rust));r['unexecuted']=r['discovered']-sum(r[k] for k in ('passed','failed','ignored','filtered_out'))
groups=('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause')
reports={g:json.loads((src/g/'python-test-inventory.json').read_text()) for g in groups}
ids=set(next(iter(reports.values()))['discovered_test_ids']);assigned=[];executed=[];counts=dict(success=0,failed_test_cases=0,test_errors=0,fixture_errors=0,runner_errors=0,skipped=0,expected_failures=0,unexpected_successes=0,unexecuted=0)
for g,d in reports.items():
 assert set(d['discovered_test_ids'])==ids
 names=d['all_group_test_ids'][g];assigned+=names;actual=list(d['test_elapsed_seconds']);executed+=actual
 assert set(actual)<=set(names) and len(actual)==d['tests_run']
 log=(src/g/'python-tests.log').read_text();fixture=len(re.findall(r'^ERROR: (?:setUpClass|tearDownClass|setUpModule|tearDownModule)\b',log,re.M));testerrors=d['errors']-fixture
 success=d['tests_run']-d['failures']-testerrors-d['skipped']-d['expected_failures']-d['unexpected_successes'];assert success>=0
 counts['success']+=success;counts['failed_test_cases']+=d['failures']+testerrors;counts['test_errors']+=testerrors;counts['fixture_errors']+=fixture;counts['runner_errors']+=d['errors']
 for k in ('skipped','expected_failures','unexpected_successes'):counts[k]+=d[k]
 counts['unexecuted']+=len(names)-len(actual)
assert len(assigned)==len(set(assigned)) and set(assigned)==ids
assert len(executed)==len(set(executed))
counts.update(discovered=len(ids),executed=len(executed),missing_test_ids=sorted(ids-set(executed)),all_discovered_tests_assigned_exactly_once=True,tools_tests_included_without_double_counting=True)
post=json.loads((src/'python-post-rerun-source-check.json').read_text());assert post['all_sha256_identical']
software=all(c['exit_code']==0 for c in commands) and r['unexecuted']==0 and r['failed']==r['ignored']==r['filtered_out']==0 and all(d['exit_code']==0 for d in reports.values())
software=software and counts['unexecuted']==counts['runner_errors']==counts['failed_test_cases']==counts['skipped']==counts['expected_failures']==counts['unexpected_successes']==0
preservation=json.loads((dest/'base-preservation.json').read_text());assert preservation['all_protected_numeric_evidence_unchanged']
shutil.copytree(src,dest/'validation')
shutil.copy2(root/'work/a02-powershell/hashes.sha256',dest/'validation/powershell-official-hashes.sha256')
for name in ('a02-materialization.log','a02-materialization.json','a02-python-all-launcher.log','a02-rust-inventory-launcher.log','a02-preservation-launcher.log'):
 shutil.copy2(root/'work'/name,dest/'validation'/name)
shutil.copy2(root/'work/a02-preparation/notes.txt',dest/'validation/preparation-command-notes.txt')
for p in pathlib.Path('/workspace/work/a02').glob('*.py'):shutil.copy2(p,dest/'validation'/('executed-'+p.name))
summary=dict(schema='a02-joint-birth-cause-validation/1',base_commit='7b254b1cdbe9c288989631944cd698aa861f83ac',
 tested_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),passed=False,a02_complete=False,
 software_tests_passed=software,rust=r,python=counts,rust_commands=commands,
 python_groups={g:{k:d[k] for k in ('assigned_total','tests_run','failures','errors','skipped','expected_failures','unexpected_successes','exit_code','elapsed_seconds')} for g,d in reports.items()},
 python_test_source_check=post,rust_sources_unchanged_from_base=True,
 formal_diagnostic_source_check=json.loads((root/'work/a02-formal/source-unchanged.json').read_text()),
 archive_verification=(src/'archive-verification.log').read_text().strip(),
 sealed_derived_recomputation=json.loads((src/'rebuilt-derived.json').read_text()),
 old_preservation=preservation,
 initial_attempts=json.loads((src/'retained-initial-attempts.json').read_text()),
 powershell_runtime=json.loads((src/'powershell-runtime.json').read_text()),
 prior_smoke_tests_excluded_from_unique_totals=True,discovery_only_runs_excluded_from_success_counts=True,
 skipped_to_bypass_failures=False,numerical_tolerances_relaxed=False,old_values_rewritten=False,nc_code_or_settings_changed=False,
 a01_complete=False,a03_a04_complete=False,b_items_complete=False,c01_c04_complete=False,nc01_complete=False,p6_complete=False,
 pc_executed=False,cloud_diagnostics_remaining=True,windows_ci_after_push='未確定',
 environment=dict(python=sys.version,platform=platform.platform(),numpy=importlib.metadata.version('numpy'),scipy=importlib.metadata.version('scipy')))
(dest/'validation.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
manifest={p.relative_to(dest).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((dest/'validation').rglob('*')) if p.is_file()}
(dest/'validation-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(dict(software_tests_passed=software,rust=r,python=counts),indent=2))
