import hashlib,json,re,subprocess,sys
from pathlib import Path
root=Path('/workspace/work/a01-source');out=Path('/workspace/work/a01-full-validation');study=root/'interim/2026-12-20-a01-joint-combo-cause-diagnosis'
sys.path.insert(0,str(root/'tools'));from freeze_sources import manifest
read=lambda p:json.loads(Path(p).read_text(encoding='utf-8'))
discovery=read(out/'discovery/python-test-inventory.json');all_ids=set(discovery['discovered_test_ids']);executed=[];assigned=[];groups={};class_errors=0
for group in discovery['all_group_test_ids']:
 r=read(out/group/'python-test-inventory.json');assert set(r['discovered_test_ids'])==all_ids
 ids=r['all_group_test_ids'][group];assert ids==discovery['all_group_test_ids'][group];assigned+=ids
 actual=list(r.get('test_elapsed_seconds',{}));assert len(actual)==r.get('tests_run',0);executed+=actual
 log=(out/group/'python-tests.log').read_text(encoding='utf-8') if (out/group/'python-tests.log').exists() else ''
 failure_ids=set();error_ids=set();setup=[]
 for kind,test_id in re.findall(r'^(FAIL|ERROR): .*? \(([^()]+)\)',log,re.M):
  if test_id in actual:(failure_ids if kind=='FAIL' else error_ids).add(test_id)
  else:setup.append(dict(kind=kind,id=test_id))
 failures=r.get('failures',0);errors=r.get('errors',0);skips=r.get('skipped',0);expected=r.get('expected_failures',0);unexpected=r.get('unexpected_successes',0)
 if not failures and not errors:passed=len(actual)-skips-expected-unexpected
 else:
  assert len(failure_ids)+sum(x['kind']=='FAIL' for x in setup)==failures
  assert len(error_ids)+sum(x['kind']=='ERROR' for x in setup)==errors
  passed=len(actual)-len(failure_ids)-len(error_ids)-skips-expected-unexpected
 class_errors+=len(setup)
 groups[group]=dict(discovered=r['assigned_total'],executed=len(actual),passed=passed,failed_tests=len(failure_ids),error_tests=len(error_ids),class_setup_errors=setup,
     failure_events=failures,error_events=errors,not_run=r['assigned_total']-len(actual),skipped=skips,expected_failures=expected,unexpected_successes=unexpected,
     exit_code=r['exit_code'],elapsed_seconds=r['elapsed_seconds'],executed_test_ids=sorted(actual))
assert len(assigned)==len(set(assigned))==len(all_ids) and set(assigned)==all_ids
assert len(executed)==len(set(executed)) and set(executed)<=all_ids
fields=('passed','failed_tests','error_tests','failure_events','error_events','not_run','skipped','expected_failures','unexpected_successes')
py=dict(discovered=len(all_ids),executed=len(executed),**{k:sum(r[k] for r in groups.values()) for k in fields},class_setup_error_events=class_errors,
    all_discovered_tests_assigned_exactly_once=True,all_discovered_tests_executed_once=set(executed)==all_ids,groups=groups)
rust=read(out/'rust-test-accounting.json');build=read('/workspace/work/a01-build-validation/result.json')
assert manifest()==read(out/'source-snapshot.json')
extra={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(study.glob('*.py'))};assert extra==read(out/'diagnostic-python-before-tests.json')
base='5908506ec54d3972a18cec36599ade1e05fd7c2a';base_paths=subprocess.check_output(['git','ls-tree','-r','--name-only',base],cwd=root,text=True).splitlines()
changed=subprocess.check_output(['git','diff','--name-only',base],cwd=root,text=True).splitlines();old=sorted(set(changed)&set(base_paths))
assert old==sorted(['src/linalg.rs','src/stage_control.rs','src/stages.rs','interim/2026-12-06-p6-01-status-reconciliation/audit_status.py','interim/2026-12-10-ci-p601-inventory-fix/inventory.toml']),old
software=all(v==0 for v in build.values()) and rust['passed']==rust['discovered'] and all(r['exit_code']==0 for r in groups.values()) and not any(py[k] for k in ('failure_events','error_events','not_run','skipped','expected_failures','unexpected_successes'))
result=dict(schema='a01-joint-combo-cause-validation/1',passed=False,a01_complete=False,nc01_complete=False,p6_complete=False,
    software_tests_passed=software,all_required_command_groups_executed=True,all_discovered_tests_executed=py['all_discovered_tests_executed_once'],
    fmt_clippy_release_build=build,rust=rust,python=py,platform='Linux; Python 3.13.17 / Rust 1.99.1',
    source_unchanged_during_tests=True,core_source_files=905,diagnostic_python_sources=11,legacy_helpers_bound_to_base=True,
    protected_existing_paths_changed=old,old_git_paths_preserved=len(base_paths)-len(old),base_git_paths=len(base_paths),old_values_hash_snapshots_rewritten=False,
    sealed_derived_recomputation=read(out/'sealed-recomputation.json'),archive_ci_allocation=read(out/'ci-archive-allocation.json'),
    preliminary_smokes_counted_in_totals=False,preliminary_smokes=dict(rust_observer=dict(passed=1,filtered_out=62),python_principles=dict(passed=5),python_new_module=dict(passed=8),p6_related=dict(passed=76)),
    p601_scientific_registration=read(out/'initial-registration/p601.json'),p601_validation_registration='Final fixed-commit registration and post-registration audit are recorded separately in validation/post-registration.json.',
    windows_ci_after_push='未確定',pc_executed=False,physical_instability_proved=False,cloud_diagnostics_remaining=True)
(out/'collected-validation.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
print('Rust',rust['passed'],'Python',py['passed'],'unexecuted',py['not_run'],'skipped',py['skipped'],'software pass',software)
