import hashlib,json,re,subprocess,sys
from pathlib import Path
out=Path('work/c04-software-validation');study=Path('interim/2026-12-19-c04-transient-water-cause-diagnosis')
groups=('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence','python-c01-switch','python-c01-path','python-c01-cause')
missing=[g for g in groups if not (out/g/'python-test-inventory.json').exists()]
if missing:raise SystemExit('NOT_READY: '+str(missing))
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rust_ids(path,listing):
 target=None;ids=[]
 for line in path.read_text().splitlines():
  match=re.search(r'Running (?:unittests )?(\S+\.rs) \(',line)
  if match:target=match[1]
  if 'Doc-tests' in line:target='doctests'
  match=re.match(r'(.+): test$',line.strip()) if listing else re.search(r'\btest (\S+) \.\.\. (ok|FAILED|ignored)(?:\s|$)',line)
  if match:
   if target is None:raise ValueError('Rust target absent')
   ids.append(target+'::'+match[1])
 return ids
listed=rust_ids(out/'rust-list.log',True)+rust_ids(out/'rust-example-list.log',True)
executed=rust_ids(out/'rust-all.log',False)+rust_ids(out/'rust-example-all.log',False)
assert len(set(listed))==len(listed)==466,(len(listed),len(set(listed)))
assert len(set(executed))==len(executed)==466,(len(executed),len(set(executed)))
assert sorted(listed)==sorted(executed),set(listed)-set(executed)
rust_summaries=[]
for name in ('rust-all','rust-example-all'):
 for row in re.findall(r'test result: (\w+)\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out',(out/(name+'.log')).read_text()):
  rust_summaries.append(dict(zip(('result','passed','failed','ignored','measured','filtered_out'),[row[0],*map(int,row[1:])])))
assert sum(r['passed'] for r in rust_summaries)==466
assert all(r['result']=='ok' and all(r[k]==0 for k in ('failed','ignored','measured','filtered_out')) for r in rust_summaries)
python={};all_executed=[];discovered=None
for group in groups:
 d=read(out/group/'python-test-inventory.json')
 if discovered is None:discovered=d['discovered_test_ids']
 assert d['discovered_test_ids']==discovered
 assigned=d['all_group_test_ids'][group]
 assert d['tests_run']==d['assigned_total']==len(assigned)
 assert set(d['test_elapsed_seconds'])==set(assigned)
 for key in ('failures','errors','skipped','expected_failures','unexpected_successes','exit_code'):assert d[key]==0,(group,key,d[key])
 all_executed+=assigned
 python[group]={k:d[k] for k in ('tests_run','assigned_total','failures','errors','skipped','expected_failures','unexpected_successes','exit_code','elapsed_seconds')}
assert len(all_executed)==len(set(all_executed))==len(discovered)==777
assert sorted(all_executed)==discovered
root_frozen=read(out/'core-before/source-snapshot.json')
root_changed=[n for n,h in root_frozen.items() if sha(Path(n))!=h]
assert root_changed==[],root_changed
interim_frozen=read(out/'core-before/interim-sources.json')
interim_changed=[n for n,h in interim_frozen.items() if sha(Path(n))!=h]
assert interim_changed==[],interim_changed
rust_frozen=read(out/'rust-before/source-snapshot.json')
rust_names=[n for n in rust_frozen if n.endswith('.rs') or n in ('Cargo.toml','Cargo.lock','rust-toolchain.toml')]
assert not [n for n in rust_names if sha(Path(n))!=rust_frozen[n]]
coverage=dict(rust_discovered_ids=sorted(listed),rust_executed_ids=sorted(executed),python_discovered_ids=discovered,python_executed_ids=sorted(all_executed),
 rust_exact_unique_coverage=True,python_exact_disjoint_coverage=True)
(out/'all-test-id-coverage.json').write_text(json.dumps(coverage,indent=2)+'\n')
report=dict(schema='c04-cause-validation/1',platform='Linux managed cloud',starting_commit='06695fccd549e3e1a89b18a8f60975da3359ed24',
 scientific_payload_commit='6a23bb9706b901ff36f63368c03b670b9b8b3c17',initial_inventory_commit=(out/'core-before/head.txt').read_text().strip(),
 software_validation_passed=True,passed=False,c04_complete=False,a04_complete=False,c01_complete=False,c02_complete=False,c03_complete=False,
 nc01_complete=False,p6_complete=False,ab_items_completed=False,pc_executed=False,cloud_diagnostics_remaining=True,
 windows_ci='push後の新先端Windows CIは未確定',
 rust=dict(discovered=466,executed=466,success=466,failure=0,unexecuted=0,skip=0,ignored=0,filtered_out=0,doc_tests=0,standard=463,examples=3,unique_id_coverage_verified=True),
 python=dict(discovered=777,executed=777,success=777,failure=0,error=0,unexecuted=0,skip=0,expected_failure=0,unexpected_success=0,
 new_tests=8,tools_tests_included=True,disjoint_complete_inventory_verified=True,groups=python),
 checks=read(out/'rust-processes.json')+read(out/'rust-example-processes.json'),
 source_checks=dict(rust_sources_and_manifests=len(rust_names),rust_unchanged=True,core_root_snapshot_files=len(root_frozen),core_root_snapshot_unchanged=True,
 core_interim_snapshot_files=len(interim_frozen),core_interim_snapshot_unchanged=True,
 dedicated_groups_started_before_initial_inventory_registration=True,
 added_after_full_test_start_observation_helper='c04cause_displacement.py; executed against SHA-verified sealed histories at all 3 levels, identity defect zero; no existing tested code changed',
 formal_execution_sources_bound_before_start=True,preliminary_all_10_field_files_byte_identical_to_formal=True),
 test_id_coverage='validation/all-test-id-coverage.json',preservation='validation/preservation-before-validation.json',
 scientific_solve_calls=dict(formal=12,formal_accepted_states=7008,preliminary_retained=10,preliminary_accepted_states=4832,total=22,total_accepted_states=11840,refused_cli_calls=3),
 sealed_verification=read(Path('work/c04-sealed-verification/verification.json')),
 new_test_smoke=dict(success=8,failure=0,error=0,skip=0,elapsed_seconds=.073,counted_again_in_full_total=False),
 initial_pre_fixture_tests=dict(success=3,failure=0,error=5,skip=0,counted_in_final_total=False,log='validation/preparation/new-tests-pre-fixture.log'),
 preparation_failures_retained=True,original_criteria_or_references_changed=False,
 displacement_contributions=dict(report_sha256=sha(study/'displacement-contributions.json'),source_sha256=sha(study/'c04cause_displacement.py'),
 maximum_identity_defect=max(x['maximum_identity_defect'] for x in read(study/'displacement-contributions.json')['levels'].values())),
 validation_snapshot_precedes_its_pinned_registration=True,final_registration_proof='validation/post-registration.json',
 retained_initial_full_core_attempts=dict(core0=dict(executed=366,success=365,failure=0,error=0,skip=1,exit_code=1),core1=dict(executed=387,success=386,failure=0,error=0,skip=1,exit_code=1),counted_in_final_totals=False),
 powershell=read(out/'powershell-runtime.json'))
with (study/'validation.json').open('x',encoding='utf-8',newline='\n') as f:json.dump(report,f,ensure_ascii=False,indent=2);f.write('\n')
print('FINAL',report['rust'],report['python']['discovered'])
