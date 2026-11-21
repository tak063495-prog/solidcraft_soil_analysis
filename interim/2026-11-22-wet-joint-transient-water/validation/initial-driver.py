import datetime,hashlib,json,os,re,subprocess,sys,time
from pathlib import Path
root=Path.cwd();out=root/'interim/2026-11-22-wet-joint-transient-water/validation';out.mkdir(exist_ok=False)
tracked=subprocess.check_output(['git','ls-files','-z']).decode().split('\0')
paths=sorted(set([p for p in tracked if p and (p.endswith(('.py','.rs','.toml','.lock','.csv')) or p.startswith(('examples/','tests/data/','.github/')) and p.endswith(('.json','.yml','.yaml')))] +[str(p.relative_to(root)) for p in (root/'interim/2026-11-22-wet-joint-transient-water').glob('*.py')]+['examples/wet_joint_transient_probe.rs','tests/test_wet_joint_transient_water.py']))
def hashes():return {p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in paths}
frozen=hashes();(out/'source-before.json').write_text(json.dumps(frozen,indent=2)+'\n')
commands=[('fmt',['cargo','fmt','--all','--check']),('clippy',['cargo','clippy','--locked','--all-targets','--','-D','warnings']),('build',['cargo','build','--locked','--release']),('rust-tests',['cargo','test','--locked','--release'])]
for group in ['python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence']:
 commands.append((group,[sys.executable,'tools/run_ci_python_tests.py','--group',group,'--output',str(out/group)]))
rows=[]
for name,cmd in commands:
 start=time.monotonic();print('start',name,flush=True)
 with (out/(name+'.console.log')).open('x') as stream:
  child=subprocess.Popen(cmd,stdout=stream,stderr=subprocess.STDOUT)
  code=child.wait()
 row=dict(name=name,command=cmd,pid=child.pid,exit_code=code,elapsed_seconds=time.monotonic()-start);rows.append(row)
 print('end',row,flush=True)
 (out/'command-progress.json').write_text(json.dumps(rows,indent=2)+'\n')
after=hashes();(out/'source-after.json').write_text(json.dumps(after,indent=2)+'\n')
text=(out/'rust-tests.console.log').read_text();rust=[tuple(map(int,x)) for x in re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out',text)]
inventories=[json.loads((out/g/'python-test-inventory.json').read_text()) for g in ['python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence']]
ids=[x for v in inventories for x in v['all_group_test_ids'][v['group']]]
summary=dict(complete=all(r['exit_code']==0 for r in rows),commands=rows,source_files=len(paths),source_hashes_unchanged=frozen==after,
 rust_total_passed=sum(x[0] for x in rust),rust_summary_groups=len(rust),rust_failed=sum(x[1] for x in rust),rust_ignored=sum(x[2] for x in rust),rust_filtered_out=sum(x[4] for x in rust),
 python_total=sum(v.get('tests_run',0) for v in inventories),python_group_totals={v['group']:v.get('tests_run') for v in inventories},python_discovered_total=inventories[0]['discovered_total'],python_id_union_exact=len(ids)==len(set(ids)) and sorted(ids)==inventories[0]['discovered_test_ids'],
 python_failures=sum(v.get('failures',0) for v in inventories),python_errors=sum(v.get('errors',0) for v in inventories),python_skipped=sum(v.get('skipped',0) for v in inventories),python_expected_failures=sum(v.get('expected_failures',0) for v in inventories),python_unexpected_successes=sum(v.get('unexpected_successes',0) for v in inventories),
 numerical_C04_accepted=False,Windows_CI_status='未確定')
(out.parent/'validation.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)
assert summary['complete'] and summary['source_hashes_unchanged'] and summary['python_id_union_exact']
