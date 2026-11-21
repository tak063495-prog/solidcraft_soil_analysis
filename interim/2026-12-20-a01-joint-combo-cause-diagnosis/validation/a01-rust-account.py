import json,re
from pathlib import Path
root=Path('/workspace/work/a01-full-validation')
def parse(path,listed):
 scope=None;ids=[];outcomes=[]
 for line in path.read_text().splitlines():
  m=re.match(r'\s*Running (.+?) \(',line)
  if m:scope=m[1]
  if listed:m=re.fullmatch(r'([^\s]+): test',line)
  else:m=re.match(r'test ([^\s]+) \.\.\. (ok|FAILED|ignored)\s*$',line)
  if m:
   assert scope,(path,line)
   ids.append(scope+'::'+m[1])
   if not listed:outcomes.append(m[2])
 return ids,outcomes
reports=[]
for run,listing in [('rust-tests.log','rust-discovered.log'),('rust-example-tests.log','rust-examples-discovered.log')]:
 actual,outcomes=parse(root/run,False);discovered,_=parse(root/listing,True)
 assert len(actual)==len(set(actual)) and len(discovered)==len(set(discovered))
 assert sorted(actual)==sorted(discovered),(len(actual),len(discovered),set(discovered)-set(actual))
 summaries=[dict(result=m[0],passed=int(m[1]),failed=int(m[2]),ignored=int(m[3]),measured=int(m[4]),filtered_out=int(m[5])) for m in re.findall(r'test result: (\w+)\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out',(root/run).read_text())]
 assert sum(x['passed'] for x in summaries)==outcomes.count('ok')
 reports.append(dict(log=run,discovered_test_ids=sorted(discovered),executed_test_ids=sorted(actual),discovered=len(discovered),passed=outcomes.count('ok'),failed=outcomes.count('FAILED'),ignored=outcomes.count('ignored'),not_run=0,skipped=0,filtered_out=sum(x['filtered_out'] for x in summaries),result_groups=summaries))
assert sum(x['passed'] for x in reports)==467
report=dict(discovered=467,executed=467,passed=467,failed=0,ignored=0,not_run=0,skipped=0,filtered_out=0,all_discovered_ids_executed_once=True,normal_tests=464,example_tests=3,groups=reports)
(root/'rust-test-accounting.json').write_text(json.dumps(report,indent=2)+'\n')
print('Rust',report['passed'],'unique tests, discovery/execution identical')
