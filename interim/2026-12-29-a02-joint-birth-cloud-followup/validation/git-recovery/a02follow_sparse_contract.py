import json,sys,subprocess,hashlib
from pathlib import Path
r=Path.cwd();sys.path.insert(0,str(r/'tools'));import run_ci_python_evidence_shards as ci
by_id,groups=ci.discover(r);new={g:[i for i in ids if i.startswith('test_a02_joint_birth_cloud_followup.')] for g,ids in groups.items()}
d=r/'interim/2026-12-29-a02-joint-birth-cloud-followup';idx=json.loads((d/'archive-index.json').read_text());parts=[(d/p['name']).relative_to(r).as_posix() for p in idx['parts']]
old=json.loads((r/'interim/2026-12-28-a01-joint-combo-cloud-followup/validation/sparse-checkout-contract.json').read_text());proof={}
for group in ('python-core-0','python-core-1'):
 rules=old['actual_core_sparse_rules'][group]['rules'];file=r/'work'/('a02follow-'+group+'-sparse-rules.txt');file.write_text('\n'.join(rules)+'\n')
 run=subprocess.run(['git','sparse-checkout','check-rules','--no-cone','--rules-file',str(file)],input='\n'.join(parts)+'\n',text=True,capture_output=True,check=True)
 selected=run.stdout.splitlines();assert selected==parts
 proof[group]=dict(rules=rules,parts=parts,all_selected=True)
assert sum(map(len,new.values()))==10 and not any(new[g] for g in ci.GROUPS if g not in ('python-core-0','python-core-1'))
v=dict(discovered_total=len(by_id),new_test_count=10,assignment=new,actual_core_sparse_rules=proof,rust_integration_1_added_heavy_tests=0,workflow_sha256=hashlib.sha256((r/'.github/workflows/fem-ci.yml').read_bytes()).hexdigest(),discovery_is_not_execution=True)
(d/'validation/sparse-checkout-contract.json').write_text(json.dumps(v,indent=2)+'\n');print('discovered',len(by_id),'new allocation',{g:len(x) for g,x in new.items()})
