"""Run existing CI groups plus explicit new checks; never skip or alter discovery.
This validation harness calls existing tests/CLI; reference generators do not.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

GROUPS=('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence')


def snapshot(root):
    paths=subprocess.check_output(['git','ls-tree','-r','--name-only','HEAD'],cwd=root,text=True).splitlines()
    selected=[p for p in paths if p.endswith(('.py','.rs','.toml','.yml','.yaml')) or p.startswith('requirements')]
    selected += [p.relative_to(root).as_posix() for p in Path(__file__).resolve().parent.glob('*.py')]
    missing=[p for p in selected if not (root/p).is_file()]
    if missing:
        raise ValueError('missing validation sources: '+str(missing))
    return {p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in sorted(set(selected))}


def command(args,root,log,env):
    t=time.monotonic()
    with log.open('x',encoding='utf-8') as f:
        r=subprocess.run(args,cwd=root,env=env,stdout=f,stderr=subprocess.STDOUT,check=False)
    return dict(command=list(map(str,args)),exit_code=r.returncode,elapsed_seconds=time.monotonic()-t,log=log.name)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[2];here=Path(__file__).resolve().parent
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONUTF8='1',PYTHONIOENCODING='utf-8')
    before=snapshot(root)
    report=dict(status='running',reference_acceptance=None,source_before=before,groups={},new_checks={})
    def save():
        (out/'validation.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    save()
    for group in GROUPS:
        print('Running '+group,flush=True)
        result=command([sys.executable,'tools/run_ci_python_tests.py','--group',group,'--output',str(out/group)],root,out/(group+'.log'),env)
        inventory=json.loads((out/group/'python-test-inventory.json').read_text())
        result['inventory']=inventory
        report['groups'][group]=result
        save()
    print('Running explicit independent checks',flush=True)
    report['new_checks']=command([sys.executable,str(here/'checks.py')],root,out/'new-checks.log',env)
    # Count tests through unittest loading without altering its default discovery.
    spec="import sys,unittest;sys.path.insert(0,sys.argv[1]);import checks;print(unittest.defaultTestLoader.loadTestsFromModule(checks).countTestCases())"
    report['new_checks']['test_count']=int(subprocess.check_output([sys.executable,'-c',spec,str(here)],cwd=root,env=env,text=True))
    after=snapshot(root)
    report['source_after']=after;report['sources_unchanged']=before==after
    invs=[report['groups'][g]['inventory'] for g in GROUPS]
    ids=[x for g,i in zip(GROUPS,invs) for x in i.get('all_group_test_ids',{}).get(g,[])]
    inventory=invs[0].get('discovered_test_ids',[])
    report['same_inventory_all_groups']=all(i.get('discovered_test_ids')==inventory for i in invs)
    report['all_executed_exactly_once']=sorted(ids)==sorted(inventory) and len(ids)==len(set(ids))
    report['existing_total']=sum(i.get('tests_run',0) for i in invs)
    report['existing_failures']=sum(i.get('failures',0) for i in invs)
    report['existing_errors']=sum(i.get('errors',0) for i in invs)
    report['existing_skipped']=sum(i.get('skipped',0) for i in invs)
    report['existing_expected_failures']=sum(i.get('expected_failures',0) for i in invs)
    report['existing_unexpected_successes']=sum(i.get('unexpected_successes',0) for i in invs)
    success=(report['sources_unchanged'] and report['same_inventory_all_groups'] and report['all_executed_exactly_once']
             and report['existing_total']==len(inventory)
             and all(report['groups'][g]['exit_code']==0 for g in GROUPS) and report['new_checks']['exit_code']==0)
    report['status']='software_checks_passed' if success else 'failed'
    save()
    print(json.dumps({k:report[k] for k in ('status','existing_total','existing_failures','existing_errors','existing_skipped')},ensure_ascii=False))
    return 0 if success else 1


if __name__=='__main__':
    sys.exit(main())
