"""Unique test-ID accounting; retains initial errors and records preparation failures."""
import argparse
import json
from pathlib import Path
import re


def read(p):return json.loads(p.read_text())


def summarize(raw):
    processes=read(raw/'validation-rust/processes.json')
    text=(raw/'validation-rust/rust-test.log').read_text()
    groups=[tuple(map(int,x)) for x in re.findall(r'test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out',text)]
    rust=dict(processes=processes,passed_tests=sum(x[0] for x in groups),result_groups=len(groups),
        failed=sum(x[1] for x in groups),ignored=sum(x[2] for x in groups),filtered=sum(x[4] for x in groups))
    rust['passed']=set(processes)=={'fmt','clippy','build','rust-test','rust-doc'} and all(x['exit_code']==0 for x in processes.values()) and rust['passed_tests']==461 and not any(rust[k] for k in ('failed','ignored','filtered'))
    initial={g:read(raw/'validation-python'/g/'python-test-inventory.json') for g in
        ['python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence']}
    discovered=initial['python-core-0']['discovered_test_ids']
    assigned=[x for g in initial.values() for x in g['all_group_test_ids'][g['group']]]
    assert len(assigned)==len(set(assigned)) and sorted(assigned)==discovered
    success=set();not_executed=[];final_groups={};initial_errors={}
    for name,report in initial.items():
        initial_errors[name]=dict(tests_run=report['tests_run'],errors=report['errors'],failures=report['failures'],exit_code=report['exit_code'])
        ids=set(report['all_group_test_ids'][name])
        if name in ('python-core-0','python-core-1'):
            log=(raw/'validation-python'/name/'python-tests.log').read_text()
            class_errors=re.findall(r'ERROR: setUpClass \(([^)]+)\)',log)
            test_errors=[c if c.endswith('.'+m) else c+'.'+m
                         for m,c in re.findall(r'ERROR: (test\S+) \(([^)]+)\)',log)]
            absent={x for x in ids if any(x.startswith(c+'.') for c in class_errors)}
            before=ids-absent-set(test_errors)
            assert len(before)==report['tests_run']-len(test_errors) and report['failures']==report['skipped']==0
            repaired=set();repair_reports=[]
            for folder in (['core0-small','core0-step'] if name=='python-core-0' else ['core1-ablation']):
                r=read(raw/'python-repairs'/folder/'result.json');repair_reports.append(r)
                repaired.update(r['success'])
            assert repaired<=ids
            final=before|repaired;success.update(final)
            final_groups[name]=dict(assigned=len(ids),successful_unique_ids=sorted(final),
                not_successful_ids=sorted(ids-final),initial_class_error_ids=sorted(absent),initial_test_errors=test_errors,
                repair_reports=repair_reports,rerun_successes_already_successful=len(repaired&before),
                no_initial_error_logs_deleted=True)
        elif name in ['python-mean-evidence','python-spatial-evidence']:
            folder='mean-final' if name=='python-mean-evidence' else 'spatial-final'
            latest=read(raw/'python-repairs'/folder/'python-test-inventory.json')
            assert set(latest['discovered_test_ids'])==set(discovered)
            assert latest['assigned_total']==4 and latest['skipped']==0
            if latest['tests_run']==4 and latest['errors']==latest['failures']==0:success.update(ids)
            else:
                assert latest['tests_run']==0 and latest['errors']==1
                not_executed.extend(sorted(ids))
            final_groups[name]=latest
        else:raise ValueError('unexpected test group')
    repairs=final_groups['python-core-0']['repair_reports']+final_groups['python-core-1']['repair_reports']
    latest_nc=[final_groups[n] for n in ('python-mean-evidence','python-spatial-evidence')]
    python=dict(discovered_total=len(discovered),successful_unique_tests=len(success),
        not_executed_ids=sorted(not_executed),not_successful_ids=sorted(set(discovered)-success),
        all_ids_partitioned_exactly_once=True,all_discovered_ids_accounted_for=success|set(not_executed)==set(discovered),
        initial_errors_retained=initial_errors,final_groups=final_groups,
        preparation_errors=sum(r['errors'] for r in latest_nc if r['tests_run']==0),
        failures=sum(len(r['failures']) for r in repairs)+sum(r['failures'] for r in latest_nc),
        errors=sum(len(r['errors']) for r in repairs)+sum(r['errors'] for r in latest_nc),
        skipped=sum(len(r['skipped']) for r in repairs)+sum(r['skipped'] for r in latest_nc),
        passed=len(success)==len(discovered))
    return dict(passed=rust['passed'] and python['passed'],rust=rust,python=python,p6_complete=False,
        numerical_acceptance_is_separate=True,windows_ci='unconfirmed')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('raw',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    result=summarize(a.raw)
    with a.output.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
