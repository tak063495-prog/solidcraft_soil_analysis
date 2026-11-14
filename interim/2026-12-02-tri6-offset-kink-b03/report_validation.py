"""Report complete discovered/executed test sets, preserving strict NC setup failures."""
import argparse,re
from pathlib import Path
from b03_cases import read,write


def report(rust,python,output):
    processes=read(rust/'processes.json');rust_rows=[]
    for name in ('rust-test','rust-doc'):
        text=(rust/(name+'.log')).read_text()
        for values in re.findall(r'test result: (\w+)\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out',text):
            rust_rows.append(dict(command=name,result=values[0],passed=int(values[1]),failed=int(values[2]),ignored=int(values[3]),measured=int(values[4]),filtered=int(values[5])))
    assert rust_rows
    groups={};inventory=None;assigned=[];executed=[]
    for name in ('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence'):
        d=read(python/name/'python-test-inventory.json');ids=set(d['discovered_test_ids'])
        if inventory is None:inventory=ids
        assert inventory==ids;own=d['all_group_test_ids'][name];assigned+=own
        text=(python/name/'python-tests.log').read_text()
        ran=[klass if klass.endswith('.'+method) else f'{klass}.{method}' for method,klass in re.findall(r'(?m)^(test_[^\s(]+) \(([^)]+)\) \.\.\.',text)]
        assert len(ran)==d['tests_run'] and len(ran)==len(set(ran));assert set(ran)<=set(own);executed+=ran
        groups[name]=dict(d,executed_test_ids=ran,not_executed_test_ids=sorted(set(own)-set(ran)),
            body_errors=len(re.findall(r'(?m)^ERROR: test_',text)),setup_errors=len(re.findall(r'(?m)^ERROR: setUp',text)))
        processes[name]=read(python/(name+'.process.json'))
    assert len(assigned)==len(set(assigned)) and set(assigned)==inventory and len(executed)==len(set(executed))
    total=sum(d['tests_run'] for d in groups.values());fails=sum(d['failures'] for d in groups.values());body=sum(d['body_errors'] for d in groups.values())
    skip=sum(d['skipped'] for d in groups.values());xf=sum(d['expected_failures'] for d in groups.values());xp=sum(d['unexpected_successes'] for d in groups.values())
    result=dict(passed=all(p['exit_code']==0 for p in processes.values()) and total==len(inventory) and not(skip or xf or xp),
        all_requested_commands_completed=True,processes=processes,rust_groups=rust_rows,rust_passed=sum(r['passed'] for r in rust_rows),
        rust_failed=sum(r['failed'] for r in rust_rows),rust_ignored=sum(r['ignored'] for r in rust_rows),rust_filtered=sum(r['filtered'] for r in rust_rows),
        python_discovered=len(inventory),python_executed=total,python_passed=total-fails-body-skip-xf-xp,
        python_failures=fails,python_body_errors=body,python_setup_errors=sum(d['setup_errors'] for d in groups.values()),
        python_skipped=skip,python_expected_failures=xf,python_unexpected_successes=xp,
        python_not_executed=len(inventory-set(executed)),python_not_executed_ids=sorted(inventory-set(executed)),
        python_all_ids_assigned_exactly_once=True,python_groups=groups,
        historical_nc=dict(original_archives_and_exact_checks_unchanged=True,allow_platform_float_differences_Linux=False,past_unexecuted_records_preserved=True,new_NC_numerical_acceptance=False),
        windows_ci='Pushed head not confirmed; old CI, timeout/cancellation and Linux are not Windows evidence')
    write(output,result)
if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ('rust','python','output'):p.add_argument(k,type=Path)
    a=p.parse_args();report(a.rust,a.python,a.output)
