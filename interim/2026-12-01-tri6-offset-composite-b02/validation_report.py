"""Final, nonduplicated Rust/Python counts and unexecuted IDs; no skip fixes."""
import argparse
from pathlib import Path
import re
from cases import read,write


def report(rust,python,retry,output):
    process=read(rust/'processes.json');rust_groups=[]
    for name in ('rust-test','rust-doc'):
        text=(rust/(name+'.log')).read_text()
        for values in re.findall(r'test result: (\w+)\. (\d+) passed; (\d+) failed; (\d+) ignored; (\d+) measured; (\d+) filtered out',text):
            rust_groups.append(dict(result=values[0],passed=int(values[1]),failed=int(values[2]),ignored=int(values[3]),measured=int(values[4]),filtered_out=int(values[5]),command=name))
    assert rust_groups
    groups={};inventory=None;assigned=[];executed=[]
    for name in ('python-core-0','python-core-1','python-mean-evidence','python-spatial-evidence'):
        folder=(retry if name=='python-core-0' else python)
        data=read(folder/name/'python-test-inventory.json')
        ids=set(data['discovered_test_ids'])
        if inventory is None:inventory=ids
        assert inventory==ids
        own=data['all_group_test_ids'][name];assigned+=own
        text=(folder/name/'python-tests.log').read_text()
        # Python 3.13 includes the method in the parenthesized full ID;
        # older unittest versions show only module.Class there.
        ran=[klass if klass.endswith('.'+method) else f'{klass}.{method}'
             for method,klass in re.findall(r'(?m)^(test_[^\s(]+) \(([^)]+)\) \.\.\.',text)]
        assert len(ran)==data['tests_run'] and len(ran)==len(set(ran)),(name,len(ran),data['tests_run'])
        assert set(ran)<=set(own);executed+=ran
        group=dict(data,executed_test_ids=ran,not_executed_test_ids=sorted(set(own)-set(ran)),
                   body_errors=len(re.findall(r'(?m)^ERROR: test_',text)),
                   setup_errors=len(re.findall(r'(?m)^ERROR: setUp',text)))
        groups[name]=group
        process[name]=read(folder/(name+'.process.json'))
    assert len(assigned)==len(set(assigned)) and set(assigned)==inventory and len(executed)==len(set(executed))
    skipped=sum(x['skipped'] for x in groups.values());failures=sum(x['failures'] for x in groups.values());errors=sum(x['errors'] for x in groups.values())
    expected=sum(x['expected_failures'] for x in groups.values());unexpected=sum(x['unexpected_successes'] for x in groups.values())
    total=sum(g['tests_run'] for g in groups.values())
    ok=all(x['exit_code']==0 for x in process.values()) and not(skipped or failures or errors or expected or unexpected) and total==len(inventory)
    result=dict(passed=ok,all_requested_commands_completed=True,processes=process,rust_groups=rust_groups,
                rust_passed=sum(x['passed'] for x in rust_groups),rust_failed=sum(x['failed'] for x in rust_groups),
                rust_ignored=sum(x['ignored'] for x in rust_groups),rust_filtered=sum(x['filtered_out'] for x in rust_groups),
                python_discovered=len(inventory),python_executed=total,python_passed=total-failures-sum(x['body_errors'] for x in groups.values())-skipped-expected-unexpected,
                python_failures=failures,python_errors=errors,python_skipped=skipped,python_expected_failures=expected,python_unexpected_successes=unexpected,
                python_body_errors=sum(x['body_errors'] for x in groups.values()),python_setup_errors=sum(x['setup_errors'] for x in groups.values()),
                python_not_executed=len(inventory-set(executed)),python_not_executed_ids=sorted(inventory-set(executed)),
                python_all_ids_assigned_exactly_once=True,python_executed_ids=executed,python_groups=groups,
                retained_first_core0=read(python/'python-core-0/python-test-inventory.json'),
                historical_nc=dict(original_archives_and_exact_checks_unchanged=True,allow_platform_float_differences_Linux=False,past_unexecuted_records_preserved=True,new_NC_numerical_acceptance=False),
                first_core0_note='Started before all artifact downloads completed; 300 executed, one missing ZIP error. Retained, final retry whole 300 used without duplicate counts.',
                windows_ci='New pushed commit not confirmed; old success, Linux, timeout and cancellation are not transferred')
    write(output,result)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--rust',type=Path,required=True);p.add_argument('--python',type=Path,required=True);p.add_argument('--retry',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();report(a.rust,a.python,a.retry,a.output)
