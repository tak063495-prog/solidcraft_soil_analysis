"""Final unique-ID outcomes, retaining every original and retry attempt."""
import argparse
import importlib.util
import re
from pathlib import Path
from common import HERE,read,write

def run(rust,python,retry,output,scratch):
    spec=importlib.util.spec_from_file_location('b02_count_aggregator',HERE/'validation_report.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    module.report(rust,python,python,scratch)
    result=read(scratch);statuses={}
    for name,group in result['python_groups'].items():
        for tid in group['not_executed_test_ids']:statuses[tid]='not_executed'
        for tid in group['executed_test_ids']:statuses[tid]='passed'
        text=(python/name/'python-tests.log').read_text()
        for label,method,klass in re.findall(r'(?m)^(FAIL|ERROR): (test_[^\s(]+) \(([^)]+)\)',text):
            tid=klass if klass.endswith('.'+method) else klass+'.'+method
            assert tid in statuses
            statuses[tid]='failed' if label=='FAIL' else 'error'
        for method,klass in re.findall(r'(?m)^(test_[^\s(]+) \(([^)]+)\) \.\.\. skipped',text):
            tid=klass if klass.endswith('.'+method) else klass+'.'+method
            assert tid in statuses;statuses[tid]='skipped'
    assert len(statuses)==result['python_discovered']
    initial=statuses.copy();recovery=read(retry/'retry.json')
    assert set(recovery['outcomes'])==set(recovery['executed_test_ids'])
    for tid,outcome in recovery['outcomes'].items():
        assert tid in statuses;statuses[tid]=outcome
    result.update(passed=False,python_passed=sum(v=='passed' for v in statuses.values()),python_failures=sum(v=='failed' for v in statuses.values()),
                  python_body_errors=sum(v=='error' for v in statuses.values()),python_skipped=sum(v=='skipped' for v in statuses.values()),
                  python_not_executed=sum(v=='not_executed' for v in statuses.values()),python_executed=sum(v!='not_executed' for v in statuses.values()),
                  python_final_failed_ids=sorted(k for k,v in statuses.items() if v in ('failed','error')),
                  python_not_executed_ids=sorted(k for k,v in statuses.items() if v=='not_executed'),
                  python_executed_ids=sorted(k for k,v in statuses.items() if v!='not_executed'),
                  python_final_outcomes=statuses,python_initial_outcomes=initial,python_recovery=recovery,
                  python_initial_skipped=sum(v=='skipped' for v in initial.values()),
                  python_method_attempts_including_recovery=sum(g['tests_run'] for g in result['python_groups'].values())+recovery['tests_run'],
                  python_errors=sum(v=='error' for v in statuses.values())+result['python_setup_errors'],
                  outcome_policy='Final unique original IDs; affected bodies rerun after exact fixture recovery / real PowerShell installation. First errors/skip and all process exits retained. Baseline ledger failures remain failures.',
                  baseline_audit='validation/baseline/inherited-audit-result.json',
                  historical_nc_records_preserved=True,production_sources_changed=False,
                  source_snapshot_timing='875 sources fixed before all Python groups; after fmt/clippy/build and during Rust compilation. Present original Git bytes and source fingerprints checked again at completion.')
    assert sum(result[k] for k in ('python_passed','python_failures','python_body_errors','python_skipped','python_not_executed'))==result['python_discovered']
    result['passed']=bool(all(v['exit_code']==0 for k,v in result['processes'].items() if k in ('fmt','clippy','build','rust-test','rust-doc')) and result['rust_failed']==0 and result['rust_ignored']==0 and result['rust_filtered']==0 and result['python_failures']==0 and result['python_body_errors']==0 and result['python_skipped']==0 and result['python_not_executed']==0)
    write(output,result)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--rust',type=Path,required=True);p.add_argument('--python',type=Path,required=True);p.add_argument('--retry',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--scratch',type=Path,required=True);a=p.parse_args();run(a.rust,a.python,a.retry,a.output,a.scratch)
