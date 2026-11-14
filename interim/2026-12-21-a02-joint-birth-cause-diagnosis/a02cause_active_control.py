"""Input-only global control with an effective contact subdivision perturbation."""
import argparse
from pathlib import Path
import sys
from a02cause_binding import study,fields
from a02cause_driver import CASES,write
from a02cause_metrics import compare


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    driver=Path(__file__).with_name('a02cause_driver.py').resolve()
    write(a.output/'source-binding.json',{p.name:study.sha(p) for p in [Path(__file__),driver,driver.with_name('a02cause_binding.py'),driver.with_name('a02cause_metrics.py')]})
    result={}
    for name,policy,cycle in CASES:
        folder=a.output/name;study.make_inputs(folder,study.fixture(policy,4,cycle))
        command=[sys.executable,driver,'--output',a.output.resolve(),'--reference',folder.resolve(),'--steps',4,'--maximum',.0005]
        process=study.invoke(command,folder/'independent.log',timeout=900);write(folder/'independent-process.json',process)
        if process['exit_code'] or process['timed_out']:raise RuntimeError(process)
        reference=study.read(folder/'independent-result.json')
        counts=[n for row in reference['rows'] for n in row['contact_partitions']]
        baseline=a.run/'n4'/name/'independent-fields.jsonl'
        result[name]=dict(contact_maximum_m=.0005,maximum_actual_partitions=max(counts),
            states_with_partitions_gt_one=sum(max(row['contact_partitions'])>1 for row in reference['rows']),
            input_only_fresh_global_equilibrium=True,comparison=compare(list(fields.load_fields(baseline)),
            list(fields.load_fields(folder/'independent-fields.jsonl')),fields.CATEGORIES),
            maximum_force_residual=reference['maximum_force_residual'])
    assert result['stress_free-cycle']['states_with_partitions_gt_one']>0
    write(a.output/'active-contact-controls.json',result)


if __name__=='__main__':main()
