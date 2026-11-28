"""Fresh independent global equilibrium with explicit minimum contact partitions."""
import argparse
from pathlib import Path
from a02cause_binding import study,reference,fields,OLD
from a02cause_driver import CASES,write
from a02cause_reference import controlled_run
from a02cause_metrics import compare


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--minimum-parts',type=int,default=2);a=p.parse_args()
    run=controlled_run(a.minimum_parts);a.output.mkdir(parents=True,exist_ok=False)
    write(a.output/'source-binding.json',{str(p):study.sha(p) for p in [Path(__file__),Path(__file__).with_name('a02cause_reference.py'),OLD/'reference.py']})
    result={}
    for name,policy,cycle in CASES:
        folder=a.output/name;study.make_inputs(folder,study.fixture(policy,4,cycle))
        m=reference.Model(*study.fixture(policy,4,cycle))
        with (folder/'independent-fields.jsonl').open('x',encoding='utf-8',newline='\n') as stream:
            r=run(m,4,fields.Sink(stream))
        write(folder/'independent-result.json',r)
        assert min(n for row in r['rows'] for n in row['contact_partitions'])>=a.minimum_parts
        result[name]=dict(minimum_parts=a.minimum_parts,input_only_fresh_global_equilibrium=True,
            maximum_force_residual=r['maximum_force_residual'],comparison=compare(
            list(fields.load_fields(a.run/'n4'/name/'independent-fields.jsonl')),
            list(fields.load_fields(folder/'independent-fields.jsonl')),fields.CATEGORIES))
    write(a.output/'partition-global-controls.json',result)


if __name__=='__main__':main()
