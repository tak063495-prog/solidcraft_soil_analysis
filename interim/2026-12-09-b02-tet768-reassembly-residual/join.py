"""Combine completed fresh cases by copying bytes; never manufactures states."""
import argparse
import shutil
from pathlib import Path
from common import read,write

def join(output,sources):
    output.mkdir(parents=True,exist_ok=False);records=[];bindings=[]
    for index,source in enumerate(sources):
        before=read(source/'source-binding-before.json');after=read(source/'source-binding-after.json');assert before==after
        for case in read(source/'fresh-solves.json')['cases']:
            assert not any(r['name']==case['name'] for r in records)
            shutil.copytree(source/case['name'],output/case['name'])
            records.append(dict(case,source_directory=str(source)))
        binding=output/'source-bindings'/str(index);binding.mkdir(parents=True)
        for name in ('source-binding-before.json','source-binding-after.json'):
            shutil.copy2(source/name,binding/name)
        bindings.append(dict(source=str(source),frozen_files=len(before),before_after_equal=True))
    records.sort(key=lambda r:(r['level'],r['steps']))
    write(output/'fresh-solves.json',dict(passed=False,b02_complete=False,cases=records,source_bindings=bindings,
          native_fresh_histories=len(records),independent_fresh_histories=len(records),
          native_full_accepted_increments=sum(r.get('native_accepted_increments',r['comparison']['accepted_points']) for r in records),
          native_field_observation_pairs=sum(r['comparison']['accepted_points'] for r in records),
          independent_full_accepted_increments=sum(r.get('independent_accepted_increments',r['comparison']['accepted_points']) for r in records),
          supplemental_unobserved_native_histories=sum('plain' in r for r in records),
          note='768/8,16 have five saved native frame observations; 768/4 has all 20 states. No unobserved intermediate field is invented.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);p.add_argument('--sources',nargs='+',type=Path,required=True);a=p.parse_args();join(a.output.resolve(),[p.resolve() for p in a.sources])
