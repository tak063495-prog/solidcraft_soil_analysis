"""Verify unmodified Python roundtrip/resealing preserves the original Run hash."""
import argparse
import os
from pathlib import Path
from study import read,write,invoke,sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);a=p.parse_args()
    root=a.run.resolve();out=root/'roundtrip-controls';out.mkdir()
    exe=root/'bin'/('rehash.exe' if os.name=='nt' else 'rehash');rows=[]
    for name in ('transfer_force','stress_free','stress_free-cycle'):
        folder=root/name;case=read(folder/'case.json');steps=round(1/case['settings']['initial_increment'])
        for count in (steps+steps//2,4*steps+steps//2):
            original=folder/f'native/checkpoint-{count:04}.json';saved=read(original)
            label=name+'-'+str(count);raw=out/(label+'-raw.json');sealed=out/(label+'-rehashed.json')
            write(raw,saved);result=invoke([exe,raw,sealed],out/(label+'.log'));assert result['exit_code']==0
            restored=read(sealed)
            assert restored['run_sha256']==saved['run_sha256']
            assert restored['identity']==saved['identity']
            rows.append(dict(case=name,accepted_number=count,original_checkpoint_sha256=sha(original),
                canonical_run_hash_identical=True,identity_identical=True,process=result))
    write(root/'roundtrip-controls.json',dict(passed=True,hash_control_not_new_solve=True,checks=rows))


if __name__=='__main__':main()
