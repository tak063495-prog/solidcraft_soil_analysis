"""Use the checkpoint producer binary; hash-mismatch rejection is not evidence."""
import argparse
import json
import os
from pathlib import Path
import sys
from study import read, write, sha, invoke, input_paths


def check(folder, probe):
    out=folder/'tamper-authenticated';out.mkdir();rows=[]
    original=folder/'tamper'
    requests=sorted(original.glob('*-rehashed.json'))
    assert len(requests)==6
    for sealed in requests:
        saved=read(sealed)
        assert saved['identity']['solver_binary_sha256']==sha(probe),'wrong producer binary'
        mode=sealed.name.removesuffix('-rehashed.json')
        result=invoke([probe,*input_paths(folder),out/(mode+'-resume'),'--resume',sealed],out/(mode+'.log'))
        text=(out/(mode+'.log')).read_text()
        assert result['exit_code']==1 and not result['timed_out']
        assert 'checkpoint accepted-path replay mismatch' in text,text
        assert 'input/solver hash mismatch' not in text,text
        assert not (out/(mode+'-resume/run.json')).exists()
        rows.append(dict(mode=mode,outer_hash_recomputed=True,producer_binary_sha256=sha(probe),
            checkpoint_identity_matches_producer=True,process=result,
            expected_error='checkpoint accepted-path replay mismatch',error=text))
    write(out/'verification.json',dict(passed=True,software_authentication_only=True,checks=rows))
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True)
    p.add_argument('--cases',nargs='+',default=['transfer_force','stress_free','stress_free-cycle']);a=p.parse_args()
    probe=(a.run/'bin'/('probe.exe' if os.name=='nt' else 'probe')).resolve()
    results={name:check((a.run/name).resolve(),probe) for name in a.cases}
    write(a.run/('authenticated-tampering-'+str(os.getpid())+'.json'),results)
    print('authenticated refusals',sum(map(len,results.values())),flush=True)


if __name__=='__main__':main()
