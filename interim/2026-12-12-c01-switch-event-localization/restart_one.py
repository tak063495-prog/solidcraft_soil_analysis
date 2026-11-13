"""Resume one actual finest-path checkpoint in a different CLI process."""
import argparse
import importlib.util
import itertools
from pathlib import Path
import subprocess
import time
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('c01_path_actual_restart',HERE/'study.py')
st=importlib.util.module_from_spec(spec);spec.loader.exec_module(st)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--steps',type=int,default=256)
    p.add_argument('--timeout',type=int,default=1800);a=p.parse_args()
    out=a.output.resolve();folder=out/f'g{a.steps}';raw=folder/'native-raw'
    checkpoint=max(raw.glob('checkpoint-*.json'),key=lambda p:int(p.stem.split('-')[-1]))
    cp=st.read(checkpoint);phase=cp['run']['last_committed']['stage_index']
    target=folder/('separate-'+checkpoint.stem)
    binary=out/('probe.exe' if (out/'probe.exe').exists() else 'probe')
    assert st.sha(binary)==cp['run']['provenance']['solver_binary_sha256']
    assert cp['run']['case']['joint_system']['global_step_control']['path_error']==st.PATH_ERROR
    cmd=[str(binary),str(out/'mesh-h0.json'),str(out/'ground.json'),str(folder/'native-stages.json'),str(target),'--resume',str(checkpoint)]
    started=time.monotonic()
    with (folder/(target.name+'.log')).open('x') as log:
        result=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=a.timeout)
    assert result.returncode==0
    assert st.without_clock(st.read(target/'run.json'))==st.without_clock(st.read(raw/'run.json'))
    n=int(checkpoint.stem.split('-')[-1]);count=0
    with (raw/'accepted.jsonl').open('rb') as left,(target/'accepted.jsonl').open('rb') as right:
        for _ in range(n): next(left)
        for x,y in itertools.zip_longest(left,right):
            assert x is not None and x==y,'restart committed suffix bytes differ'
            count+=1
    verification=st.read(target/'verification.json');direct=st.read(raw/'verification.json')
    assert verification['resumed_from_file'] and verification['process_id']!=direct['process_id']
    st.write(folder/'restart-verification.json',dict(steps=a.steps,new_option_enabled=True,
        all_seven_phases=False,actual_separate_process=True,uninterrupted_extra_run_performed=False,
        checks=[dict(stage=st.STAGES[phase],checkpoint=checkpoint.name,checkpoint_sha256=st.sha(checkpoint),
            process_id=verification['process_id'],direct_process_id=direct['process_id'],
            suffix_states_exact=count,suffix_state_bytes_exact=True,entire_run_exact_except_elapsed=True,
            elapsed_seconds=time.monotonic()-started,command=cmd)]))
    print('actual separate-process restart exact',count,flush=True)

if __name__=='__main__': main()
