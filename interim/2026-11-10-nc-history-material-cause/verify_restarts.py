"""Compare whole Runs through their numerical stop, using actual sealed checkpoints."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from probe_low_preload_nc import deterministic, sha


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def verify(data, binary, out):
    out.mkdir(parents=True, exist_ok=False)
    identity = sha(binary)
    direct = {}
    for n in (128,256):
        trace, reference = [read(data/f'{kind}{n}/run.json') for kind in ('trace','direct')]
        if deterministic(trace) != deterministic(reference):
            raise ValueError('continuous run differs from every-increment checkpoint trace')
        direct[str(n)] = dict(entire_run_exact_except_elapsed_seconds=True,
                             solver_binary_sha256=trace['provenance']['solver_binary_sha256'])
        if direct[str(n)]['solver_binary_sha256'] != identity:
            raise ValueError('binary changed after trace')
    def resume(item):
        n, stage = item
        source = data/f'trace{n}'/f'checkpoint-stage-{stage}.json'
        cp = read(source); dest = out/f'n{n}-stage{stage}'
        command = [str(binary), 'resume', str(data/f'input{n}/mesh.json'), str(data/f'input{n}/ground.json'), str(source), str(dest)]
        t=time.monotonic()
        with (out/f'n{n}-stage{stage}.log').open('x') as log:
            child = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
        if child.returncode: raise ValueError('separate checkpoint process failed')
        restored = read(dest/'run.json'); expected = read(data/f'direct{n}/run.json')
        if deterministic(restored) != deterministic(expected):
            raise ValueError('separate checkpoint process differs from continuous run')
        return dict(case=n, starting_stage=stage, checkpoint_sha256=sha(source),
            starting_progress=cp['run']['last_committed']['stage_progress'], command=command,
            elapsed_seconds=time.monotonic()-t, exit_code=child.returncode,
            entire_run_exact_except_elapsed_seconds=True, terminal_progress=restored['last_committed']['stage_progress'],
            full_removal_complete=restored['complete'], completion_restart_verified=False)
    with ThreadPoolExecutor(max_workers=2) as pool:
        results=list(pool.map(resume, [(n,s) for n in (128,256) for s in (2,3)]))
    if sha(binary) != identity: raise ValueError('binary changed during restart verification')
    result=dict(continuous_vs_serialized_trace=direct, separate_process_restarts=results,
        solver_binary_sha256=identity, excluded_recursive_fields=['elapsed_seconds'],
        numerical_acceptance=False, completion_restart_verified=False)
    (out/'verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('data','binary','output'): p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();verify(a.data.resolve(),a.binary.resolve(),a.output.resolve())
