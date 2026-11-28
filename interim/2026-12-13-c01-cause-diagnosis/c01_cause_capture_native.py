"""Observe unmodified native solver partitions, bind every trial to saved u bytes."""
import argparse
from collections import Counter
from pathlib import Path
import subprocess
import time
import numpy as np
import c01_cause_diagnostics as diag


def capture(saved, output, solver, control, steps):
    old,_ = diag.dependencies()
    output.mkdir(parents=True,exist_ok=False)
    case = old.configuration(diag.read(saved/'stages.json'),steps,5e-5)
    case['joint_system']['global_step_control'] = ({'switch_bisections':4} if control=='switch' else {'path_error':{'slip_scale_m':1e-4,'dissipation_scale_kpa_m':.01}})
    case['settings']['minimum_increment'] = min(case['settings']['minimum_increment'],1/(steps*(32 if control=='switch' else 65536)))
    diag.write(output/'case.json',case)
    command = [str(solver.resolve()),str(saved/'mesh-h0.json'),str(saved/'ground.json'),str(output/'case.json'),str(output/'raw')]
    started = time.monotonic()
    with (output/'process.log').open('x') as log:
        result = subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,timeout=1800)
    if result.returncode:
        raise ValueError('native diagnostic solve failed: '+str(result.returncode))
    trials,accepted = {},[]
    with (output/'raw/trials.jsonl').open() as stream:
        for line in stream:
            row = __import__('json').loads(line)
            if row['kind']=='accepted':
                accepted.append((row['stage_index'],row['progress']))
            else:
                key = (row['stage_index'],row['progress'],np.asarray(row['old_u'],dtype=np.float64).tobytes(),np.asarray(row['u'],dtype=np.float64).tobytes())
                trials.setdefault(key,set()).add(tuple(row['contact_partition_counts']))
    before = np.zeros((54,3),dtype=np.float64)
    rows,hist = [],Counter()
    for row,fields in old.states(saved/f'g{steps}/native'):
        stage = diag.STAGES.index(row['stage'])
        key = (stage,row['progress'],before.tobytes(),fields['u'].tobytes())
        found = trials.get(key,set())
        if len(found)!=1:
            raise ValueError('saved accepted displacement not uniquely matched: '+str((row['stage'],row['progress'],len(found))))
        counts = list(next(iter(found)))
        rows.append(dict(stage=row['stage'],progress=row['progress'],contact_partition_counts=counts))
        hist.update(counts)
        before = fields['u']
    if accepted != [(diag.STAGES.index(r['stage']),r['progress']) for r in rows]:
        raise ValueError('native accepted coordinate sequence changed')
    diag.write(output/'counts.json',dict(rows=rows,counts_bound_to_saved_old_and_new_displacement_bytes=True,accepted_coordinates_exact=True,partition_histogram=dict(sorted(hist.items())),control=control,steps=steps,solver_sha256=diag.sha(solver),elapsed_seconds=time.monotonic()-started))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--saved',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--solver',type=Path,required=True)
    p.add_argument('--control',choices=('switch','path'),required=True)
    p.add_argument('--steps',type=int,choices=(16,32,64,128,256),required=True)
    a=p.parse_args()
    capture(a.saved.resolve(),a.output.resolve(),a.solver,a.control,a.steps)


if __name__=='__main__':
    main()
