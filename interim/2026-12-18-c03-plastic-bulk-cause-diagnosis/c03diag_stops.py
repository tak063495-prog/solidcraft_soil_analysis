"""Classify actual rejected attempts; a global residual has no single failing GP."""
import argparse
from collections import Counter
from pathlib import Path
import json
import numpy as np
from c03diag_study import read,write,Model


def main():
    p=argparse.ArgumentParser();p.add_argument('--case',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run=read(a.case/'raw/run.json');base=a.case.parents[1]/'matrix'
    model=Model(read(base/'mesh.json'),read(a.case/'ground.json'),read(a.case/'stages.json'))
    progress=[json.loads(s)['event'] for s in (a.case/'raw/progress.jsonl').open()]
    # Diagnostic iterate events are already evaluated; no replay is counted as equilibrium.
    grouped={};last=None;errors=[]
    for line in (a.case/'raw/cause-events.jsonl').open():
        event=json.loads(line)
        if event['kind']=='iterate':
            last=(event['stage'],event['progress']);grouped.setdefault(last,[]).append(event)
        if 'error' in event:errors.append(dict(attempt=last,error=event['error']))
    rows=[]
    rejected=[x for x in run['attempts'] if not x['accepted']]
    for attempt in rejected:
        index=next(i+1 for i,s in enumerate(run['case']['stages']) if s['name']==attempt['stage'])
        # Retries at the same endpoint have distinct old committed progress.
        candidates=[v for v in grouped[index,attempt['progress_to']]
                    if v['committed']['stage_progress']==attempt['progress_from']]
        event=candidates[-1];states=event['states'];strain=np.asarray(event['strain_increment'])
        residual=np.asarray(event['reduced_residual']);point=np.max(np.abs(strain),axis=1)
        jump=np.einsum('qag,g->qa',model.system.c,np.asarray(event['displacement_m']).ravel())
        gap=np.einsum('qa,qa->q',jump,model.system.normal)
        reason=attempt['reason'];classification='numerical_global_nonconvergence'
        if 'invalid input' in reason:classification='scope_rejection'
        rows.append(dict(stage=attempt['stage'],progress_from=attempt['progress_from'],progress_to=attempt['progress_to'],
             classification=classification,reason=reason,relative_force_residual=event['relative_force_residual'],
             force_gate=run['case']['settings']['relative_force_tolerance'],iteration=event['iteration'],
             rejected_global_residual_not_single_gp=True,failing_gp=None,
             maximum_reduced_residual_dof=int(np.argmax(np.abs(residual))),
             maximum_reduced_residual=float(np.max(np.abs(residual))),
             material_gp_last_returns={str(i):s.get('last_return') for i,s in enumerate(states)},
             maximum_strain_increment_gp=int(np.argmax(point)),maximum_strain_increment=float(point.max()),
             material_partitions=event['partitions'][:len(states)],contact_partitions=event['partitions'][len(states):],
             contact_gp_jump_m=jump.tolist(),contact_gp_gap_m=gap.tolist(),
             physical_failure_proven=False))
    write(a.output,dict(complete=run['complete'],failure_reason=run.get('failure_reason'),
           accepted=sum(e['kind']=='accepted' for e in progress),rejected=rows,
           local_evaluation_error_count=len(errors),local_evaluation_errors=errors,
           input_scope_rejection=False,physical_failure_proven=False,passed=False,
           note='An equilibrium rejection is assessed on the reduced nodal residual, not a GP failure criterion.'))

if __name__=='__main__':main()
