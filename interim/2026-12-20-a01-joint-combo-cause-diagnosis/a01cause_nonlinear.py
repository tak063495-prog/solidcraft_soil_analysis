"""Postprocess complete raw Newton observations and old local-work audits only."""
import argparse
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
from a01cause_binding import fixtures,legacy,OLD,original_reference
from a01cause_driver import write


def old_local_auditor():
    saved={key:sys.modules.get(key) for key in ('fixtures','reference','fields')}
    try:
        for key in saved:sys.modules[key]=legacy[key]
        spec=importlib.util.spec_from_file_location('a01cause_legacy_local_auditor',OLD/'diagnose.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    finally:
        for key,value in saved.items():
            if value is None:sys.modules.pop(key,None)
            else:sys.modules[key]=value
    return module


def main():
    p=argparse.ArgumentParser();p.add_argument('--formal',type=Path,required=True);p.add_argument('--old',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    folder=a.formal/'stacked';model=original_reference.Model(*(fixtures.read(folder/(k+'.json')) for k in ('mesh','ground','case')));s=model.system
    rows=[]
    for line in (folder/'native/trace.jsonl').read_text(encoding='utf-8').splitlines():
        event=json.loads(line)
        if event['kind']!='iterate':continue
        u=np.asarray(event['displacement_m']).ravel();jump=s.c@u;old=event['committed'];points=old['joints']['points']
        rest=np.asarray([v.get('birth',{}).get('elastic_offset_m',[0.,0.,0.]) for v in points])
        gap=np.sum((jump-rest)*s.normal,axis=1)
        live=all(v['active'] for v in old['joints']['lifecycle']);support=None
        if live:
            from types import SimpleNamespace
            traction=np.zeros_like(jump);counts=np.asarray(event['partitions'][len(event['states']):])
            before=np.asarray([v['jump_m'] for v in points]);oldq=np.asarray([v['traction_kpa'] for v in points])
            for count in np.unique(counts):
                mask=counts==count;sub=SimpleNamespace(**{k:getattr(s,k)[mask] for k in ('normal','kn','kt','mu')})
                traction[mask]=original_reference.contact(sub,before[mask]-rest[mask],oldq[mask],jump[mask]-rest[mask],int(count))[0][:,0]
            state=dict(active_elements=event['active'],displacement_m=event['displacement_m'],material_states=event['states'],boundary_factor=old['boundary_factor'],joints=dict(points=[dict(traction_kpa=v.tolist()) for v in traction]))
            support=s.rigid_support(model.case,state)
        rows.append(dict(stage=event['stage'],progress=event['progress'],projected_load=event['projected_load'],iteration=event['iteration'],
            relative_force_residual=event['relative_force_residual'],contact_active=live,open_gp=int(np.sum(gap>0)),closed_gp=int(np.sum(gap<=0)),
            gap_min=float(gap.min()),gap_max=float(gap.max()),independent_rigid_support_diagnostic=support))
    ref=[json.loads(x) for x in (folder/'independent/trace.jsonl').read_text(encoding='utf-8').splitlines()]
    last=next(x for x in reversed(ref) if x['kind']=='linear')
    cuts=[x for x in ref if x['kind']=='line_search' and (x['stage'],x['progress'],x['iteration'])==(last['stage'],last['progress'],last['iteration'])]
    r0=float(np.linalg.norm(last['rhs']))
    write(a.output/'nonlinear-path-diagnosis.json',dict(passed=False,local_replay_is_global_equilibrium=False,
        native_iterates=rows,last_independent_line_search=dict(stage=last['stage'],progress=last['progress'],iteration=last['iteration'],
            starting_residual_norm=r0,trial_count=len(cuts),minimum_trial_residual_norm=min(v['residual_norm'] for v in cuts),trials=cuts),
        physical_instability_proved=False,supports_or_constraints_changed=False))
    auditor=old_local_auditor();results={}
    for folder in sorted((a.old/'selected').iterdir()):
        if (folder/'comparison.json').exists():results[folder.name]=auditor.diagnose(folder,8)
    write(a.output/'actual-partition-local-audit.json',dict(passed=False,saved_state_audit_only=True,counted_as_fresh_3d_solve=False,cases=results))

if __name__=='__main__':main()
