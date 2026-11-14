"""Observe all existing RK4 accepted steps without changing its numerical work."""
import argparse,gzip,hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'tools'))
import low_preload_nc_reference as reference
import mcc_rate_fem as rate

def capture(out,steps):
    history=out.parent/(out.name+'-history.jsonl.gz');assert not history.exists()
    mesh,ground,case=reference.expected_inputs();full=rate.seed(mesh,ground);frozen_contact=None
    original=reference.advance_rk4;segment=0;rows=0;counts={}
    labels=[('initial_equilibrated',0.,1.),('load',0.,1.),('unload',0.,1.),('excavate_upper',0.,.875),('excavate_upper',.875,.9375),('excavate_upper',.9375,1.)]
    with history.open('xb') as raw,gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as compressed:
        def integrator(system,state,displacement,force_rate,count,progress=None,accepted_callback=None):
            nonlocal segment,rows,full,frozen_contact
            label,a,b=labels[segment];segment+=1
            def observe(t,material,u):
                nonlocal rows,full,frozen_contact
                if accepted_callback:accepted_callback(t,material,u)
                full[system.active]=material.reshape(full[system.active].shape)
                contact=system.contacts(u)
                if frozen_contact is None:frozen_contact=contact
                else:
                    for key,value in contact.items():
                        old=np.array(frozen_contact[key]);old[system.live]=np.array(value)[system.live];frozen_contact[key]=old.tolist()
                row=dict(stage=label,progress=a+(b-a)*t,active_elements=system.active.tolist(),material_state=full.tolist(),displacement_m=u.reshape(-1,3).tolist(),contacts=frozen_contact)
                compressed.write((json.dumps(row,allow_nan=False,separators=(',',':'))+'\n').encode());rows+=1;counts[label]=counts.get(label,0)+1
            return original(system,state,displacement,force_rate,count,progress=progress,accepted_callback=observe)
        reference.advance_rk4=integrator
        try:report=reference.solve(out,steps,'rk4',180.)
        finally:reference.advance_rk4=original
    result=dict(history=str(history.name),sha256=hashlib.sha256(history.read_bytes()).hexdigest(),accepted_rows=rows,by_stage=counts,numerical_integrator_unchanged=True,passed=False)
    (out/'history-binding.json').write_text(json.dumps(result,indent=2)+'\n');return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--steps',type=int,required=True);a=p.parse_args();print(json.dumps(capture(a.output.resolve(),a.steps),indent=2))
