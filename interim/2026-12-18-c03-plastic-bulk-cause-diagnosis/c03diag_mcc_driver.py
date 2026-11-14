"""Fresh input-only 3D MCC reference for a CLI-accepted case, all history."""
import os
for _key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[_key]='1'
import argparse
from pathlib import Path
import json
import time
import numpy as np
from c03diag_study import read,write
from c03diag_mcc_reference import Model,run,STAGES,CONTACT


def main():
    p=argparse.ArgumentParser();p.add_argument('--case',type=Path,required=True)
    p.add_argument('--mesh',type=Path,required=True);p.add_argument('--steps',type=int,default=32)
    a=p.parse_args();folder=a.case.resolve();case=read(folder/'stages.json')
    model=Model(read(a.mesh),read(folder/'ground.json'),case)
    out=folder/'independent';out.mkdir();start=time.monotonic()
    try:
        with (out/'fields.jsonl').open('x',encoding='utf-8') as f:
            def sink(row,fields):
                f.write(json.dumps(dict(row=row,fields={k:v.tolist() for k,v in fields.items()}),allow_nan=False)+'\n');f.flush()
                print(row['stage'],row['progress'],flush=True)
            result=run(model,a.steps,case['settings']['material_maximum_substep'],case['joint_system']['maximum_increment_m'],sink)
        write(out/'result.json',result)
        write(out/'process.json',dict(pid=os.getpid(),elapsed_seconds=time.monotonic()-start,input_only=True))
        # Native fields are extracted only after the fresh reference has finished.
        with (folder/'native-fields.jsonl').open('x',encoding='utf-8') as f:
            for i,line in enumerate((folder/'raw/accepted.jsonl').open()):
                state=json.loads(line);j=state['joints']['points'];v=lambda k:np.asarray([p[k] for p in j])
                fields=model.fields(np.asarray(state['displacement_m']).ravel(),state['material_states'],v('jump_m'),v('traction_kpa'),v('plastic_slip_m'),v('accumulated_plastic_slip_m'),v('plastic_dissipation_kpa_m'),v('numerical_dissipation_kpa_m'),v('endpoint_work_kpa_m'),np.asarray(state['applied_nodal_load_kn']),[CONTACT.index(p['status']) for p in j])
                fields['physical_pressure']=np.asarray(state['joint_water']['physical_pressure_kpa']);fields['water_load']=np.asarray(state['joint_water']['nodal_load_kn']).reshape(-1,3)
                row=dict(stage=STAGES[i//a.steps],progress=(i%a.steps+1)/a.steps,factor=state['boundary_factor'])
                f.write(json.dumps(dict(row=row,fields={k:v.tolist() for k,v in fields.items()}),allow_nan=False)+'\n')
            assert i+1==7*a.steps
    except Exception as exc:
        write(out/'failure.json',dict(passed=False,complete=False,error=repr(exc)));raise

if __name__=='__main__':main()
