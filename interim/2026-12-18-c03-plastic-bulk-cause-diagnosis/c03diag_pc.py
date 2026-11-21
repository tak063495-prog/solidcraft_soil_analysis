"""Explicit PC-only procedure for larger Tet10/history; never run by CI.

No constitutive/control changes. Diagnostic Newton copies are disabled beyond
the existing 512-GP bound. Dense independent reference can be prohibitively slow.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[_key]='1'
import argparse
from pathlib import Path
import json
import time
import sys
import numpy as np
from c03diag_study import initialize,read,write,sha,binary,configuration,invoke,schedule,restart
from c03diag_reference import Model,run,STAGES,CONTACT
from c03diag_metrics import load,compare
from refine_tet10 import refine


def export_native(out,folder,steps):
    model=Model(read(out/'mesh.json'),read(out/'ground.json'),read(folder/'stages.json'),user_pc=True)
    with (folder/'fields.jsonl').open('x',encoding='utf-8') as f:
        for i,line in enumerate((folder/'raw/accepted.jsonl').open()):
            state=json.loads(line);j=state['joints']['points'];v=lambda k:np.asarray([p[k] for p in j])
            fields=model.fields(np.asarray(state['displacement_m']).ravel(),state['material_states'],v('jump_m'),v('traction_kpa'),v('plastic_slip_m'),v('accumulated_plastic_slip_m'),v('plastic_dissipation_kpa_m'),v('numerical_dissipation_kpa_m'),v('endpoint_work_kpa_m'),np.asarray(state['applied_nodal_load_kn']),[CONTACT.index(p['status']) for p in j])
            fields['physical_pressure']=np.asarray(state['joint_water']['physical_pressure_kpa']);fields['water_load']=np.asarray(state['joint_water']['nodal_load_kn']).reshape(-1,3)
            row=dict(stage=STAGES[i//steps],progress=(i%steps+1)/steps,factor=state['boundary_factor'])
            f.write(json.dumps(dict(row=row,fields={k:v.tolist() for k,v in fields.items()}),allow_nan=False)+'\n')
        assert i+1==7*steps


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--user-pc',action='store_true',required=True)
    p.add_argument('--mode',choices=['init','native','independent','compare','restart','independent-worker'],required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--solver',type=Path)
    p.add_argument('--refine',type=int,choices=range(4),default=0)
    p.add_argument('--folder',default='baseline');p.add_argument('--other',default='independent')
    p.add_argument('--steps',type=int,default=32);p.add_argument('--material',type=float,default=1e-4)
    p.add_argument('--contact',type=float,default=1e-4);p.add_argument('--timeout',type=int,default=14400)
    a=p.parse_args();os.environ['C03_USER_PC']='1';out=a.output.resolve()
    try:
        if a.mode=='init':
            initialize(out,a.solver.resolve(),'tet10',30.,open_factor=.12)
            if a.refine:
                (out/'mesh.json').rename(out/'mesh-seed.json');mesh=read(out/'mesh-seed.json')
                for level in range(a.refine):
                    mesh,audit=refine(mesh);write(out/f'refinement-{level+1}.json',audit)
                write(out/'mesh.json',mesh)
                identity=read(out/'identity.json');(out/'identity.json').rename(out/'initial-identity.json')
                identity['inputs']['mesh.json']=sha(out/'mesh.json');identity['refinement_levels']=a.refine
                write(out/'identity.json',identity)
        elif a.mode=='compare':write(out/(a.folder+'--'+a.other+'-comparison.json'),compare(load(out/a.folder/'fields.jsonl'),load(out/a.other/'fields.jsonl')))
        elif a.mode=='restart':restart(out,out/a.folder,out/(a.folder+'-restarts'),a.timeout,a.steps)
        elif a.mode=='independent':
            command=sys.argv[1:].copy();command[command.index('--mode')+1]='independent-worker'
            process=invoke([sys.executable,Path(__file__),*command],out/(a.folder+'-driver.log'),a.timeout)
            write(out/(a.folder+'-terminal.json'),process)
            if process['exit_code']:raise RuntimeError('independent incomplete; retain partial evidence')
        else:
            folder=out/a.folder;folder.mkdir();case=configuration(out,a.steps,a.material,a.contact);write(folder/'stages.json',case)
            if a.mode=='native':
                process=invoke([binary(out),out/'mesh.json',out/'ground.json',folder/'stages.json',folder/'raw'],folder/'solve.log',a.timeout)
                write(folder/'process.json',process)
                if process['exit_code']:raise RuntimeError('native incomplete; retain all partial files')
                write(folder/'schedule.json',schedule(read(folder/'raw/run.json'),a.steps));export_native(out,folder,a.steps)
            else:
                model=Model(read(out/'mesh.json'),read(out/'ground.json'),case,user_pc=True);start=time.monotonic()
                with (folder/'fields.jsonl').open('x',encoding='utf-8') as f:
                    def sink(row,fields):
                        f.write(json.dumps(dict(row=row,fields={k:v.tolist() for k,v in fields.items()}),allow_nan=False)+'\n');f.flush()
                    result=run(model,a.steps,a.material,a.contact,sink)
                write(folder/'result.json',result);write(folder/'process.json',dict(pid=os.getpid(),elapsed_seconds=time.monotonic()-start,input_only=True))
    except Exception as exc:
        if out.exists():write(out/f'failure-pc-{time.time_ns()}.json',dict(passed=False,complete=False,error=repr(exc)))
        raise

if __name__=='__main__':main()
