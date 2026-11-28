"""Actual CLI scope checks, including retained open=1.0 numerical stops."""
from pathlib import Path
import argparse
import copy
import sys
from c03diag_study import read,write,invoke,configuration,sha


def main():
    p=argparse.ArgumentParser();p.add_argument('--matrix',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--cli',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True);base=a.matrix.resolve();out=a.output.resolve();rows=[]
    for label,pc,open_factor,scope in [('mc-open1',None,1.,'standard'),('mcc-pc200',200.,.12,'standard'),
                                      ('mcc-pc105',105.,.12,'standard'),('mcc-mixed',200.,.12,'mixed')]:
        folder=out/label;folder.mkdir();ground=read(base/'ground.json');case=configuration(base,32,1e-4,1e-4)
        next(s for s in case['stages'] if s['name']=='open')['target_factor']=open_factor
        if pc is not None:
            for material in ground['materials']:
                del material['mohr_coulomb']
                material['mcc']=dict(critical_state_m=1.2,**{'lambda':.2},kappa=.04,poisson=.3,
                                     reference_pressure_kpa=1.,specific_volume_ncl_at_reference=3.5)
            for seed in ground['initialization']['geostatic']['seed']['states'].values():
                for key in ('equivalent_plastic_strain','plastic_work_kpa','last_return'):del seed[key]
                seed['preconsolidation_pressure_kpa']=pc
        if scope=='mixed':case['schema']='soil-fem-mixed-mcc-stages/1'
        write(folder/'ground.json',ground);write(folder/'stages.json',case)
        command=[a.cli.resolve(),'solve-stages',base/'mesh.json',folder/'ground.json',folder/'stages.json',folder/'cli-run.json']
        process=invoke(command,folder/'cli.log',900)
        row=dict(label=label,process=process,input_sha256={p.name:sha(p) for p in (folder/'ground.json',folder/'stages.json')},
                 cli_sha256=sha(a.cli),physics_or_threshold_modified_to_avoid_stop=False)
        if (folder/'cli-run.json').exists():
            result=read(folder/'cli-run.json');row.update(complete=result['complete'],failure_reason=result.get('failure_reason'))
        rows.append(row)
        # Trace even incomplete MC/MCC solves without asserting their schedule.
        trace=invoke([base/'probe',base/'mesh.json',folder/'ground.json',folder/'stages.json',folder/'raw','--diagnose'],folder/'trace.log',900)
        row['trace_process']=trace
    write(out/'scope.json',dict(cases=rows,new_physics=False,passed=False))

if __name__=='__main__':main()
