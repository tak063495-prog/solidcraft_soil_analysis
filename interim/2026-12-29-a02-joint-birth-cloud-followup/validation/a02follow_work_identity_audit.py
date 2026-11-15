"""Read NEW accepted fields only; retain the existing absolute work gate."""
import argparse,json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from a02follow_binding import fields,study
from a02follow_run import write,sha

def main(run,extra,out):
    result={};total=0
    for prefix,folder in [('formal',run),('contact-exchange',extra)]:
        for p in study.read(folder/'paths.json'):
            maximum=0.;count=0
            for _,f in fields.load_fields(folder/p['path']/(p['engine']+'-fields.jsonl')):
                balance=f['elastic_energy']-f['initial_energy']+f['plastic_dissipation']+f['numerical_dissipation']-f['endpoint_work']
                maximum=max(maximum,float(np.max(np.abs(balance),initial=0.)));count+=balance.size
            assert maximum<=1e-10,(prefix,p['path'],maximum)
            total+=count;result[prefix+'/'+p['path']]=dict(maximum_absolute_balance_kpa_m=maximum,gp_states=count)
    write(out,dict(passed=False,a02_complete=False,diagnostic_only=True,existing_absolute_work_gate_kpa_m=1e-10,
        limited_work_identity_check_passed=True,gp_states=total,paths=result,
        postsolve_audit_not_fresh_solve=True,source_sha256=sha(Path(__file__))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--extra',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();main(a.run,a.extra,a.output)
