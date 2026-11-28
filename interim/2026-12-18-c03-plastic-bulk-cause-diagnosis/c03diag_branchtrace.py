"""Read-only material substep replay at native accepted states; no equilibrium solve."""
import argparse
from collections import Counter
import json
from pathlib import Path
import zipfile
import numpy as np
from c03diag_study import read,write,sha,invoke,Model


def main():
    p=argparse.ArgumentParser();p.add_argument('--matrix',type=Path,required=True)
    p.add_argument('--solver',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--publish',type=Path,required=True);a=p.parse_args();a.output.mkdir();requests=[];actual=[];inputs={}
    for route in ('m0c0-native','m0c1-native','m0c2-native','m2c2-native'):
        folder=a.matrix/route;model=Model(read(a.matrix/'mesh.json'),read(a.matrix/'ground.json'),read(folder/'stages.json'))
        old=model.seeds
        paths=[folder/'raw/accepted.jsonl',folder/'raw/constitutive.jsonl'];inputs.update({route+'/'+p.name:sha(p) for p in paths})
        states=[json.loads(s) for s in paths[0].open()];diagnostics=[json.loads(s) for s in paths[1].open()]
        assert len(states)==len(diagnostics)==224
        for i,(state,d) in enumerate(zip(states,diagnostics)):
            for gp,(mat,before,e,n,after) in enumerate(zip(model.mats,old,d['strain_increment'],d['partitions'][:12],state['material_states'])):
                requests.append(dict(route=route,increment=i,gp=gp,parameters=mat,old=before,strain=e,count=n));actual.append(after)
            old=state['material_states']
    write(a.output/'request.json',requests)
    process=invoke([a.solver.resolve(),a.output.resolve()/'request.json',a.output.resolve()/'result.json'],a.output/'process.log',900)
    assert process['exit_code']==0
    result=read(a.output/'result.json');assert len(result)==len(actual)
    categories=Counter();byroute={};maxima={};mismatch=0
    for req,record,target in zip(requests,result,actual):
        assert (req['route'],req['increment'],req['gp'],req['count'])==(record['route'],record['increment'],record['gp'],record['count'])
        assert sum(record['returns'].values())==req['count'];categories.update(record['returns']);byroute.setdefault(req['route'],Counter()).update(record['returns'])
        for key,value in target.items():
            if key=='last_return':mismatch+=record['state'][key]!=value
            else:maxima[key]=max(maxima.get(key,0.),float(np.max(np.abs(np.asarray(value)-record['state'][key]))))
    summary=dict(diagnostic_only=True,global_equilibrium_solved=False,accepted_point_updates=len(actual),
          observed_local_substeps=sum(categories.values()),return_histogram=dict(categories),
          by_route={k:dict(v) for k,v in byroute.items()},maximum_final_state_difference=maxima,final_category_mismatches=mismatch,
          inputs_sha256=inputs,solver_sha256=sha(a.solver),process=process,passed=False)
    write(a.output/'summary.json',summary)
    archive=a.publish/'return-trace-evidence.zip'
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for name in ('request.json','result.json','summary.json','process.log'):z.write(a.output/name,name)
        z.write(a.solver,'return-trace-executable')
        z.write(Path(__file__),'c03diag_branchtrace.py')
        z.write(Path(__file__).resolve().parents[2]/'examples/c03_bulk_return_trace.rs','c03_bulk_return_trace.rs')
    summary['archive_sha256']=sha(archive);write(a.publish/'return-trace-summary.json',summary)

if __name__=='__main__':main()
