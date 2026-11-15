"""Fully exchanged active C-04 predictor/repartition policies, global solves."""
import argparse
import sys
from a04follow_common import *
sys.path.insert(0,str(ROOT/'tools'))
OLD=ROOT/'interim/2026-11-22-wet-joint-transient-water'
sys.path.insert(0,str(OLD))
import a04follow_c04_reference as reference

def legacy():
    return load('a04follow_legacy_c04_study',OLD/'study.py')

def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--engine',choices=('native','independent'),required=True)
    p.add_argument('--p',choices=('N','I'),required=True);p.add_argument('--r',choices=('N','I'),required=True)
    a=p.parse_args();s=legacy();folder=a.run/(a.engine+'-p'+a.p+'r'+a.r)
    if a.engine=='native':
        s.native(a.run,folder,32,32,3600)
        run=read(folder/'raw/run.json')
        attempts=run['geostatic']['attempts']+run['attempts']
        predictors=[q for v in attempts for q in v.get('contact_predictors',[])]
        if predictors:raise RuntimeError(('auxiliary predictor exercised: exchange incomplete',predictors))
        write(folder/'predictor-audit.json',dict(auxiliary_predictors=0,accepted_attempts=len(attempts),supported_case=True))
    else:
        folder.mkdir();case=s.mechanical_case(32,32);write(folder/'case.json',case)
        model=reference.Model(read(a.run/'mesh.json'),read(a.run/'ground.json'),case)
        f,emit=sink(folder/'fields.jsonl')
        try:result=reference.run(model,32,case['settings']['material_maximum_substep'],case['joint_system']['maximum_increment_m'],emit,predictor=a.p,repartition=a.r)
        finally:f.close()
        write(folder/'result.json',result)

if __name__=='__main__':main()
