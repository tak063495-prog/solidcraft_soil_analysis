"""Read-only displacement/traction accounting on SHA-verified saved histories."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import c04cause_metrics as metrics
from c04cause_driver import reference,write


def prefix(values):
    return np.asarray([[math.fsum(map(float,values[:i+1,g])) for g in range(values.shape[1])] for i in range(len(values))])


def account(native,independent,system):
    a=native[1];b=independent[1]
    delta=lambda v:np.diff(v,axis=0,prepend=np.zeros_like(v[:1]))
    tangent=lambda v:v-np.sum(v*system.normal,axis=2)[:,:,None]*system.normal
    ja,jb=tangent(a['jump']),tangent(b['jump']);qa,qb=tangent(a['traction']),tangent(b['traction'])
    ra=ja-qa/system.kt[:,None]-a['plastic_slip'];rb=jb-qb/system.kt[:,None]-b['plastic_slip']
    meanq=(qa+qb)*.5
    pieces=dict(traction_difference=prefix(np.sum((qa-qb)*(delta(a['plastic_slip'])+delta(b['plastic_slip']))*.5,axis=2)),
      displacement_jump_difference=prefix(np.sum(meanq*delta(ja-jb),axis=2)),
      elastic_traction_difference=prefix(np.sum(-meanq*delta(qa-qb)/system.kt[:,None],axis=2)),
      contact_reference_offset_difference=prefix(np.sum(-meanq*delta(ra-rb),axis=2)))
    actual=a['plastic_dissipation']-b['plastic_dissipation']
    endpoint=pieces['traction_difference']+pieces['displacement_jump_difference']+pieces['elastic_traction_difference']+pieces['contact_reference_offset_difference']
    pieces['within_increment_remainder']=actual-endpoint
    return dict(actual_absolute_l2=metrics.norm(actual),contributions={k:metrics.budget(actual,v) for k,v in pieces.items()},
      maximum_identity_defect=float(np.max(abs(actual-sum(pieces.values())))),
      reference_offset_definition='jump_t - traction_t/kt - stored_plastic_slip; an accounting variable, not a new physical law',
      vectors={**pieces,'actual':actual})


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--baseline',type=Path,required=True);p.add_argument('--report',type=Path,required=True)
    a=p.parse_args();read=lambda p:json.loads(p.read_text(encoding='utf-8'))
    report=dict(schema='c04-displacement-traction-accounting/1',passed=False,c04_complete=False,new_solve=False,
                new_acceptance_gate=False,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),levels={})
    for water in (8,16,32):
        n=a.baseline/f'native-w{water}';b=a.baseline/f'independent-w{water}'
        model=reference.Model(read(a.baseline/'mesh.json'),read(a.baseline/'ground.json'),read(n/'case.json'))
        report['levels'][str(water)]=account(metrics.load(n/'fields.jsonl'),metrics.load(b/'fields.jsonl'),model.system)
    write(a.report,report)

if __name__=='__main__':main()
