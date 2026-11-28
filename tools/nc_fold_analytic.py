"""A known nonsymmetric 2-DOF fold using the existing NC continuation driver."""
import argparse
import json
from pathlib import Path
import numpy as np
from fold_stability import point,judge
from nc_fold_path import load_driver,control_binding
from analyze_nc_fold import linear_audit


class Analytic:
    kind='analytic'
    R=np.eye(2)
    L=np.array([-1.,0.])
    contacts=[]
    def __init__(self):
        self.u=np.array([-.025,.0175]);self.old=self.u.copy();self.lam=.95-self.u[0]**2
        k=np.array([[2*self.u[0],0.],[.7,1.]])
        self.t=np.r_[np.linalg.solve(k,self.L),1.];self.t/=np.linalg.norm(self.t);self.m=self
    def load(self,lam):return np.array([.95-lam,0.])
    def evaluate(self,du,old,label):
        x=np.array(old)+du
        return np.array([x[0]**2,x[1]+.7*x[0]]),np.array([[2*x[0],0.],[.7,1.]]),x
    def serial(self,state):return np.asarray(state).tolist()
    def state(self,state):return np.asarray(state)
    def audit(self,*args):return dict(diagnostic_only=True)


def run(output):
    output.mkdir(parents=True,exist_ok=False);driver=load_driver()
    binding=control_binding(dict(equation='R=[u^2+lambda-0.95,y+0.7u]',known_fold=[0.,0.,.95],step_size_m=.005))
    for name in ('whole','first','resumed'):(output/name).mkdir()
    driver.follow(Analytic(),output/'whole',.005,12,30.,binding)
    driver.follow(Analytic(),output/'first',.005,6,30.,binding)
    driver.follow(Analytic(),output/'resumed',.005,6,30.,binding,output/'first/restart.json')
    readrows=lambda path:[json.loads(line) for line in path.read_text(encoding='utf-8').splitlines()]
    samples=readrows(output/'whole/samples.jsonl');systems=readrows(output/'whole/linear-systems.jsonl');points=[]
    for row in samples:
        actual=[s for s in systems if s['step']==row['step']];tangent=actual[-1];a=np.asarray(tangent['matrix'])
        d=point(a[:-1,:-1],a[:-1,-1],a[-1]/1000.)
        audits=linear_audit(actual)
        points.append(dict(accepted=True,**{'lambda':row['progress']},diagnostic=d,relative_force_residual=row['relative_force_residual'],
                           force_block_linear_residual=max(s['force_block'] for s in audits),augmented_linear_residual=max(s['augmented'] for s in audits)))
    result=dict(known_fold=dict(u=0.,y=0.,lambda_value=.95),equations=binding['equation'],
                verdict=judge(points),points=points,
                actual_checkpoint_restart_all_saved_fields_exact=samples==readrows(output/'first/samples.jsonl')+readrows(output/'resumed/samples.jsonl'),
                regular_prefix_verdict=judge(points[:3]),physical_nc_acceptance=False)
    (output/'analytic-fold.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    r=run(p.parse_args().output);print(json.dumps(dict(verdict=r['verdict'],restart=r['actual_checkpoint_restart_all_saved_fields_exact']),indent=2))
