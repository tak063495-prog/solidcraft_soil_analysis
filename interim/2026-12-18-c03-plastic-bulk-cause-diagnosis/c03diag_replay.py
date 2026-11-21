"""Frozen-path local controls. These fields do not solve global equilibrium."""
import os
for _key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[_key]='1'
import argparse
from pathlib import Path
import copy
from collections import Counter
from types import SimpleNamespace
import numpy as np
from c03diag_study import read,write,Model
from c03diag_reference import integrate
from c03diag_metrics import load,norm,budget
from verify_joint_coupling import contact


def contact_replay(model, fields, count):
    system=model.system;sub=SimpleNamespace(**{k:getattr(system,k) for k in ('normal','kn','kt','mu')})
    old=np.zeros_like(fields['jump'][0]);q=old.copy();slip=old.copy();values=[];tractions=[]
    for target in fields['jump']:
        result,dp,_=contact(sub,old,q,target,count);q=result[:,0];slip=slip+dp[:,0]
        old=target;values.append(slip.copy());tractions.append(q.copy())
    return dict(plastic_slip=np.asarray(values),traction=np.asarray(tractions))


def material_replay(model, fields, count):
    states=copy.deepcopy(model.seeds);previous=np.zeros_like(fields['strain'][0]);stress=[];plastic=[];returns=Counter()
    for strain in fields['strain']:
        de=-(strain-previous)
        states=[integrate(mat,s,e,count) for mat,s,e in zip(model.mats,states,de)]
        returns.update(s['last_return'] for s in states)
        stress.append([s['effective_stress_compression_positive_kpa'] for s in states])
        plastic.append([s['plastic_strain_compression_positive'] for s in states]);previous=strain
    return dict(stress=np.asarray(stress),plastic_strain=np.asarray(plastic)),dict(returns)


def main():
    p=argparse.ArgumentParser();p.add_argument('--matrix',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir();base=a.matrix
    model=Model(read(base/'mesh.json'),read(base/'ground.json'),read(base/'m0c0-native/stages.json'))
    data={};replayed={};materials={};report=dict(local_only=True,substitute_for_global_equilibrium=False,contact={},material={},returns={})
    for kind in ('native','independent'):
        for c in range(4):
            name=f'm0c{c}-{kind}';_,fields=load(base/name/'fields.jsonl');data[name]=fields
            for count in (16,32):
                values=contact_replay(model,fields,count);replayed[name,count]=values
                np.savez(a.output/f'{name}-contact{count}.npz',**values)
        for c in range(3):
            x=f'm0c{c}-{kind}';y=f'm0c{c+1}-{kind}'
            original=data[x]['plastic_slip']-data[y]['plastic_slip']
            for count in (16,32):
                delta=replayed[x,count]['plastic_slip']-replayed[y,count]['plastic_slip']
                report['contact'][f'{x}--{y}-fixed{count}']=budget(original,delta)
            report['contact'][f'{x}--{y}-replay16--32']=budget(original,
                 (replayed[x,16]['plastic_slip']-replayed[y,16]['plastic_slip'])-
                 (replayed[x,32]['plastic_slip']-replayed[y,32]['plastic_slip']))
    for c in range(3):
        name=f'm0c{c}-native'
        for count in (8,32):
            print('material replay',name,count,flush=True)
            values,branches=material_replay(model,data[name],count);materials[name,count]=values
            np.savez(a.output/f'{name}-material{count}.npz',**values);report['returns'][f'{name}-fixed{count}']=branches
    for c in range(2):
        x=f'm0c{c}-native';y=f'm0c{c+1}-native'
        for field in ('stress','plastic_strain'):
            original=data[x][field]-data[y][field]
            for count in (8,32):
                report['material'][f'{x}--{y}-{field}-fixed{count}']=budget(original,materials[x,count][field]-materials[y,count][field])
    write(a.output/'summary.json',report)

if __name__=='__main__':main()
