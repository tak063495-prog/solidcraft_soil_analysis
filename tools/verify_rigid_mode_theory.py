"""Analytic free-mode counts checked by independent NumPy contact SVD.

Synthetic contact fields certify only kinematics, not constitutive history or
the ability to initialize a free block from a zero-pressure contact.
"""
import hashlib
import json
from pathlib import Path
import numpy as np
import fem_reference as ref
from joint_reference import System


def verify():
    root=Path(__file__).resolve().parents[1];records=[];inputs={}
    configuration=json.loads((root/'examples/stages-joint-column.json').read_text())['joint_system']
    for filename in ('joint-column-mesh.json','joint-column-tet10-mesh.json'):
        path=root/'examples'/filename;original=json.loads(path.read_text());inputs[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
        for rotated in (False,True):
            mesh=json.loads(json.dumps(original));xyz=np.asarray(mesh['nodes'])
            angle=.47;c,s=np.cos(angle),np.sin(angle)
            rotation=np.array([[c,-s,0],[s,c,0],[0,0,1]])@np.array([[c,0,s],[0,1,0],[-s,0,c]]) if rotated else np.eye(3)
            xyz=xyz@rotation.T;mesh['nodes']=xyz.tolist()
            lower=set(int(n) for e,r in enumerate(mesh['region_ids']) if r==1 for n in ref.cells(mesh)[e])
            mesh['node_sets']['fixed_lower']=sorted(lower)
            case={'settings':{'kinematics':'explicit','relative_force_tolerance':1e-8},'boundary':{'constraints':[
                {'node_set':'fixed_lower','direction':axis,'value_m':0.} for axis in np.eye(3).tolist()]}}
            for mu,jump,missing in ((.3,[0,0,0],6),(.3,[0,0,.001],6),(0,[0,0,-.0001],3),
                                    (.3,[0,0,-.0001],0),(.3,[.00015,0,-.0001],1),(.3,[.001,0,-.0001],1),(.3,[0,0,-1e-15],6)):
                config=json.loads(json.dumps(configuration));config['joints'][0]['material']['friction_coefficient']=mu
                system=System(mesh,config);d=rotation@np.array(jump);points=[]
                for p in system.points:
                    n=np.asarray(p['normal']);gn=d@n;pressure=-1e6*min(gn,0);tau=2e5*(d-gn*n)
                    length=np.linalg.norm(tau)
                    if length>mu*pressure:tau*=mu*pressure/length
                    points.append({'traction_kpa':(tau-pressure*n).tolist()})
                u=np.asarray([[0,0,0] if n in lower else d for n in range(len(xyz))]).tolist()
                count=len(mesh['tetrahedra'])*len(ref.volume_rule(mesh)[1])
                state={'active_elements':[True]*len(mesh['tetrahedra']),'displacement_m':u,'joints':{'points':points},
                       'material_states':[{'effective_stress_compression_positive_kpa':[100,100,100,0,0,0]} for _ in range(count)]}
                report=system.rigid_support(case,state)
                assert report['unrestrained_modes']==missing,(filename,rotated,mu,jump,report)
                records.append({'mesh':filename,'rotated':rotated,'friction':mu,'unrotated_jump_m':jump,
                                'expected_free_modes':missing,**report})
    for path in (Path(__file__),root/'tools/joint_reference.py',root/'tools/fem_reference.py',root/'examples/stages-joint-column.json'):
        inputs[str(path.resolve())]=hashlib.sha256(path.read_bytes()).hexdigest()
    return {'schema':'soil-fem-rigid-mode-theory/1','passed':True,'cases':records,'sha256':inputs}


if __name__=='__main__':print(json.dumps(verify(),indent=2))
