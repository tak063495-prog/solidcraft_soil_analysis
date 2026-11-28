"""Supported and free two-face shell patch, numerical drilling selection evidence."""
import argparse
import json
from pathlib import Path
import subprocess
import numpy as np
from scipy.linalg import block_diag,eigh,null_space,orth,subspace_angles
from scipy.spatial.transform import Rotation


def make(axes,t,alpha,distorted):
    param=np.array([[0,0],[1,0],[0,1],[.5,0],[.5,.5],[0,.5],[1,1],[1,.5],[.5,1]],float)
    if distorted: param[6]+=[-.12,.07]
    for i in [4,7,8]:
        u,v={4:(1,2),7:(1,6),8:(2,6)}[i]
        param[i]=(param[u]+param[v])/2
    x,y=param.T
    nodes=np.column_stack((x,y,.15*x*x-.10*y*y))@axes.T
    faces=[dict(id='a',nodes=[0,1,2,3,4,5],offset_m=0.,group='s'),
           dict(id='b',nodes=[2,1,6,4,7,8],offset_m=0.,group='s')]
    return dict(nodes=nodes.tolist(),faces=faces,directors=[],steps=[],
                section=dict(id='dry',young_kpa=20000.,poisson=.25,thickness_m=t,unit_weight_kn_m3=0.,drilling_factor=alpha))


def measure(binary,out):
    out.mkdir(parents=True,exist_ok=False)
    rows=[]
    for distorted in [False, True]:
        basis=np.eye(3)
        rotated=Rotation.from_rotvec([.21,-.41,.65]).as_matrix()
        for ir,axes in enumerate([basis,rotated]):
            for t in [.1,.01,.001]:
                for mode in ['tied','raw-membrane','raw-shear']:
                    for alpha in ([1e-5,1e-4,1e-3] if mode=='tied' else [1e-4]):
                        case=make(axes,t,alpha,distorted)
                        tag=f'{distorted}-{ir}-{t}-{alpha}-{mode}'
                        source=out/(tag+'.input.json')
                        source.write_text(json.dumps(case),encoding='utf-8')
                        dest=out/(tag+'.matrix.json')
                        subprocess.run([binary,source,mode,dest],check=True)
                        x=np.array(case['nodes'])
                        k=np.zeros((54,54));drill=np.zeros_like(k);d=np.zeros_like(x)
                        for elem,face in zip(json.loads(dest.read_bytes()),case['faces']):
                            ids=np.array(face['nodes']);idx=(6*ids[:,None]+np.arange(6)).ravel()
                            k[np.ix_(idx,idx)]+=np.reshape(elem['stiffness_columns'],(36,36),order='F')
                            drill[np.ix_(idx,idx)]+=np.reshape(elem['drill_columns'],(36,36),order='F')
                            d[ids]=elem['directors']
                        rigid=np.empty((54,6))
                        for i,omega in enumerate(np.eye(3)):
                            rigid[:,i]=np.tile(np.r_[omega,[0.]*3],9)
                            rigid[:,3+i]=np.column_stack((np.cross(omega,x),np.tile(omega,(9,1)))).ravel()
                        values,vectors=eigh(k)
                        physical=k-drill
                        free=np.array([j for j in range(54) if j//6 not in [0,1,3]])
                        # Five physical DOFs/node. Rotations live in each node's tangent plane.
                        lift=block_diag(*[block_diag(np.eye(3),null_space(d[i][None,:])) for i in range(9)])
                        active=[j for j in range(45) if j//5 not in [0,1,3]]
                        L=lift[:,active]
                        supported=k[np.ix_(free,free)]
                        # Apply the same load after rotating the whole patch.
                        loads=[]
                        for direction in [axes[:,0],axes[:,2]]:
                            f=np.zeros(54);f[36:39]=direction
                            q=np.zeros(54);q[free]=np.linalg.solve(supported,f[free])
                            ref=L@np.linalg.solve(L.T@physical@L,L.T@f)
                            loads.append(dict(displacement6=float(f@q),displacement5=float(f@ref),
                                physical_energy=float(q@physical@q/2),five_dof_energy=float(ref@physical@ref/2),
                                drilling_work=float(q@drill@q/2),free_residual=float(np.linalg.norm((k@q-f)[free]))))
                        rows.append(dict(distorted=distorted,rotation=ir,t_over_L=t,alpha=alpha,tying=mode,
                            eigenvalues=values.tolist(),first_positive_mode=vectors[:,6].tolist(),
                            principal_angles_rad=subspace_angles(orth(rigid),vectors[:,:6]).tolist(),
                            Kr_norm=np.linalg.norm(k@rigid,axis=0).tolist(),
                            supported_condition=float(np.linalg.cond(supported)),loads=loads))
    (out/'modes.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('observer',type=Path)
    p.add_argument('output',type=Path)
    args=p.parse_args()
    measure(args.observer.resolve(),args.output)
