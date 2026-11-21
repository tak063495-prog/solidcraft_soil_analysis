"""Fixed-order binary64 reassembly controls; no global equilibrium solve."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
from a01hp_common import write
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'interim/2026-12-28-a01-joint-combo-cloud-followup'))
from a01follow_contact import point


def norm(a):
    import math
    return math.sqrt(math.fsum(float(x)**2 for x in np.asarray(a).ravel(order='C')))


def run(source, output):
    folder=source/'formal/stacked__native__Pnative__Rnative__g8'
    read=lambda path:json.loads(path.read_text(encoding='utf-8'))
    mesh=read(folder/'mesh.json');ground=read(folder/'ground.json')
    kernel=read(source/'separation/native-kernel-output.json')
    events=[json.loads(line) for line in (folder/'native/progress.jsonl').read_text().splitlines()]
    tangent=[v for v in events if v['kind']=='tangent'][-1]
    t=np.zeros((len(tangent['dof_rows']),48))
    for i,row in enumerate(tangent['dof_rows']):
        for j,c in row:t[i,j]=c
    full=np.zeros((t.shape[0],t.shape[0]));materials={m['id']:m['mohr_coulomb'] for m in ground['materials']}
    for e,cell in enumerate(mesh['tetrahedra']):
        x=np.c_[np.ones(4),np.asarray(mesh['nodes'])[cell]];inv=np.linalg.inv(x);volume=abs(np.linalg.det(x))/6
        b=np.zeros((6,12))
        for a in range(4):
            dx,dy,dz=inv[1:,a];b[0,3*a]=dx;b[1,3*a+1]=dy;b[2,3*a+2]=dz
            b[3,3*a:3*a+2]=dy,dx;b[4,3*a+1:3*a+3]=dz,dy;b[5,3*a]=dz;b[5,3*a+2]=dx
        m=materials[mesh['material_ids'][e]];young,nu=m['young_kpa'],m['poisson'];lam=young*nu/((1+nu)*(1-2*nu));shear=young/(2*(1+nu))
        d=np.zeros((6,6));d[:3,:3]=lam
        for a in range(3):d[a,a]+=2*shear;d[a+3,a+3]=shear
        k=b.T@d@b*volume;ids=[3*n+a for n in cell for a in range(3)]
        for i,ii in enumerate(ids):
            for j,jj in enumerate(ids):full[ii,jj]+=k[i,j]
    bulk=t.T@full@t;mats={name:bulk.copy() for name in ('analytic','fd_h_div_1','fd_h_div_2','fd_h_div_4','saved_native_local_fd')}
    for q,loc in enumerate(kernel['local_inputs']):
        c=np.zeros((3,t.shape[0]))
        for n,coefficient in kernel['jump_coefficients'][q]:
            for a in range(3):c[a,3*n+a]=coefficient
        cr=c@t;m=loc['material'];old=loc['old'];args=(old['normal'],m['normal_stiffness_kpa_m'],m['tangential_stiffness_kpa_m'],m['friction_coefficient'],old['jump_m'],old['traction_kpa'],old['birth']['elastic_offset_m'])
        target=np.array(loc['target']);count=loc['count'];derivatives={'analytic':point(*args,target,count)[1], 'saved_native_local_fd':np.array(kernel['records'][q]['tangent_columns']).T}
        for divisor in (1,2,4):
            h=loc['h']/divisor;d=np.zeros((3,3))
            for a in range(3):
                step=np.eye(3)[a]*h
                d[:,a]=(point(*args,target+step,count)[0]-point(*args,target-step,count)[0])/(2*h)
            derivatives['fd_h_div_'+str(divisor)]=d
        for name,d in derivatives.items():mats[name]+=cr.T@d@cr*kernel['layout'][q]['weight_m2']
    saved=np.array(tangent['matrix_columns']).T;mats['saved_pre_predictor']=saved
    write(output,dict(numpy_version=np.__version__,passed=False,a01_complete=False,diagnostic_only=True,
        matrix_product_threads=1,aggregation_order='element/GP then row/column, norm C-order math.fsum',
        records={name:dict(matrix=a.tolist(),spectrum_descending=np.linalg.svd(a,compute_uv=False).tolist(),
            spectrum_with_vectors_descending=np.linalg.svd(a,compute_uv=True)[1].tolist(),
            saved_pre_predictor_difference=norm(a-saved)) for name,a in mats.items()}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.source,a.output)
