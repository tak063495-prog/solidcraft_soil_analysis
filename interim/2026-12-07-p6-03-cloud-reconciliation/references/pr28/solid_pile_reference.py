"""Small standalone 3D Hex8 elasticity reference for the historical square pile.

Own tensor-product quadrature, assembly, constraints, loads and sparse solve.
No Rust/old Python solver imports. Large refinements require --pc-series.
Sources and small-strain increment interpretation are documented in README.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

SIGNS = np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],
                  [-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]], dtype=float)


def tensor(young, poisson):
    mu=young/(2*(1+poisson))
    lam=young*poisson/((1+poisson)*(1-2*poisson))
    d=np.zeros((6,6))
    d[:3,:3]=lam
    d[np.arange(3),np.arange(3)]+=2*mu
    d[3:,3:]=np.eye(3)*mu
    return d


def shape(point, widths):
    factors=1+SIGNS*np.asarray(point)
    n=np.prod(factors,axis=1)/8
    grad=np.column_stack([SIGNS[:,i]*np.prod(factors[:,[j for j in range(3) if j!=i]],axis=1)/(4*widths[i])
                          for i in range(3)])
    b=np.zeros((6,24))
    b[0,0::3]=grad[:,0];b[1,1::3]=grad[:,1];b[2,2::3]=grad[:,2]
    b[3,0::3]=grad[:,1];b[3,1::3]=grad[:,0]
    b[4,1::3]=grad[:,2];b[4,2::3]=grad[:,1]
    b[5,0::3]=grad[:,2];b[5,2::3]=grad[:,0]
    return n,b


def element(widths, d):
    g,w=leggauss(2)
    k=np.zeros((24,24))
    for i in range(2):
        for j in range(2):
            for l in range(2):
                _,b=shape([g[i],g[j],g[l]],widths)
                k+=b.T@d@b*np.prod(widths)/8*w[i]*w[j]*w[l]
    return k


def axis(points, n):
    return np.concatenate([np.linspace(a,b,n,endpoint=False) for a,b in zip(points[:-1],points[1:])]+[np.array([points[-1]])])


def mesh(n):
    axes=[axis([0.,4.5,5.5,10.],n),axis([0.,4.5,5.5,10.],n),axis([0.,2.,8.],n)]
    x,y,z=axes
    xyz=np.array(np.meshgrid(x,y,z,indexing='ij')).reshape(3,-1).T
    index=np.arange(len(xyz)).reshape(len(x),len(y),len(z))
    cells=[];pile=[];widths=[]
    for i in range(len(x)-1):
        for j in range(len(y)-1):
            for k in range(len(z)-1):
                cells.append([index[i,j,k],index[i+1,j,k],index[i+1,j+1,k],index[i,j+1,k],
                              index[i,j,k+1],index[i+1,j,k+1],index[i+1,j+1,k+1],index[i,j+1,k+1]])
                pile.append(n <= i < 2*n and n <= j < 2*n and k >= n)
                widths.append([x[i+1]-x[i],y[j+1]-y[j],z[k+1]-z[k]])
    return xyz,np.array(cells),np.array(pile),np.array(widths)


def solve(n):
    if isinstance(n,bool) or not isinstance(n,int) or not 1 <= n <= 16:
        raise ValueError('integer subdivision 1..16 required')
    started=time.perf_counter()
    xyz,cells,pile,widths=mesh(n)
    size=xyz.size
    rows=[];cols=[];values=[]
    gravity=np.zeros(size);increment=np.zeros(size);head_weights=np.zeros(len(xyz))
    volume=np.zeros(2)
    ds=[tensor(30_000.,.3),tensor(30_000_000.,.2)]
    for nodes,is_pile,dimensions in zip(cells,pile,widths):
        material=int(is_pile)
        ids=(nodes[:,None]*3+np.arange(3)).ravel()
        ke=element(dimensions,ds[material])
        rows.append(np.repeat(ids,24));cols.append(np.tile(ids,24));values.append(ke.ravel())
        v=float(np.prod(dimensions));volume[material]+=v
        gravity[3*nodes+2]-=[20.,24.][material]*v/8
        if is_pile and np.all(xyz[nodes[4:],2]==8.):
            area=dimensions[0]*dimensions[1]
            head_weights[nodes[4:]]+=area/4
            increment[3*nodes[4:]]+=10.*area/4
            increment[3*nodes[4:]+2]-=100.*area/4
    matrix=coo_matrix((np.concatenate(values),(np.concatenate(rows),np.concatenate(cols))),shape=(size,size)).tocsc()
    constrained=np.zeros_like(xyz,dtype=bool)
    constrained[:,0]=(xyz[:,0]==0.)|(xyz[:,0]==10.)
    constrained[:,1]=(xyz[:,1]==0.)|(xyz[:,1]==10.)
    constrained[:,2]=xyz[:,2]==0.
    free=np.flatnonzero(~constrained.ravel())
    lu=splu(matrix[free][:,free])
    forces=np.column_stack((gravity,increment))
    u=np.zeros_like(forces)
    u[free]=lu.solve(forces[free])
    residual=matrix@u-forces
    rel=np.linalg.norm(residual[free],axis=0)/np.linalg.norm(forces[free],axis=0)
    reaction=residual.reshape(-1,3,2).sum(axis=0)
    head_increment=head_weights@u[:,1].reshape(-1,3)/head_weights.sum()
    head_gravity=head_weights@u[:,0].reshape(-1,3)/head_weights.sum()
    return dict(subdivisions=n, elements=len(cells),nodes=len(xyz),free_dofs=len(free),
                h_max_m=float(widths.max()),pile_cross_section_h_m=1/n,
                volumes_m3=volume.tolist(), head_area_m2=float(head_weights.sum()),
                applied_gravity_kn=gravity.reshape(-1,3).sum(axis=0).tolist(),
                applied_increment_kn=increment.reshape(-1,3).sum(axis=0).tolist(),
                head_gravity_displacement_m=head_gravity.tolist(),
                head_increment_displacement_m=head_increment.tolist(),
                head_additional_settlement_mm=float(-1000*head_increment[2]),
                free_relative_residual=rel.tolist(), support_reaction_kn=reaction.tolist(),
                support_reaction_columns=['gravity','load_increment'],
                support_reaction_rows=['global_x','global_y','global_z'],
                elapsed_seconds=time.perf_counter()-started,accepted=None,tolerance='未定義')


def report(ns=(1,2,4), progress=None):
    if len(ns)<3 or any(b!=2*a for a,b in zip(ns[:-1],ns[1:])):
        raise ValueError('at least three successive doubling subdivisions required')
    results=[]
    for n in ns:
        result=solve(n)
        results.append(result)
        if progress is not None:
            progress(result)
    settlements=[r['head_additional_settlement_mm'] for r in results]
    orders=[]
    for a,b,c in zip(settlements[:-2],settlements[1:-1],settlements[2:]):
        ratio=(b-a)/(c-b) if c!=b else 0.
        orders.append(float(np.log2(ratio)) if ratio>0 else None)
    return dict(schema='pile-independent-small-3d/1', complete=True,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                physical_conditions=dict(units='m-kN-kPa',soil_box_m=[[0,10],[0,10],[0,8]],
                    pile_box_m=[[4.5,5.5],[4.5,5.5],[2,8]],soil_E_kpa=30000.,soil_poisson=.3,
                    pile_E_kpa=30000000.,pile_poisson=.2,soil_gamma_kn_m3=20.,pile_gamma_kn_m3=24.,
                    interface='perfect bond shared nodes',water='dry',initial_stress='zero before gravity',
                    boundary='x sides ux=0; y sides uy=0; bottom uz=0',head_traction_kpa=[10.,0.,-100.]),
                model='3d-solid-elastic-bonded',
                scope='same finite box, square pile, elasticity, rollers and traction; own Hex8, no project solvers',
                reference_state='gravity equilibrium followed by linear load increment; unload returns to gravity by linearity',
                series=results,observed_head_settlement_order=orders,
                order_definition='log2((s_2h-s_h)/(s_h-s_h/2)); asymptotic regime NOT established',
                accepted=None,tolerance='未定義',p6_complete=False,
                original_record_position='判定不能: approximate Hex8 series is not an exact/error-bounded 3D reference')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--subdivisions',type=int,nargs='+',default=[1,2,4])
    parser.add_argument('--pc-series',action='store_true')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if max(args.subdivisions)>4 and not args.pc_series:
        parser.error('subdivisions >4 reserved for user PC; --pc-series required')
    result=report(tuple(args.subdivisions), progress=lambda r: print(json.dumps(r,allow_nan=False),flush=True))
    with args.output.open('x',encoding='utf-8') as stream:
        json.dump(result,stream,ensure_ascii=False,indent=2,allow_nan=False);stream.write('\n')


if __name__=='__main__':main()
