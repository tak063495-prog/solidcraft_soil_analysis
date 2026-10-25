"""Input-only column continuum/Q2 reference and independent linear 3D solve."""
import copy
from pathlib import Path
import json
import numpy as np
import fem_reference as fem
from verify_high_order import physical, constraints


def column_mesh(layers,quadratic=True,order=5):
    assert layers>=2 and layers%2==0
    root=Path(__file__).resolve().parents[1];mesh=json.loads((root/'examples/mc-column-mesh.json').read_text())
    points=[[float(x),float(y),k/layers] for k in range(layers+1) for y in (0,1) for x in (0,1)];tets=[];regions=[]
    for k in range(layers):
        corners=[4*k+i for i in (0,1,3,2)]+[4*(k+1)+i for i in (0,1,3,2)]
        for local in ((0,1,2,6),(0,2,3,6),(0,3,7,6),(0,7,4,6),(0,4,5,6),(0,5,1,6)):
            tets.append([corners[i] for i in local]);regions.append(1 if k<layers/2 else 2)
    mesh.update(nodes=points,tetrahedra=tets,region_ids=regions,material_ids=regions.copy())
    if quadratic:
        mesh['schema']='soil-fem-mesh/3'
        mids={};high=[]
        for tet in tets:
            row=tet.copy()
            for i,j in fem.EDGES:
                key=tuple(sorted((tet[i],tet[j])))
                if key not in mids:mids[key]=len(points);points.append(((np.array(points[key[0]])+points[key[1]])/2).tolist())
                row.append(mids[key])
            high.append(row)
        mesh.update(tetrahedra10=high,quadrature_order=order)
    else:mesh.pop('tetrahedra10',None);mesh.pop('quadrature_order',None)
    mesh['node_sets']={name:[i for i,p in enumerate(points) if name=='all' or abs(p[2]-(0 if name=='bottom' else 1))<1e-12] for name in ('all','bottom','top')}
    return mesh


def continuum(z,level,modulus=12000.,gamma=10.):
    z=np.asarray(z);wet=np.minimum(max(level,0.),z)
    return -gamma/modulus*(z-z*z/2-level*wet+wet*wet/2+max(level-1.,0.)*z)


def q2_column(layers,level,modulus=12000.,gamma=10.):
    count=2*layers+1;matrix=np.zeros((count,count));load=np.zeros(count);length=1/layers
    local=np.array([[7.,-8.,1.],[-8.,16.,-8.],[1.,-8.,7.]])*modulus/(3*length)
    x,w=np.polynomial.legendre.leggauss(3)
    for i in range(layers):
        ids=np.arange(2*i,2*i+3);matrix[np.ix_(ids,ids)]+=local
        left=i/layers;right=(i+1)/layers;cuts=sorted(set([left,right]+([level] if left<level<right else [])))
        for a,b in zip(cuts,cuts[1:]):
            z=(a+b)/2+(b-a)*x/2;xi=(z-left)/length
            derivative=np.array([-3+4*xi,4-8*xi,-1+4*xi])/length
            delta=gamma*(np.maximum(level-z,0)-(1-z))
            load[ids]+=derivative@(w*delta)*(b-a)/2
    load[-1]-=gamma*max(level-1.,0.)
    u=np.zeros(count);u[1:]=np.linalg.solve(matrix[1:,1:],load[1:]);return np.linspace(0,1,count),u


def linear_3d(mesh,ground,case,active,pressure,level,pore_pressure_kpa=None):
    assert all('mohr_coulomb' in m for m in ground['materials']) and case['settings']['kinematics']=='side_rollers'
    xyz=np.asarray(mesh['nodes']);materials={m['id']:m for m in ground['materials']};layers={r['body']:r for r in ground['layers']}
    frame_ground=copy.deepcopy(ground);frame_ground['initialization']['water_level_m']=level
    load,geometry=physical(mesh,frame_ground,np.asarray(active),pressure,case.get('pressure_on_original_surface_only',False))
    used,fixed=constraints(mesh,np.asarray(active),'side_rollers');free=np.flatnonzero(~fixed.reshape(-1));lookup=np.full(xyz.size,-1);lookup[free]=np.arange(len(free))
    k=np.zeros((len(free),len(free)));internal=np.zeros(xyz.size);operators=[]
    point_count=len(fem.volume(mesh,0)['weight'])
    if pore_pressure_kpa is not None:assert len(pore_pressure_kpa)==point_count*len(mesh['tetrahedra'])
    for ei in np.flatnonzero(active):
        q=fem.volume(mesh,ei);sigma,old_pore,material,layer=fem.initial_fields(mesh,ground,ei,q['position'])
        par=material['mohr_coulomb'];d=fem.elastic_tensor(par['young_kpa'],par['poisson']);total=sigma.copy()
        if pore_pressure_kpa is not None:total[:,:3]+=np.asarray(pore_pressure_kpa)[ei*point_count:(ei+1)*point_count,None]
        elif level is not None:total[:,:3]+=np.maximum(level-q['position'][:,2],0)[:,None]*ground['initialization']['water_unit_weight_kn_m3']
        force=-np.einsum('qai,qa,q->i',q['b'],total,q['weight']);ids=(q['nodes'][:,None]*3+np.arange(3)).reshape(-1)
        np.add.at(internal,ids,force)
        ke=np.einsum('qai,ab,qbj,q->ij',q['b'],d,q['b'],q['weight'],optimize=True)
        mapped=lookup[ids];live=np.flatnonzero(mapped>=0);k[np.ix_(mapped[live],mapped[live])]+=ke[np.ix_(live,live)]
        operators.append((ei,q,d,sigma,layer))
    rhs=(load.reshape(-1)-internal)[free];u=np.zeros(xyz.size);u[free]=np.linalg.solve(k,rhs)
    residual=np.linalg.norm(k@u[free]-rhs)/max(np.linalg.norm(rhs),1.);assert residual<1e-10
    stress={};void={}
    for ei,q,d,sigma,layer in operators:
        de=-np.einsum('qai,i->qa',q['b'],u.reshape(-1,3)[q['nodes']].reshape(-1));stress[ei]=sigma+de@d.T
        assert layer['void_ratio']['mode']=='explicit';v=1+layer['void_ratio']['value'];void[ei]=v*np.exp(-de[:,:3].sum(axis=1))-1
    return {'displacement_m':u.reshape(-1,3).tolist(),'stress':stress,'void_ratio':void,'true_linear_residual':float(residual),'geometry':geometry}
