"""Bounded Tet10 nodal LP proposal; only full-GP verifiers certify a witness.

Depth-tied corner and midside values remain independent quadratic DOFs.
The polygon is inside the MC dual cone. Sampling the four reference vertices
is sufficient for the affine strain of straight Tet10, but can overrestrict
the discrete quadrature problem. The optional Gauss mode uses the actual
positive Jacobi point family. Neither a zero objective nor LP failure means
physical stability, collapse, or nonlinear FEM convergence.
"""
from pathlib import Path
import argparse,copy,datetime,hashlib,json,math,os,time
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np
from scipy import sparse
from scipy.optimize import linprog
from scipy.special import roots_jacobi

EDGES=[(0,1),(1,2),(2,0),(0,3),(1,3),(2,3)]
DL=np.array([[-1.,-1.,-1.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
def strict(q):
    def constant(v):raise ValueError('nonfinite JSON '+v)
    def pairs(items):
        result={}
        for k,v in items:
            if k in result:raise ValueError('duplicate key '+k)
            result[k]=v
        return result
    return json.loads(Path(q).read_text('utf-8-sig'),parse_constant=constant,object_pairs_hook=pairs)
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
def derivatives(natural):
    r,s,t=natural;L=np.array([1-r-s-t,r,s,t])
    return np.vstack([(4*L[:,None]-1)*DL,
        np.array([4*(L[i]*DL[j]+L[j]*DL[i]) for i,j in EDGES])])
def sampling(mode):
    if mode=='vertices':return np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    nodes=[(roots_jacobi(3,alpha,0)[0]+1)/2 for alpha in [2,1,0]]
    return np.array([[a,(1-a)*b,(1-a)*(1-b)*c] for a in nodes[0] for b in nodes[1] for c in nodes[2]])
def main():
    parser=argparse.ArgumentParser();parser.add_argument('mesh',type=Path);parser.add_argument('problem',type=Path);parser.add_argument('output',type=Path)
    parser.add_argument('--sides',type=int,default=64);parser.add_argument('--sampling',choices=['vertices','gauss'],default='vertices')
    parser.add_argument('--time-limit',type=float,default=45.);args=parser.parse_args()
    assert args.sides in [16,32,64,128] and 1.<=args.time_limit<=45.
    meta=args.output.with_suffix('.generation.json');assert not args.output.exists() and not meta.exists()
    args.output.parent.mkdir(parents=True,exist_ok=True);started=time.monotonic()
    mesh=strict(args.mesh);problem=strict(args.problem)
    assert problem['schema']=='soil-fem-mc-virtual-work-bound/1' and problem['water_policy']=='dry'
    xyz=np.asarray(mesh['nodes'],dtype=float);tets=np.asarray(mesh['tetrahedra'],dtype=int);high=np.asarray(mesh['tetrahedra10'],dtype=int)
    assert xyz.ndim==2 and xyz.shape[1]==3 and np.isfinite(xyz).all()
    assert high.shape==(len(tets),10) and np.array_equal(high[:,:4],tets)
    assert args.sampling!='gauss' or mesh.get('quadrature_order')==3
    active=problem['active_elements'];assert len(active)==len(tets) and all(type(v)is bool for v in active)
    selected=np.flatnonzero(active);assert 0<len(selected)<=20000
    natural=sampling(args.sampling);dref=np.asarray([derivatives(q) for q in natural]);dcenter=derivatives([.25,.25,.25])
    row_count=len(selected)*len(natural)*args.sides;assert row_count<=3_000_000,'bounded proposal matrix too large'
    materials={}
    for m in problem['materials']:
        assert type(m['id'])is int and m['id']>0 and m['id'] not in materials
        params=m['parameters'];c=float(params['cohesion_kpa']);phi=math.radians(params['friction_deg']);psi=math.radians(params['dilation_deg'])
        assert all(math.isfinite(v) for v in [c,phi,psi]) and c>=0 and 1e-8<phi<math.radians(89.) and 0<=psi<=phi
        if params['flow']=='davis':
            beta=math.cos(phi)*math.cos(psi)/(1-math.sin(phi)*math.sin(psi));c*=beta;phi=math.atan(beta*math.tan(phi))
        else:assert params['flow'] in ['associated','non_associated']
        materials[m['id']]=(c,phi)
    force=np.asarray(problem['physical_nodal_force_kn'],dtype=float);assert force.shape==(3*len(xyz),) and np.isfinite(force).all()
    used=sorted(set(high[selected].reshape(-1).tolist()));used_set=set(used)
    assert all(np.all(force[3*n:3*n+3]==0) for n in range(len(xyz)) if n not in used_set)
    groups={};group_of={}
    for n in used:
        key=tuple((0. if xyz[n,i]==0 else float(xyz[n,i])).hex() for i in [0,2])
        if key not in groups:groups[key]=len(groups)
        group_of[n]=groups[key]
    fixed=set()
    for c in problem['constraints']:
        assert c['value_m']==0. and c.get('gradient',[0.,0.,0.])==[0.,0.,0.]
        direction=c['direction'];assert len(direction)==3 and all(math.isfinite(v) for v in direction)
        axes=[i for i,v in enumerate(direction) if v!=0];assert len(axes)==1
        if axes[0]!=1:fixed.update((group_of[n],axes[0]) for n in mesh['node_sets'][c['node_set']] if n in used_set)
    variables={item:i for i,item in enumerate((g,axis) for g in range(len(groups)) for axis in [0,2] if (g,axis) not in fixed)}
    count=len(variables);assert count>0
    external=np.zeros(count);capacity=np.zeros(count)
    for n in used:
        for axis in [0,2]:
            index=variables.get((group_of[n],axis))
            if index is not None:external[index]+=force[3*n+axis]
    angles=np.arange(args.sides)*2*math.pi/args.sides;cosine=np.cos(angles);sine=np.sin(angles)
    inflation=(1+1e-4)/math.cos(math.pi/args.sides);row_parts=[];col_parts=[];data_parts=[];cells=[];offset=0
    geometry_scale=max(1.,float(np.max(np.abs(xyz))));maximum_depth_derivative=0.
    for ei in selected:
        corners=tets[ei];nodes=high[ei];x=xyz[corners];J=(x[1:]-x[0]).T;det=float(np.linalg.det(J));assert det>0
        for k,(i,j) in enumerate(EDGES):assert np.max(np.abs(xyz[nodes[k+4]]-(x[i]+x[j])/2))<=1e-12*geometry_scale
        inverse=np.linalg.inv(J);grads=dref@inverse;center=dcenter@inverse;volume=det/6
        c,phi=materials[mesh['material_ids'][ei]];apex=c/math.tan(phi)
        indices=np.array(sorted({variables[(group_of[int(n)],axis)] for n in nodes for axis in [0,2] if (group_of[int(n)],axis) in variables}),dtype=np.int32)
        lookup={int(index):i for i,index in enumerate(indices)}
        trace=np.zeros((len(natural),len(indices)));difference=trace.copy();shear=trace.copy();mean_trace=np.zeros(len(indices));dy={}
        for local,n in enumerate(nodes):
            group=group_of[int(n)];dy.setdefault(group,np.zeros(len(natural)));dy[group]+=grads[:,local,1]
            for axis in [0,2]:
                global_index=variables.get((group,axis))
                if global_index is None:continue
                j=lookup[global_index];trace[:,j]+=grads[:,local,axis];mean_trace[j]+=center[local,axis]
                difference[:,j]+=(1 if axis==0 else -1)*grads[:,local,axis]
                shear[:,j]+=grads[:,local,2 if axis==0 else 0]
        derivative=max((float(np.max(np.abs(v))) for v in dy.values()),default=0.)
        maximum_depth_derivative=max(maximum_depth_derivative,derivative)
        assert derivative<=1e-12*max(1.,float(np.max(np.abs(grads)))),'quadratic nodal tying is not plane strain'
        capacity[indices]+=volume*apex*mean_trace
        local=(math.sin(phi)*inflation*(difference[:,None,:]*cosine[None,:,None]+shear[:,None,:]*sine[None,:,None])-trace[:,None,:]).reshape(len(natural)*args.sides,len(indices))
        rr,cc=np.nonzero(local)
        row_parts.append((rr+offset).astype(np.int32));col_parts.append(indices[cc]);data_parts.append(local[rr,cc])
        offset+=len(natural)*args.sides;cells.append((int(ei),nodes,grads,center,volume,apex,phi))
    assert offset==row_count
    A=sparse.coo_matrix((np.concatenate(data_parts),(np.concatenate(row_parts),np.concatenate(col_parts))),shape=(row_count,count)).tocsr()
    del row_parts,col_parts,data_parts
    objective=external-capacity;scale=max(float(np.max(np.abs(objective))),1e-100)
    begin=time.monotonic();result=linprog(-objective/scale,A_ub=A,b_ub=np.zeros(row_count),bounds=[(-1.,1.)]*count,method='highs',
        options=dict(time_limit=args.time_limit,primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9))
    LP=dict(status=int(result.status),message=result.message,elapsed_seconds=time.monotonic()-begin,
        preparation_elapsed_seconds=begin-started,variables=count,inequalities=row_count,nonzeros=A.nnz,
        csr_storage_bytes=A.data.nbytes+A.indices.nbytes+A.indptr.nbytes,polygon_sides=args.sides,
        inner_polygon_inflation=inflation,sampling=args.sampling,samples_per_element=len(natural),
        maximum_tied_depth_derivative=maximum_depth_derivative,objective_scale=scale,
        best_scaled_work_excess=float(-result.fun) if result.status==0 else None)
    candidates=[]
    if result.status==0 and result.x is not None:
        zmin=float(np.min(xyz[used,2]));height=max(1.,float(np.max(xyz[used,2]))-zmin)
        can_add=all(c['direction'][2]==0 or all(xyz[n,2]==zmin for n in mesh['node_sets'][c['node_set']] if n in used_set) for c in problem['constraints'])
        for label,solution in [('raw',result.x),('rounded12',np.round(result.x,12)),('rounded10',np.round(result.x,10))]:
            for epsilon in ([1e-6,1e-5,0.] if can_add else [0.]):
                velocity=np.zeros_like(xyz)
                for n in used:
                    for axis in [0,2]:
                        index=variables.get((group_of[n],axis))
                        if index is not None:velocity[n,axis]=solution[index]
                    velocity[n,2]+=epsilon*(xyz[n,2]-zmin)/height
                for constraint in problem['constraints']:
                    axis=next(i for i,v in enumerate(constraint['direction']) if v!=0)
                    assert all(velocity[n,axis]==0 for n in mesh['node_sets'][constraint['node_set']] if n in used_set)
                valid=True;minimum=math.inf;cap=[];nonzero=0;failed=None
                for ei,nodes,grads,center,volume,apex,phi in cells:
                    local=velocity[nodes]-velocity[nodes[0]]
                    gradient=np.einsum('ni,snj->sij',local,grads);D=(gradient+gradient.swapaxes(1,2))/2
                    ratio=(1+math.sin(phi))/(1-math.sin(phi))
                    for point,tensor in enumerate(D):
                        size=float(np.max(np.abs(tensor)))
                        if size==0:continue
                        nonzero+=1;normal=tensor/size;eigen=np.linalg.eigvalsh(normal)
                        margin=float(np.trace(normal)+(ratio-1)*np.minimum(eigen,0).sum())/((1+ratio)*float(np.linalg.norm(normal)))
                        if not math.isfinite(margin) or margin<=1e-8:valid=False;failed=dict(element=ei,sample=point,relative_dual_margin=margin if math.isfinite(margin) else None);break
                        minimum=min(minimum,margin)
                    if not valid:break
                    cap.append(volume*apex*float(np.trace(local.T@center)))
                work=float(force@velocity.reshape(-1));total=math.fsum(cap)
                clearance=1e-5*max(float(np.abs(force*velocity.reshape(-1)).sum()),abs(total),1e-100)
                row=dict(kind=label,upward_epsilon_m=epsilon,proposal_samples_dual_valid=valid,
                    nonzero_sample_tensors=nonzero,minimum_relative_dual_margin=minimum if math.isfinite(minimum) else None,
                    external_work=work,integrated_capacity=total if valid else None,work_excess=work-total if valid else None,
                    first_failed_sample=failed,generator_clearance=clearance)
                candidates.append(row)
                if valid and nonzero and work-total>clearance:
                    output=copy.deepcopy(problem);output['virtual_displacement_m']=velocity.tolist();dump(args.output,output)
                    dump(meta,dict(proposed=True,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                        mesh_file_sha256=hashlib.sha256(args.mesh.read_bytes()).hexdigest(),
                        input_problem_file_sha256=hashlib.sha256(args.problem.read_bytes()).hexdigest(),
                        candidate_file_sha256=hashlib.sha256(args.output.read_bytes()).hexdigest(),
                        driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),LP=LP,candidates=candidates,
                        selected=row,independent_midside_quadratic_freedom=True,full_actual_GP_Rust_and_independent_verification_required=True,
                        not_a_nonlinear_FEM_or_Fs_or_stability_acceptance=True,p6_complete=False))
                    print(json.dumps(dict(proposed=True,selected=row,variables=count,sampling=args.sampling)),flush=True);return 0
    dump(meta,dict(proposed=False,LP=LP,candidates=candidates,
        no_proposal_is_not_equilibrium_or_stability_or_collapse_classification=True,p6_complete=False))
    print('No accepted quadratic proposal; no stability/collapse inference',flush=True);return 2
if __name__=='__main__':raise SystemExit(main())
