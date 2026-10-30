"""Bounded LP proposal for a plane-strain nodal separating field.

The polygon is an INNER approximation of the MC dual SOC. Every generated
candidate still needs the independent full-GP Rust verifier. No LP status is an
Fs/collapse classification. This reference generator does not change FEM inputs.
"""
from pathlib import Path
import argparse,copy,datetime,hashlib,json,math,os,time
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np
from scipy import sparse
from scipy.optimize import linprog

def strict(q):
    def constant(v):raise ValueError('nonfinite JSON '+v)
    def pairs(items):
        out={}
        for k,v in items:
            if k in out:raise ValueError('duplicate key '+k)
            out[k]=v
        return out
    return json.loads(Path(q).read_text('utf-8-sig'),parse_constant=constant,object_pairs_hook=pairs)
def finite(x):return type(x) in (int,float) and math.isfinite(x)
def dump(q,v):
    with q.open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
def main():
    a=argparse.ArgumentParser();a.add_argument('mesh',type=Path);a.add_argument('problem',type=Path);a.add_argument('output',type=Path)
    a.add_argument('--sides',type=int,default=16);a.add_argument('--time-limit',type=float,default=45.)
    args=a.parse_args();assert args.sides in [8,16,32,64] and 1.<=args.time_limit<=45.
    meta=args.output.with_suffix('.generation.json')
    assert not args.output.exists() and not meta.exists()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    mesh=strict(args.mesh);problem=strict(args.problem)
    assert problem['schema']=='soil-fem-mc-virtual-work-bound/1' and problem['water_policy']=='dry'
    xyz=np.asarray(mesh['nodes'],dtype=float);tets=np.asarray(mesh['tetrahedra'],dtype=int)
    active=problem['active_elements'];assert len(active)==len(tets) and all(type(x)is bool for x in active)
    selected=np.flatnonzero(active);assert 0<len(selected)<=20000
    assert xyz.ndim==2 and xyz.shape[1]==3 and np.isfinite(xyz).all()
    force=np.asarray(problem['physical_nodal_force_kn'],dtype=float)
    assert force.shape==(3*len(xyz),) and np.isfinite(force).all()
    high=mesh.get('tetrahedra10');high=np.asarray(high,dtype=int) if high else None
    if high is not None:assert high.shape==(len(tets),10) and np.array_equal(high[:,:4],tets)
    materials={}
    for m in problem['materials']:
        assert type(m['id'])is int and m['id']>0 and m['id'] not in materials
        p=m['parameters'];assert all(finite(p[k]) for k in ['cohesion_kpa','friction_deg','dilation_deg'])
        c=float(p['cohesion_kpa']);phi=math.radians(p['friction_deg']);psi=math.radians(p['dilation_deg'])
        assert c>=0. and 1e-8<phi<math.radians(89.) and 0.<=psi<=phi
        if p['flow']=='davis':
            beta=math.cos(phi)*math.cos(psi)/(1.-math.sin(phi)*math.sin(psi));c*=beta;phi=math.atan(beta*math.tan(phi))
        else:assert p['flow'] in ['associated','non_associated']
        materials[m['id']]=(c,phi)
    all_corners=set(tets.reshape(-1).tolist());mid_edges={};edges=[(0,1),(1,2),(2,0),(0,3),(1,3),(2,3)]
    geometry_scale=max(1.,float(np.max(np.abs(xyz))))
    if high is not None:
        for cell,row in zip(tets,high):
            for k,(i,j) in enumerate(edges):
                n=int(row[4+k]);pair=tuple(sorted((int(cell[i]),int(cell[j]))))
                assert n not in all_corners and (n not in mid_edges or mid_edges[n]==pair)
                assert np.max(np.abs(xyz[n]-(xyz[pair[0]]+xyz[pair[1]])/2.))<=1e-12*geometry_scale
                mid_edges[n]=pair
    used=set((high[selected] if high is not None else tets[selected]).reshape(-1).tolist())
    for n in range(len(xyz)):
        if n not in used:assert np.all(force[3*n:3*n+3]==0.)
    corners=sorted(set(tets[selected].reshape(-1).tolist()));groups={};node_group={}
    for n in corners:
        key=(float(xyz[n,0]).hex(),float(xyz[n,2]).hex())
        if key not in groups:groups[key]=len(groups)
        node_group[n]=groups[key]
    fixed=set()
    for constraint in problem['constraints']:
        assert constraint['value_m']==0. and constraint.get('gradient',[0.,0.,0.])==[0.,0.,0.]
        d=constraint['direction'];assert len(d)==3 and all(finite(v) for v in d)
        axes=[i for i,v in enumerate(d) if v!=0.];assert len(axes)==1
        axis=axes[0]
        if axis==1:continue
        for n in mesh['node_sets'][constraint['node_set']]:
            if n not in used:continue
            parents=[n] if n in node_group else mid_edges[n]
            # Overrestricting a constrained midpoint is conservative only for
            # proposing fields; the exact full verifier remains authoritative.
            fixed.update((node_group[k],axis) for k in parents)
    variables={(g,axis):i for i,(g,axis) in enumerate((g,axis) for g in range(len(groups)) for axis in [0,2] if (g,axis)not in fixed)}
    count=len(variables);assert count>0
    external=np.zeros(count);capacity=np.zeros(count);gradients=[]
    def scatter(row,group,axis,value):
        index=variables.get((group,axis))
        if index is not None:row[index]+=value
    for n in used:
        parents=[(n,1.)] if n in node_group else [(k,.5) for k in mid_edges[n]]
        for k,weight in parents:
            for axis in [0,2]:scatter(external,node_group[k],axis,force[3*n+axis]*weight)
    rr=[];cc=[];vv=[];rhs=[]
    inflation=(1.+1e-4)/math.cos(math.pi/args.sides)
    for element in selected:
        cell=tets[element];x=xyz[cell];J=(x[1:]-x[0]).T;det=float(np.linalg.det(J));assert det>0.
        inv=np.linalg.inv(J);grad=np.vstack([-inv.sum(axis=0),inv]);volume=det/6.
        c,phi=materials[mesh['material_ids'][element]];apex=c/math.tan(phi)
        trace=np.zeros(count);difference=np.zeros(count);shear=np.zeros(count);dy={}
        for n,g in zip(cell,grad):
            group=node_group[int(n)];dy[group]=dy.get(group,0.)+g[1]
            scatter(trace,group,0,g[0]);scatter(trace,group,2,g[2])
            scatter(difference,group,0,g[0]);scatter(difference,group,2,-g[2])
            scatter(shear,group,0,g[2]);scatter(shear,group,2,g[0])
        assert max(map(abs,dy.values()),default=0.)<=1e-12*max(1.,float(np.max(np.abs(grad)))),'element is not a plane-strain extrusion'
        capacity+=volume*apex*trace
        for k in range(args.sides):
            angle=2.*math.pi*k/args.sides
            row=math.sin(phi)*inflation*(math.cos(angle)*difference+math.sin(angle)*shear)-trace
            index=len(rhs);nz=np.flatnonzero(row);rr.extend([index]*len(nz));cc.extend(nz.tolist());vv.extend(row[nz].tolist());rhs.append(0.)
        gradients.append((int(element),cell,grad,volume,apex,phi))
    A=sparse.csr_matrix((vv,(rr,cc)),shape=(len(rhs),count));b=np.asarray(rhs)
    objective=external-capacity;scale=max(float(np.max(np.abs(objective))),1e-100);objective/=scale
    start=time.monotonic()
    result=linprog(-objective,A_ub=A,b_ub=b,bounds=[(-1.,1.)]*count,method='highs',
        options=dict(time_limit=args.time_limit,primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9))
    lp=dict(status=int(result.status),message=result.message,elapsed_seconds=time.monotonic()-start,
        variables=count,inequalities=len(rhs),nonzeros=A.nnz,polygon_sides=args.sides,
        inner_polygon_inflation=inflation,objective_scale=scale,
        best_scaled_work_excess=float(-result.fun) if result.status==0 else None)
    candidates=[]
    if result.status==0 and result.x is not None:
        z_reference=float(np.min(xyz[list(used),2]));height=max(1.,float(np.max(xyz[list(used),2]))-z_reference)
        can_add_upward=all(
            constraint['direction'][2]==0. or all(xyz[n,2]==z_reference for n in mesh['node_sets'][constraint['node_set']] if n in used)
            for constraint in problem['constraints'])
        variants=[]
        for base_label,solution in [('raw',result.x),('rounded12',np.round(result.x,12)),('rounded10',np.round(result.x,10)),('rounded8',np.round(result.x,8))]:
            variants.append((base_label,solution,0.))
            if can_add_upward:
                for epsilon in [1e-8,1e-7,1e-6,1e-5]:variants.append((base_label+'+upward'+str(epsilon),solution,epsilon))
        for label,solution,epsilon in variants:
            velocity=np.zeros((len(xyz),3))
            for n in corners:
                for axis in [0,2]:
                    index=variables.get((node_group[n],axis))
                    if index is not None:velocity[n,axis]=solution[index]
            for n,(i,j) in mid_edges.items():
                if n in used:velocity[n]=(velocity[i]+velocity[j])/2.
            if epsilon:
                for n in used:velocity[n,2]+=epsilon*(xyz[n,2]-z_reference)/height
            valid=True;minimum=math.inf;cap=[];nonzero=0;failed_element=None;failed_margin=None
            for element,cell,grad,volume,apex,phi in gradients:
                local=velocity[cell]-velocity[cell[0]];gradient=local.T@grad;D=(gradient+gradient.T)/2.
                size=float(np.max(np.abs(D)))
                if size==0.:cap.append(0.);continue
                normal=D/size;eigen=np.linalg.eigvalsh(normal);R=(1.+math.sin(phi))/(1.-math.sin(phi))
                margin=float(np.trace(normal)+(R-1.)*np.minimum(eigen,0.).sum())/((1.+R)*float(np.linalg.norm(normal)))
                if not math.isfinite(margin) or margin<=1e-8:
                    valid=False;failed_element=element;failed_margin=margin if math.isfinite(margin) else None;break
                minimum=min(minimum,margin);cap.append(volume*apex*float(np.trace(D)));nonzero+=1
            work=float(force@velocity.reshape(-1));total=math.fsum(cap)
            clearance=1e-5*max(float(np.abs(force*velocity.reshape(-1)).sum()),abs(total),1e-100)
            candidate=dict(kind=label,linearized_all_cells_dual_valid=valid,nonzero_strain_elements=nonzero,
                minimum_relative_dual_margin=minimum if math.isfinite(minimum) else None,
                external_work=work,integrated_capacity=total if valid else None,work_excess=work-total if valid else None,
                checked_prefix_capacity_not_a_full_bound=None if valid else total,
                first_failed_element=failed_element,first_failed_relative_dual_margin=failed_margin,
                additional_upward_virtual_displacement_m=epsilon,
                admissible_global_upward_field_preserves_all_axis_constraints=can_add_upward,
                generator_clearance=clearance)
            candidates.append(candidate)
            if valid and nonzero and work-total>clearance:
                for constraint in problem['constraints']:
                    axis=next(i for i,v in enumerate(constraint['direction']) if v!=0.)
                    assert all(velocity[n,axis]==0. for n in mesh['node_sets'][constraint['node_set']] if n in used)
                output=copy.deepcopy(problem);output['virtual_displacement_m']=velocity.tolist()
                dump(args.output,output)
                dump(meta,dict(proposed=True,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    mesh_file_sha256=hashlib.sha256(args.mesh.read_bytes()).hexdigest(),
                    input_problem_file_sha256=hashlib.sha256(args.problem.read_bytes()).hexdigest(),
                    candidate_file_sha256=hashlib.sha256(args.output.read_bytes()).hexdigest(),
                    driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),LP=lp,candidates=candidates,
                    selected_candidate=label,full_actual_GP_Rust_verification_required=True,
                    not_a_nonlinear_FEM_solution_or_Fs_or_stability_acceptance=True,p6_complete=False))
                print(json.dumps(dict(proposed=True,selected=label,work_excess=work-total,requires_full_Rust_GP_verification=True)),flush=True)
                return 0
    dump(meta,dict(proposed=False,LP=lp,candidates=candidates,not_found_is_not_physical_exclusion_or_solver_collapse=True,
        original_inputs_unchanged=True,p6_complete=False))
    print('No accepted proposal; no collapse or Fs classification',flush=True)
    return 2
if __name__=='__main__':raise SystemExit(main())
