"""Actual frozen650 CLI plus independent all-GP slope work certificate audit.

This checks a discrete MC equilibrium exclusion, not a continuum collapse Fs.
The independent calculation uses own shape derivatives and SciPy quadrature.
"""
from pathlib import Path
import datetime,hashlib,itertools,json,math,os,subprocess,sys,time,traceback
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np
from scipy.special import roots_jacobi
from virtual_work_references import point_stream_hash
p=Path(__file__).resolve().parent;candidate=p/'srm-general-virtual-work-v1'
out=candidate/'slope1200-virtual-audit-v1';out.mkdir()
read=lambda q:json.loads(Path(q).read_text('utf-8-sig'))
def sha(q):
    with Path(q).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(name,v):
    with (out/name).open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
def active():
    if (p/'pause-requested').exists():raise SystemExit(2)
active()
solver=candidate/'quality-v1/soil-fem-3d.exe'
meshfile=candidate/'slope1200-LP-v1/mesh.json'
problemfile=candidate/'slope1200-LP-v2/F2p0/candidate.json'
proof=candidate/'slope1200-LP-v1/physical-load-proof.json'
frozen=candidate/'quality-v1/source-snapshot.json';sources=read(frozen)
binding={str(q):sha(q) for q in [solver,meshfile,problemfile,proof,frozen,
    candidate/'slope1200-LP-v2/F2p0/candidate.generation.json',Path(__file__)]}
assert len(sources)==650 and all(sha(candidate/'runtime'/n)==h for n,h in sources.items())
assert read(proof)['passed'] and read(candidate/'independent-FEM-v2/terminal.json')['passed']
assert read(candidate/'slope1200-LP-v2/F2p0/candidate.generation.json')['candidate_file_sha256']==sha(problemfile)
dump('driver-child.json',dict(pid=os.getpid(),executable=sys._base_executable,
    script=str(Path(__file__).resolve()),created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
dump('protocol.json',dict(prepared=True,bindings_sha256=binding,absolute_strength_factor=2.,
    actual_1200_Tet10_all32400_GP_required=True,independent_shape_derivatives_and_SciPy_Jacobi=True,
    raw_tensor_error_budget='4096 eps cond(J) geometry_ratio Frobenius(abs(dV)^T abs(gradN))',
    eig_comparison_uses_Weyl_Frobenius_bound=True,
    dual_cone_must_remain_strict_after_Lipschitz_error_bound=True,
    existing_six_cube_comparisons_keep_original2e_minus12_limits=True,
    original5percent_displacement_acceptance_remains_FAILED=True,
    numerical_Newton_failure_never_used_as_an_upper=True,p6_complete=False))

def main():
    active();resultfile=out/'actual-Rust-all-GP-result.json'
    command=[str(solver),'mc-virtual-work-bound',str(meshfile),str(problemfile),str(resultfile)]
    start=time.monotonic()
    with (out/'Rust-CLI.log').open('x',encoding='utf-8') as log:
        child=subprocess.Popen(command,cwd=candidate/'runtime',stdout=log,stderr=subprocess.STDOUT,
            creationflags=subprocess.CREATE_NO_WINDOW|subprocess.BELOW_NORMAL_PRIORITY_CLASS)
        dump('Rust-CLI-child.json',dict(pid=child.pid,parent_pid=os.getpid(),executable=command[0],
            arguments=command[1:],created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
        code=child.wait()
    dump('Rust-CLI-terminal.json',dict(exit_code=code,elapsed_seconds=time.monotonic()-start))
    assert code==0
    raw=resultfile.read_text('utf-8');result=json.loads(raw)
    assert result['assessment']['status']=='evaluated' and result['complete_point_audit']
    cert=result['assessment']['certificate'];points=result['point_checks']
    assert cert['integration_points']==len(points)==32400 and cert['active_elements']==1200
    assert cert['excludes_admissible_discrete_equilibrium']
    assert not cert['continuum_collapse_or_stability_certified']
    digest,count=point_stream_hash(raw)
    assert digest==cert['point_checks_sha256'] and count==32400
    del raw
    mesh=read(meshfile);problem=read(problemfile)
    xyz=np.asarray(mesh['nodes']);v=np.asarray(problem['virtual_displacement_m'])
    force=np.asarray(problem['physical_nodal_force_kn'])
    assert v.shape==xyz.shape==(2583,3) and len(force)==7749
    assert len(mesh['tetrahedra10'])==1200 and all(problem['active_elements'])
    assert problem['water_policy']=='dry' and len(problem['materials'])==1
    params=problem['materials'][0]['parameters'];phi=math.radians(params['friction_deg']);c=params['cohesion_kpa']
    assert abs(c-5.)<1e-14 and abs(math.tan(phi)-math.tan(math.radians(20.))/2.)<1e-14
    assert params['flow']=='non_associated'
    ratio=(1.+math.sin(phi))/(1.-math.sin(phi));intercept=c/math.tan(phi)
    constraints=0
    for constraint in problem['constraints']:
        assert constraint['value_m']==0. and constraint.get('gradient',[0.,0.,0.])==[0.,0.,0.]
        axes=np.flatnonzero(constraint['direction']);assert len(axes)==1
        for n in mesh['node_sets'][constraint['node_set']]:
            assert v[n,axes[0]]==0.;constraints+=1
    assert cert['all_constrained_components_exactly_zero'] and cert['constrained_components_checked']==constraints
    rules=[]
    for alpha in [2.,1.,0.]:
        q,w=roots_jacobi(3,alpha,0.);rules.append(((q+1.)/2.,w/2.**(alpha+1.)))
    quadrature=[]
    for i,j,k in itertools.product(range(3),repeat=3):
        a,b,d=rules[0][0][i],rules[1][0][j],rules[2][0][k]
        quadrature.append((np.array([a,(1.-a)*b,(1.-a)*(1.-b)*d]),rules[0][1][i]*rules[1][1][j]*rules[2][1][k]))
    moments=0
    for i in range(6):
        for j in range(6-i):
            for k in range(6-i-j):
                actual=math.fsum(w*r[0]**i*r[1]**j*r[2]**k for r,w in quadrature)
                exact=math.factorial(i)*math.factorial(j)*math.factorial(k)/math.factorial(i+j+k+3)
                assert abs(actual-exact)<2e-15;moments+=1
    assert moments==56
    grad_l=np.array([[-1.,-1.,-1.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    edges=[(0,1),(1,2),(2,0),(0,3),(1,3),(2,3)]
    eps=np.finfo(float).eps;capacities=[];errors=[];weights=[]
    minimum_clear_margin=math.inf;maximum_tensor_fraction=0.;maximum_eigen_fraction=0.
    maximum_weight_scaled_error=0.;max_tensor_budget=0.;negative_control_checked=False
    analytic_force=np.zeros(len(force));index=0
    for element,cell in enumerate(mesh['tetrahedra10']):
        active();local_xyz=xyz[cell];local=local_xyz-local_xyz[0];dv=v[cell]-v[cell[0]]
        det_linear=float(np.linalg.det((local[1:4]).T));assert det_linear>0.
        volume=det_linear/6.
        for n in cell[:4]:analytic_force[3*n+2]+=volume
        for n in cell[4:]:analytic_force[3*n+2]-=4.*volume
        for gi,(rst,weight) in enumerate(quadrature):
            actual=points[index];index+=1
            assert actual['element_index']==element and actual['local_index']==gi
            assert actual['material_id']==1 and actual['accepted_finite_support']
            l=np.r_[1.-sum(rst),rst]
            shape=np.r_[l*(2.*l-1.),[4.*l[i]*l[j] for i,j in edges]]
            dn=np.vstack([(4.*l[:,None]-1.)*grad_l,
                [4.*(l[j]*grad_l[i]+l[i]*grad_l[j]) for i,j in edges]])
            J=local.T@dn;det=float(np.linalg.det(J));assert det>0.
            inverse=np.linalg.inv(J);grad=dn@inverse
            G=dv.T@grad;D=(G+G.T)/2.
            cond=float(np.linalg.cond(J));geom=max(1.,float(np.max(np.abs(local_xyz)))*float(np.linalg.norm(inverse)))
            operator=float(np.linalg.norm(np.abs(dv).T@np.abs(grad)))
            tensor_budget=4096.*eps*max(1.,cond)*geom*max(operator,np.finfo(float).tiny)
            max_tensor_budget=max(max_tensor_budget,tensor_budget)
            xx,yy,zz,xy,yz,xz=actual['normalized_strain_tensor_tension_positive']
            rust=np.array([[xx,xy,xz],[xy,yy,yz],[xz,yz,zz]])*actual['strain_tensor_scale']
            error=float(np.linalg.norm(rust-D));assert error<=tensor_budget,(element,gi,error,tensor_budget)
            maximum_tensor_fraction=max(maximum_tensor_fraction,error/tensor_budget)
            eig=np.linalg.eigvalsh(D)[::-1]
            reig=np.asarray(actual['normalized_principal_strains_descending'])*actual['strain_tensor_scale']
            eigen_error=float(np.max(np.abs(reig-eig)));assert eigen_error<=tensor_budget
            maximum_eigen_fraction=max(maximum_eigen_fraction,eigen_error/tensor_budget)
            # Weyl + Cauchy-Schwarz give Lipschitz constant sqrt(3)*R.
            dual=float(np.trace(D)+(ratio-1.)*np.minimum(eig,0.).sum())
            lower=dual-math.sqrt(3.)*ratio*tensor_budget
            required=1e-8*(1.+ratio)*(float(np.linalg.norm(D))+tensor_budget)
            assert lower>required,(element,gi,lower,required)
            minimum_clear_margin=min(minimum_clear_margin,lower-required)
            w=weight*det
            we=4096.*eps*max(1.,cond)*geom*max(abs(w),np.finfo(float).tiny)
            assert abs(actual['weight_m3']-w)<=we
            maximum_weight_scaled_error=max(maximum_weight_scaled_error,abs(actual['weight_m3']-w)/max(1.,abs(w)))
            assert np.max(np.abs(np.asarray(actual['natural_coordinates'])-rst))<2e-12
            assert np.max(np.abs(np.asarray(actual['position_m'])-shape@local_xyz))<2e-12*max(1.,float(np.max(np.abs(local_xyz))))
            density=intercept*float(np.trace(D));density_budget=intercept*math.sqrt(3.)*tensor_budget+4096.*eps*abs(density)
            assert abs(actual['support_density_kpa']-density)<=density_budget
            assert abs(actual['effective_cohesion_kpa']-c)<2e-12
            assert abs(actual['effective_friction_radians']-phi)<2e-12
            capacities.append(w*density);weights.append(w)
            errors.append(abs(w)*density_budget+abs(density)*we+we*density_budget)
            if not negative_control_checked:
                corrupted=rust.copy();corrupted[0,0]+=10.*tensor_budget
                assert np.linalg.norm(corrupted-D)>tensor_budget
                bad=-np.eye(3);bad_eig=np.linalg.eigvalsh(bad)
                assert float(np.trace(bad)+(ratio-1.)*np.minimum(bad_eig,0.).sum())<0.
                negative_control_checked=True
    assert index==32400
    work=math.fsum(float(a)*float(b) for a,b in zip(force,v.reshape(-1)))
    abs_work=math.fsum(abs(float(a)*float(b)) for a,b in zip(force,v.reshape(-1)))
    force_error=float(np.max(np.abs(force-analytic_force)))/max(1.,float(np.max(np.abs(force))))
    assert force_error<1e-11
    analytic_work=math.fsum(float(a)*float(b) for a,b in zip(analytic_force,v.reshape(-1)))
    force_work_budget=math.fsum(abs(float(a-b)*float(w)) for a,b,w in zip(force,analytic_force,v.reshape(-1)))
    work_budget=4096.*eps*max(1.,abs_work)
    assert abs(work-analytic_work)<=force_work_budget+work_budget
    capacity=math.fsum(capacities);capacity_budget=math.fsum(errors)+4096.*eps*math.fsum(map(abs,capacities))
    assert abs(work-cert['external_virtual_work_kn_m'])<=work_budget
    assert abs(capacity-cert['maximum_admissible_internal_virtual_work_kn_m'])<=capacity_budget
    excess=work-capacity
    assert excess>cert['numerical_work_clearance_kn_m']+capacity_budget+work_budget+force_work_budget
    assert abs(math.fsum(weights)-220.)<1e-10
    assert all(sha(q)==h for q,h in binding.items())
    assert all(sha(candidate/'runtime'/n)==h for n,h in sources.items())
    dump('terminal.json',dict(passed=True,actual_Rust_exit_code=0,absolute_strength_factor=2.,
        elements=1200,nodes=2583,all_GP_checked=32400,framed_point_sha256=digest,
        all56_exact_quadrature_moments_matched=True,all_axis_supports_exact_zero=True,
        independent_raw_tensor_and_Weyl_eigen_comparisons_passed=True,
        minimum_independent_dual_lower_minus_required_margin=minimum_clear_margin,
        maximum_tensor_error_fraction_of_predeclared_bound=maximum_tensor_fraction,
        maximum_eigen_error_fraction_of_predeclared_bound=maximum_eigen_fraction,
        maximum_raw_tensor_error_budget=max_tensor_budget,maximum_weight_scaled_error=maximum_weight_scaled_error,
        independent_external_work_kn_m=work,independent_capacity_kn_m=capacity,
        independent_work_excess_kn_m=excess,capacity_error_budget_kn_m=capacity_budget,
        work_error_budget_kn_m=work_budget,source_force_work_difference_budget_kn_m=force_work_budget,
        actual_Rust_clearance_kn_m=cert['numerical_work_clearance_kn_m'],
        independent_source_load_component_error=force_error,
        corrupted_tensor_and_dual_negative_controls_detected=negative_control_checked,
        excludes_admissible_discrete_MC_equilibrium_at_F2=True,
        continuum_collapse_or_stability_certified=False,exact_collapse_Fs_not_established=True,
        original_displacement5percent_acceptance_remains_FAILED=True,
        original_inputs_result_and_frozen650_sources_unchanged=True,
        result_sha256=sha(resultfile),elapsed_seconds=time.monotonic()-start,p6_complete=False))
    print('PASS actual slope1200 all32400 GP discrete MC virtual-work exclusion F2; exact Fs/P6 remain',flush=True)
try:main()
except Exception as exc:
    if not (out/'terminal.json').exists():dump('terminal.json',dict(passed=False,error=repr(exc),traceback=traceback.format_exc(),p6_complete=False))
    raise
