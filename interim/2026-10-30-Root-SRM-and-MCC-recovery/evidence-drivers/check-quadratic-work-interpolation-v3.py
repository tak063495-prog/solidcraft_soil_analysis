"""Independent polynomial interpolation and exact moments for the proposer."""
from pathlib import Path
import hashlib,importlib.util,json,math,os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np
from scipy.special import roots_jacobi
p=Path(__file__).resolve().parent;source=p/'generate-plane-strain-quadratic-work-v3.py';out=p/'slope1200-durable-SRM-v1/quadratic-proposal-interpolation-v3';out.mkdir()
def sha(q):
    with q.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
before=sha(source)
spec=importlib.util.spec_from_file_location('fixed_quadratic_proposer',source);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
corners=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
nodes=np.vstack([corners,[(corners[i]+corners[j])/2 for i,j in [(0,1),(1,2),(2,0),(0,3),(1,3),(2,3)]]])
def monomial(q):
    r,s,t=q;return np.array([1,r,s,t,r*r,s*s,t*t,r*s,s*t,r*t])
def dmonomial(q):
    r,s,t=q
    return np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[2*r,0,0],[0,2*s,0],[0,0,2*t],[s,r,0],[0,t,s],[t,0,r]])
coeff=np.linalg.solve(np.array([monomial(q) for q in nodes]),np.eye(10))
rng=np.random.default_rng(20261030);maximum_derivative_error=0.;maximum_affine_strain_error=0.;maximum_physical_gradient_error=0.
samples=np.vstack([corners,rng.uniform(-.3,1.3,size=(80,3))])
for q in samples:
    actual=m.derivatives(q);reference=coeff.T@dmonomial(q)
    maximum_derivative_error=max(maximum_derivative_error,float(np.max(np.abs(actual-reference))))
    assert np.max(np.abs(actual-reference))<=2e-14
    values=rng.normal(size=(10,3));at_vertices=np.array([values.T@m.derivatives(v) for v in corners]);bary=np.r_[1-sum(q),q]
    error=float(np.max(np.abs(np.einsum('s,sij->ij',bary,at_vertices)-values.T@actual)))
    maximum_affine_strain_error=max(maximum_affine_strain_error,error);assert error<5e-14
for _ in range(20):
    J=np.eye(3)+rng.uniform(-.2,.2,size=(3,3));shift=rng.uniform(-2,2,size=3);physical=nodes@J.T+shift
    # Genuine quadratic x/z fields. A corner-average midpoint prescription cannot reproduce these.
    x,y,z=physical.T;values=np.c_[2*x*x+3*x*z-z*z+.4*x-.5*z,np.zeros(10),-x*x+2*x*z+4*z*z+.7*x+.2*z]
    assert max(np.linalg.norm(values[4+k]-(values[i]+values[j])/2) for k,(i,j) in enumerate(m.EDGES))>1e-4
    for q in samples:
        x,y,z=J@q+shift
        expected=np.array([[4*x+3*z+.4,0,3*x-2*z-.5],[0,0,0],[-2*x+2*z+.7,0,2*x+8*z+.2]])
        got=values.T@(m.derivatives(q)@np.linalg.inv(J));error=float(np.max(np.abs(got-expected)))
        maximum_physical_gradient_error=max(maximum_physical_gradient_error,error);assert error<1e-12
quadrature=m.sampling('gauss');weighted=[roots_jacobi(3,a,0) for a in [2,1,0]]
weights=np.array([a*b*c/64 for a in weighted[0][1] for b in weighted[1][1] for c in weighted[2][1]])
moments=0;maximum_moment_error=0.
for i in range(6):
    for j in range(6-i):
        for k in range(6-i-j):
            exact=math.factorial(i)*math.factorial(j)*math.factorial(k)/math.factorial(i+j+k+3)
            actual=float(np.dot(weights,quadrature[:,0]**i*quadrature[:,1]**j*quadrature[:,2]**k))
            error=abs(actual-exact);maximum_moment_error=max(maximum_moment_error,error);assert error<2e-16;moments+=1
assert moments==56
# Independent spectral dual check on a grid between polygon normals.
dual_checks=0;min_dual=math.inf
for sides in [16,32,64,128]:
    normals=np.arange(sides)*2*math.pi/sides;inflation=(1+1e-4)/math.cos(math.pi/sides)
    for theta in np.linspace(0,2*math.pi,257):
        for phi in [.1,.35,.8]:
            difference,shear=math.cos(theta),math.sin(theta)
            trace=math.sin(phi)*inflation*max(math.cos(a)*difference+math.sin(a)*shear for a in normals)
            D=np.array([[(trace+difference)/2,0,shear/2],[0,0,0],[shear/2,0,(trace-difference)/2]])
            eigen=np.linalg.eigvalsh(D);ratio=(1+math.sin(phi))/(1-math.sin(phi))
            dual=float(np.sum(eigen[eigen>=0])+ratio*np.sum(eigen[eigen<0]))
            min_dual=min(min_dual,dual);assert dual>=-1e-14;dual_checks+=1
# The condensed rows must also reproduce a genuinely quadratic field on each
# distinct projected triangle, using a separate analytical derivative.
mesh=json.loads((p/'srm-general-virtual-work-v1/slope1200-LP-v1/mesh.json').read_text('utf-8-sig'))
xyz=np.asarray(mesh['nodes']);groups={};group_of={}
for n,q in enumerate(xyz):
    key=(float(q[0]).hex(),float(q[2]).hex())
    if key not in groups:groups[key]=len(groups)
    group_of[n]=groups[key]
variables={(g,a):i for i,(g,a) in enumerate((g,a) for g in range(len(groups)) for a in [0,2])}
values=np.zeros(len(variables))
for n,(x,y,z) in enumerate(xyz):
    values[variables[(group_of[n],0)]]=2*x*x+3*x*z-z*z+.4*x-.5*z
    values[variables[(group_of[n],2)]]=-x*x+2*x*z+4*z*z+.7*x+.2*z
maximum_condensed_polynomial_error=0.
for cell,nodes10 in zip(mesh['tetrahedra'],mesh['tetrahedra10']):
    signature,indices,T,D,S,mapping=m.projected_rows(xyz,np.asarray(nodes10),np.asarray(cell),group_of,variables)
    reps={group_of[n]:n for n in cell}
    reference=[]
    for group in signature:
        x,y,z=xyz[reps[group]]
        xx=4*x+3*z+.4;zz=2*x+8*z+.2;xz=(3*x-2*z-.5)+(-2*x+2*z+.7)
        reference.append([xx+zz,xx-zz,xz])
    actual=np.c_[T@values[indices],D@values[indices],S@values[indices]]
    error=float(np.max(np.abs(actual-reference)));maximum_condensed_polynomial_error=max(maximum_condensed_polynomial_error,error)
    assert error<2e-11

assert sha(source)==before
proof=dict(passed=True,condensed_actual_all1200_elements_analytic_quadratic_field_checked=True,
    maximum_condensed_polynomial_error=maximum_condensed_polynomial_error,source_sha256=before,independent_Vandermonde_shape_gradient_samples=len(samples),
    maximum_derivative_error=maximum_derivative_error,
    independent_physical_quadratic_field_samples=20*len(samples),maximum_physical_gradient_error=maximum_physical_gradient_error,
    affine_strain_vertex_interpolation_maximum_error=maximum_affine_strain_error,
    all56_exact_Jacobi_tetrahedral_moments_checked=True,maximum_moment_error=maximum_moment_error,
    independently_eigenchecked_polygon_inner_dual_cases=dual_checks,minimum_dual_support=min_dual,
    actual_full_GP_Rust_and_independent_work_acceptance_still_separate=True,p6_complete=False)
with (out/'terminal.json').open('x',encoding='utf-8') as f:json.dump(proof,f,indent=2,allow_nan=False)
print(json.dumps(proof),flush=True)
