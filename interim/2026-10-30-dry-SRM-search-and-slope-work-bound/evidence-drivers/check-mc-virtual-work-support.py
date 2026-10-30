"""Independent LP and rotated stress samples for a proposed general MC witness.

This derives/checks a support function, not a slope solver or accepted Rust gate.
"""
from pathlib import Path
import datetime, hashlib, json, math, os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1'
import numpy as np
from scipy.optimize import linprog
p=Path(__file__).resolve().parent
out=p/'srm-general-virtual-work-v1/independent-support-v1';out.mkdir(parents=True)
assert not (p/'pause-requested').exists()
def dump(name,v):
    with (out/name).open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
dump('source-binding.json',dict(driver=str(Path(__file__).resolve()),
    driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    numpy_version=np.__version__,scipy_version=__import__('scipy').__version__,
    MC_primary_source='https://docs.itascacg.com/itasca950/common/models/mohr/doc/modelmohr.html',
    convex_support_and_dual_cone_primary_source='https://web.stanford.edu/~boyd/cvxbook/bv_cvxbook.pdf',
    source_conventions_differ_Itasca_compression_negative_here_compression_positive=True,
    formula_is_our_derivation_not_a_quote_or_implementation_from_the_sources=True,
    no_source_code_copied=True,not_a_Rust_or_slope_or_P6_acceptance=True))
rng=np.random.default_rng(721630)
rows=[];samples=0;max_violation=-math.inf;negative_control=0.
for phi_deg in [5.,15.,30.,60.]:
    phi=math.radians(phi_deg);s=math.sin(phi);ratio=(1.+s)/(1.-s)
    for cohesion in [0.,10.,80.]:
        apex=cohesion/math.tan(phi)
        A=[]
        for i in range(3):
            for j in range(3):
                if i==j:continue
                row=np.zeros(3);row[i]=1.-s;row[j]=-(1.+s);A.append(row)
        b=np.full(6,2.*cohesion*math.cos(phi))
        for label,diagonal,finite in [
            ('interior_rank_one',[-1.,0.,ratio+.2],True),
            ('interior_three_principals',[-1.,-.25,1.25*ratio+.2],True),
            ('boundary_associated',[-(1.-s),0.,1.+s],True),
            ('outside_polar',[-1.,0.,ratio-.2],False),
            ('compression_direction',[-1.,-1.,-1.],False),
        ]:
            d=np.asarray(diagonal)
            result=linprog(d,A_ub=A,b_ub=b,bounds=[(None,None)]*3,method='highs',
                options=dict(primal_feasibility_tolerance=1e-9,dual_feasibility_tolerance=1e-9))
            dual_margin=float(np.sum(np.maximum(d,0.))+ratio*np.sum(np.minimum(d,0.)))
            predicted=apex*float(d.sum()) if finite else None
            if finite:
                assert result.status==0,(phi_deg,cohesion,label,result.message)
                assert abs(-result.fun-predicted)<=1e-8*max(1.,abs(predicted))
                # Test the full rotated 3D tensor inequality against admissible
                # stress matrices at many hydrostatic scales, not just diagonals.
                q,_=np.linalg.qr(rng.normal(size=(3,3)))
                strain=q@np.diag(d)@q.T
                for k in range(32):
                    rotation,_=np.linalg.qr(rng.normal(size=(3,3)))
                    alpha=10.**rng.uniform(-3.,6.)
                    eigen=alpha*(1.+(ratio-1.)*rng.random(3))
                    stress=-apex*np.eye(3)+rotation@np.diag(eigen)@rotation.T
                    actual=float(-np.sum(stress*strain))
                    gap=actual-predicted
                    max_violation=max(max_violation,gap)
                    assert gap<=1e-7*max(1.,abs(actual),abs(predicted))
                    samples+=1
                wrong_sign=float(np.sum((100000.*np.eye(3))*strain))
                negative_control=max(negative_control,wrong_sign-predicted)
            else:
                assert result.status==3,(phi_deg,cohesion,label,result.message)
            rows.append(dict(phi_deg=phi_deg,cohesion_kpa=cohesion,case=label,
                strain_principals=diagonal,dual_margin=dual_margin,LP_status=int(result.status),
                predicted_support=predicted,independent_LP_support=float(-result.fun) if finite else None))
assert negative_control>1000.
# Tresca is a separate singular case: tiny nonzero trace is not silently zero.
tresca=[]
for d,finite in [([-1.,0.,1.],True),([-1.,0.,1.001],False),([-1.,0.,.999],False)]:
    A=[]
    for i in range(3):
        for j in range(3):
            if i==j:continue
            row=np.zeros(3);row[i]=1.;row[j]=-1.;A.append(row)
    r=linprog(d,A_ub=A,b_ub=np.full(6,20.),bounds=[(None,None)]*3,method='highs')
    assert r.status==(0 if finite else 3)
    if finite:assert abs(-r.fun-20.)<1e-10
    tresca.append(dict(strain_principals=d,LP_status=int(r.status),support=float(-r.fun) if finite else None))
dump('terminal.json',dict(passed=True,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    independent_LP_cases=len(rows),rotated_full_tensor_samples=samples,
    all_finite_supports_and_unbounded_directions_matched=True,
    maximum_sample_work_minus_bound=max_violation,
    deliberately_wrong_stress_sign_violates_the_bound_by=negative_control,
    all3_Tresca_directions_matched=True,nonzero_trace_never_treated_as_incompressible=True,
    derived_frictional_support='h(D)=c*cot(phi)*tr(D), finite iff sumpositive(eig(D))+R*sumnegative(eig(D))>=0; R=(1+sin(phi))/(1-sin(phi))',
    next_step='Implement nodal virtual-work verifier with actual FE B, positive weights, exact supported-component zero checks and bound endpoint physical loads. Then independently validate distributed-support FEM and actual slope mechanisms.',
    Rust_implementation_or_ordinary_slope_witness_not_yet_accepted=True,
    no_new_FEM_solution_or_physical_slope_upper_or_safety_factor_from_this_test=True,
    p6_complete=False,rows=rows,tresca=tresca))
print(json.dumps(dict(passed=True,independent_LP_cases=len(rows),rotated_tensor_samples=samples,
    Tresca_cases=3,not_yet_Rust_or_slope_acceptance=True)),flush=True)
