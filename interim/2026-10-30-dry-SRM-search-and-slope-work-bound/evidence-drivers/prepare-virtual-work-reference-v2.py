"""Match quadrature families independently; retain the valid Legendre reference."""
from pathlib import Path
import hashlib,json
p=Path(__file__).resolve().parent
source=(p/'check-mc-virtual-work-fem.py').read_text('utf-8')
needle="out=p/'srm-general-virtual-work-v1/independent-FEM-v1'"
assert source.count(needle)==1
source=source.replace(needle,"out=p/'srm-general-virtual-work-v1/independent-FEM-v2'",1)
old="""gx,gw=np.polynomial.legendre.leggauss(3);gx=(gx+1.)/2.;gw=gw/2.
quadrature=[]
for i,j,k in itertools.product(range(3),repeat=3):
    a,b,c=gx[[i,j,k]]
    quadrature.append(([a,(1.-a)*b,(1.-a)*(1.-b)*c],gw[i]*gw[j]*gw[k]*(1.-a)**2*(1.-b)))
assert abs(math.fsum(w for _,w in quadrature)-1./6.)<1e-15"""
new="""from scipy.special import roots_jacobi
rules=[]
for alpha in [2.,1.,0.]:
    x,w=roots_jacobi(3,alpha,0.)
    rules.append(((x+1.)/2.,w/2.**(alpha+1.)))
quadrature=[]
for i,j,k in itertools.product(range(3),repeat=3):
    a,b,c=rules[0][0][i],rules[1][0][j],rules[2][0][k]
    quadrature.append(([a,(1.-a)*b,(1.-a)*(1.-b)*c],rules[0][1][i]*rules[1][1][j]*rules[2][1][k]))
assert abs(math.fsum(w for _,w in quadrature)-1./6.)<1e-15
moment_count=0
for i in range(6):
    for j in range(6-i):
        for k in range(6-i-j):
            actual=math.fsum(w*r[0]**i*r[1]**j*r[2]**k for r,w in quadrature)
            exact=math.factorial(i)*math.factorial(j)*math.factorial(k)/math.factorial(i+j+k+3)
            assert abs(actual-exact)<2e-15
            moment_count+=1
assert moment_count==56"""
assert source.count(old)==1
source=source.replace(old,new,1)
needle="dump('terminal.json',dict(passed=True,actual_cases=6,orders=[4,10],distributed_supports=True,"
assert source.count(needle)==1
source=source.replace(needle,needle+"\n    quadrature_family='independent scipy.special shifted Gauss-Jacobi alpha2,1,0',\n    all56_monomials_through_total_degree5_matched=True,original_Legendre_reference_v1_kept=True,",1)
with (p/'check-mc-virtual-work-fem-v2.py').open('x',encoding='utf-8',newline='\n') as f:f.write(source)
out=p/'srm-general-virtual-work-v1/reference-family-correction-v1';out.mkdir()
with (out/'preparation.json').open('x',encoding='utf-8') as f:
    json.dump(dict(prepared=True,actual_650_driver_session55398_exited1=True,
        actual_Rust12_new_plus11_prior_search_plus9_prior_mean_passed=True,
        first3_Tet4_CLI_reference_comparisons_passed=True,actual_Tet10_CLI_exit_code0=True,
        failed_comparison='point0 weight0.008770474929651058 vs0.014972747367083927',
        cause='Rust uses weighted shifted Gauss-Jacobi2,1,0; independent v1 used unweighted Legendre3 with explicit Duffy Jacobian. Both integrate the global patch but their individual GP sets differ.',
        no_Rust_sources_or_integration_rule_or_tolerances_changed=True,
        independent_SciPy_Jacobi_generator_does_not_copy_Rust_recurrence_or_tables=True,
        all56_exact_reference_tetrahedron_moments_required=True,
        original_v1_reference_and650_failed_harness_logs_retained=True,
        new_reference_driver_sha256=hashlib.sha256(source.encode()).hexdigest(),p6_complete=False),f,indent=2)
print('Prepared independent matching Jacobi reference; original records unchanged',flush=True)
