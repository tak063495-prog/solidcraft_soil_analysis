"""Independent dense 3D FEM/virtual-work reference with distributed supports.

Own Tet4/Tet10 shape derivatives and quadrature; no Rust solver code is called.
"""
from pathlib import Path
import copy, hashlib, itertools, json, math, os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1';os.environ['MKL_NUM_THREADS']='1'
import numpy as np
p=Path(__file__).resolve().parent
out=p/'srm-general-virtual-work-v1/independent-FEM-v2';out.mkdir(parents=True)
assert not (p/'pause-requested').exists()
def dump(name,v):
    with (out/name).open('x',encoding='utf-8') as f:json.dump(v,f,indent=2,allow_nan=False)
base=p/'srm-search-development-v1/runtime/examples/srm-mean-stress-mesh.json'
linear=json.loads(base.read_text('utf-8'))
edges=[(0,1),(1,2),(2,0),(0,3),(1,3),(2,3)]
grad_l=np.array([[-1.,-1.,-1.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
from scipy.special import roots_jacobi
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
assert moment_count==56
E=100000.;nu=.3;mu=E/(2.*(1.+nu));lam=E*nu/((1.+nu)*(1.-2.*nu))
elastic=np.zeros((6,6));elastic[:3,:3]=lam;np.fill_diagonal(elastic[:3,:3],lam+2.*mu)
elastic[3:,3:]=np.eye(3)*mu
def tensor(e):
    return np.array([[e[0],e[3]/2.,e[5]/2.],[e[3]/2.,e[1],e[4]/2.],[e[5]/2.,e[4]/2.,e[2]]])
rows=[]
for order in [4,10]:
    mesh=copy.deepcopy(linear);mids={};cells=[]
    for tet in mesh['tetrahedra']:
        cell=tet.copy()
        if order==10:
            for i,j in edges:
                edge=tuple(sorted((tet[i],tet[j])))
                if edge not in mids:
                    mids[edge]=len(mesh['nodes'])
                    mesh['nodes'].append(((np.array(mesh['nodes'][edge[0]])+mesh['nodes'][edge[1]])/2.).tolist())
                cell.append(mids[edge])
        cells.append(cell)
    xyz=np.asarray(mesh['nodes']);dimension=3*len(xyz)
    K=np.zeros((dimension,dimension));points=[]
    for element,cell in enumerate(cells):
        ids=np.array([[3*n+a for a in range(3)] for n in cell]).reshape(-1)
        for rst,weight in ([(np.full(3,.25),1./6.)] if order==4 else quadrature):
            l=np.r_[1.-sum(rst),rst]
            if order==4:shape=l;derivative=grad_l
            else:
                shape=np.r_[l*(2.*l-1.),[4.*l[i]*l[j] for i,j in edges]]
                derivative=np.vstack([(4.*l[:,None]-1.)*grad_l,
                    [4.*(l[j]*grad_l[i]+l[i]*grad_l[j]) for i,j in edges]])
            J=xyz[cell].T@derivative;det=float(np.linalg.det(J));assert det>0.
            grad=derivative@np.linalg.inv(J);B=np.zeros((6,3*len(cell)))
            for i,(dx,dy,dz) in enumerate(grad):
                B[0,3*i]=dx;B[1,3*i+1]=dy;B[2,3*i+2]=dz
                B[3,3*i]=dy;B[3,3*i+1]=dx
                B[4,3*i+1]=dz;B[4,3*i+2]=dy
                B[5,3*i]=dz;B[5,3*i+2]=dx
            w=weight*det;K[np.ix_(ids,ids)]+=w*(B.T@elastic@B)
            points.append((element,ids,B,w,(shape@xyz[cell]).tolist()))
    assert abs(math.fsum(q[3] for q in points)-1.)<1e-14
    fixed=np.zeros(dimension,dtype=bool);fixed[1::3]=True
    fixed[np.flatnonzero(xyz[:,0]==0.)*3]=True
    fixed[np.flatnonzero(xyz[:,2]==0.)*3+2]=True
    free=np.flatnonzero(~fixed)
    owners={}
    for tet in mesh['tetrahedra']:
        for face in itertools.combinations(tet,3):
            key=tuple(sorted(face));owners[key]=owners.get(key,0)+1
    for shear,cohesion,label in [(0.,10.,'uniaxial'),(3.,100.,'shear-convention')]:
        pressure=30.;force=np.zeros(dimension)
        for face,count in owners.items():
            if count!=1 or not all(xyz[n,2]==1. for n in face):continue
            a,b,c=xyz[list(face)];area=float(np.linalg.norm(np.cross(b-a,c-a)))/2.
            loaded=list(face) if order==4 else [mids[tuple(sorted(pair))] for pair in itertools.combinations(face,2)]
            for n in loaded:force[3*n:3*n+3]+=area/3.*np.array([shear,0.,-pressure])
        assert abs(force[2::3].sum()+pressure)<1e-12
        u=np.zeros(dimension);u[free]=np.linalg.solve(K[np.ix_(free,free)],force[free])
        residual=float(np.linalg.norm((K@u-force)[free])/max(1.,np.linalg.norm(force)))
        assert residual<1e-12
        if shear==0.:
            exact=np.zeros_like(xyz);exact[:,0]=nu*(1.+nu)*pressure/E*xyz[:,0]
            exact[:,2]=-(1.-nu*nu)*pressure/E*xyz[:,2]
            assert np.max(np.abs(u-exact.reshape(-1)))<1e-14
        data=[];states=[]
        for element,ids,B,w,position in points:
            tension=elastic@(B@u[ids]);states.append(tension)
            if shear==0.:assert np.max(np.abs(tension-np.array([0.,-nu*pressure,-pressure,0.,0.,0.])))<1e-9
            sigma=-tensor(np.r_[tension[:3],2.*tension[3:]])
            eigen=np.linalg.eigvalsh(sigma);phi=math.radians(30.)
            f=(1.-math.sin(phi))*eigen[-1]-(1.+math.sin(phi))*eigen[0]-2.*cohesion*math.cos(phi)
            assert f<0.,'the independent subcritical elastic state must be MC-admissible'
            data.append(dict(element=element,position_m=position,weight_m3=w,tension_stress_kpa=tension.tolist(),original_MC_margin_kpa=float(f)))
        for factor in ([1.,1.5] if shear==0. else [1.]):
            phi=math.atan(math.tan(math.radians(30.))/factor);c=cohesion/factor
            ratio=(1.+math.sin(phi))/(1.-math.sin(phi));a=c/math.tan(phi)
            v=np.zeros_like(xyz);v[:,0]=(ratio+.2)*xyz[:,0]+(.02 if shear else 0.)*xyz[:,0]*xyz[:,2]
            v[:,2]=-xyz[:,2];flat=v.reshape(-1)
            assert np.all(flat[fixed]==0.)
            capacities=[];internal=[];wrong=[];dual=[]
            for q,stress,entry in zip(points,states,data):
                element,ids,B,w,_=q
                eng=B@flat[ids];D=tensor(eng);e=np.linalg.eigvalsh(D)
                margin=float(np.maximum(e,0.).sum()+ratio*np.minimum(e,0.).sum())
                assert margin>1e-5
                capacities.append(w*a*float(np.trace(D)));internal.append(w*float(stress@eng));dual.append(margin)
                wrong.append(w*float(stress[:3]@eng[:3]+stress[3:]@(eng[3:]/2.)))
                entry['virtual_engineering_strain_factor'+str(factor)]=eng.tolist()
            work=float(force@flat);capacity=math.fsum(capacities);actual=math.fsum(internal)
            assert abs(actual-work)<1e-9
            if shear:
                assert abs(math.fsum(wrong)-work)>1e-5
            else:
                assert (work>capacity)==(factor==1.5)
            rows.append(dict(order=order,case=label,factor=factor,nodes=len(xyz),elements=len(cells),GP=len(points),
                free_dofs=len(free),force_residual=residual,external_work=work,integrated_capacity=capacity,
                actual_elastic_internal_virtual_work=actual,work_error=abs(actual-work),
                minimum_dual_margin=min(dual),excluded=work>capacity,
                wrong_engineering_shear_work_error=abs(math.fsum(wrong)-work),
                source_elastic_stress_field_valid_at_factor1=True))
        dump(f'full-fields-Tet{order}-{label}.json',dict(nodes=xyz.tolist(),cells=cells,
            force_kn=force.tolist(),displacement_m=u.reshape(-1,3).tolist(),points=data))
assert len(rows)==6
dump('terminal.json',dict(passed=True,actual_cases=6,orders=[4,10],distributed_supports=True,
    quadrature_family='independent scipy.special shifted Gauss-Jacobi alpha2,1,0',
    all56_monomials_through_total_degree5_matched=True,original_Legendre_reference_v1_kept=True,
    independent_dense_elastic_FEM_and_analytical_displacements_verified=True,
    every_original_GP_MC_admissible=True,every_velocity_supported_component_exact_zero=True,
    every_GP_dual_margin_strictly_positive=True,
    all_actual_internal_and_external_virtual_work_matched=True,
    weak_factor1point5_excluded_and_original_factor1_never_excluded=True,
    shear_factor2_negative_controls_detected_for_both_orders=True,
    independent_UCS_at_F1=20.*math.sqrt(3.),
    independent_critical_F=2.*math.sqrt(10.*(30.*math.tan(math.radians(30.))+10.))/30.,
    no_Rust_search_or_nonlinear_FEM_or_actual_slope_acceptance_claimed=True,
    primary_MC_source='https://docs.itascacg.com/itasca950/common/models/mohr/doc/modelmohr.html',
    our_support_formula_derived_and_separately_LP_verified=True,
    driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    mesh_source_sha256=hashlib.sha256(base.read_bytes()).hexdigest(),p6_complete=False,rows=rows))
print(json.dumps(dict(passed=True,cases=6,Tet4_and_Tet10=True,
    distributed_support_virtual_work_verified=True,not_actual_slope_or_Rust_acceptance=True)),flush=True)
