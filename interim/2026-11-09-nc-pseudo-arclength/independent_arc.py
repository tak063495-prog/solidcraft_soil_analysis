"""Input-only independent continuum MCC arclength evolution, without projection.

Starting history is generated exclusively by low_preload_nc_reference, never by
native fractions or fields. Analytical continuum tangents and RK4 yield/unloading
events are the existing independent rate equations. The bordered tangent fixes
orientation and normalizes arc speed in (u, 1 m * lambda). This is a rate trace,
not a production equilibrium corrector or a thermodynamic stability certificate.
"""
import argparse,json,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'tools'))
import mcc_rate_fem as rate
from mcc_rate_rk4 import advance_rk4
from low_preload_nc_reference import ClosedStickSystem
from probe_low_preload_nc import expected_inputs
from verify_high_order import physical

def read(p):return json.loads(p.read_text())
def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')

class ArcSystem(ClosedStickSystem):
    def __init__(self,*args):
        super().__init__(*args);self.orientation=None;self.maximum_bordered_residual=0.;self.last_matrix=None;self.last_tangent=None;self.saved_systems=[]
    def rate(self,state,arc_length,plastic):
        yf,d,n,dn,den=rate.terms(state,self.parameters)
        if np.any(den[plastic]<=0):raise ValueError('independent nonpositive consistency denominator')
        ep=d.copy();ep[plastic]-=np.einsum('qa,qb,q->qab',dn[plastic],dn[plastic],1/den[plastic])
        ke=np.einsum('eqai,eqab,eqbj,eq->eij',self.rb,ep.reshape(*self.weight.shape,6,6),self.rb,self.weight,optimize=True)
        values=np.bincount(self.inverse,weights=ke.reshape(-1)[self.mask],minlength=len(self.sparse_rows))
        K=self.kr.copy();K[self.sparse_rows,self.sparse_cols]+=values
        nd=len(self.free);A=np.zeros((nd+1,nd+1));A[:nd,:nd]=K;A[:nd,nd]=-self.force_derivative
        A[nd,:]=self.orientation;rhs=np.r_[np.zeros(nd),1.];x=np.linalg.solve(A,rhs)
        linear=float(np.linalg.norm(A@x-rhs)/max(np.linalg.norm(rhs),1e-30));self.maximum_bordered_residual=max(self.maximum_bordered_residual,linear)
        if linear>1e-10:raise ValueError('independent augmented linear residual exceeds existing 1e-10')
        x/=np.linalg.norm(x)
        if self.orientation@x<0:x=-x
        coordinates=np.zeros(self.xyz.size);coordinates[self.free]=x[:-1];du=self.lift(coordinates)
        de=-np.einsum('eqai,ei->eqa',self.b,du[self.dofs]).reshape(-1,6)
        potential=np.einsum('qa,qa->q',dn,de)/np.where(plastic,den,np.maximum(abs(den),1.))
        flow=(potential*plastic)[:,None]*n
        dy=np.zeros_like(state);dy[:,:6]=np.einsum('qab,qb->qa',d,de-flow)
        dy[:,6]=np.exp(state[:,7])/(self.parameters[:,1]-self.parameters[:,2])*flow[:,:3].sum(axis=1)
        dy[:,7]=-de[:,:3].sum(axis=1);dy[:,8:]=flow
        self.last_matrix=K;self.last_tangent=x;self.calls+=1
        # Save actual solved systems, not just reported residuals, for auditing.
        self.saved_systems.append(dict(matrix=A.tolist(),rhs=rhs.tolist(),solution=(x/(self.orientation@x)).tolist(),true_relative_residual=linear))
        return np.r_[du,x[-1]]*arc_length,dy*arc_length,potential*np.linalg.norm(n,axis=1)*arc_length,yf

def run(reference,out,ds,count,resume=None):
    out.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    mesh,ground,case=expected_inputs()
    assert [read(reference/(name+'.json')) for name in ('mesh','ground','stages')]==[mesh,ground,case]
    source=read(reference/'reference.json');seed=source['last_accepted'];assert seed['name']=='excavate_upper'
    active=np.array(seed['active_elements']);full=np.array(seed['material_state']);u=np.array(seed['displacement_m']).ravel();lam=seed['progress']
    system=ArcSystem(mesh,ground,active,case);assert not system.live.any() and not np.any(system.kjoint)
    # Independently recover original residual-release force from the input-only
    # pre-excavation field; do not obtain it from the current endpoint stresses.
    unload=read(reference/'unload.json');before=system.total_force(np.array(unload['material_state'])[active].reshape(-1,14),np.array(unload['displacement_m']).ravel())
    final=physical(mesh,ground,active,0.,True,joint_system=system.joint)[0].reshape(-1)
    system.force_derivative=system.reduce(final-before)
    frozen=full[~active].copy();contacts=seed['contacts'];state=full[active].reshape(-1,14);z=np.r_[u,lam];frames=[];start_index=0
    # Orientation is selected independently from the continuum loading tangent.
    fy,d,n,dn,den=rate.terms(state,system.parameters)
    elastic=np.einsum('eqai,eqab,eqbj,eq->eij',system.rb,d.reshape(*system.weight.shape,6,6),system.rb,system.weight,optimize=True)
    v=np.bincount(system.inverse,weights=elastic.reshape(-1)[system.mask],minlength=len(system.sparse_rows));K=system.kr.copy();K[system.sparse_rows,system.sparse_cols]+=v
    x=np.r_[np.linalg.solve(K,system.force_derivative),1.];system.orientation=x/np.linalg.norm(x)
    if resume:
        s=read(resume);assert s['input_only'] and s['step_size_m']==ds
        state=np.array(s['material_state'])[active].reshape(-1,14);full=np.array(s['material_state']);z=np.r_[np.array(s['displacement_m']).ravel(),s['progress']]
        system.orientation=np.array(s['tangent']);start_index=s['step']+1
    zerores=float(np.linalg.norm(system.reduce(system.total_force(state,z[:-1])-(before+z[-1]*(final-before))))/max(np.linalg.norm(before+z[-1]*(final-before)),1.))
    write(out/'controls.json',dict(input_only=True,native_fields_used=False,native_progress_used=False,
        initial_reference_failure=source['failure_reason'],seed_progress=float(z[-1]),seed_relative_force_residual=zerores,
        step_size_m=ds,count=count,arc_lambda_scale_m=1.,method='normalized bordered continuum rate + existing RK4 events',
        existing_yield_drift_limit=2e-6,existing_true_linear_tolerance=1e-10,production_adopted=False,passed=False))
    previous_state=state.copy();previous_z=z.copy();stop=None;segments=[]
    try:
        for index in range(start_index,start_index+count):
            last_callback=None
            def observe(t,s,disp):
                nonlocal last_callback
                last_callback=(s.copy(),disp.copy())
            state,z,diag=advance_rk4(system,state,z,ds,2,accepted_callback=observe)
            segments.append(diag)
            p=system.select(state,ds,np.zeros(len(state),dtype=bool));system.rate(state,ds,p)
            tangent=system.last_tangent.copy();K=system.last_matrix;system.orientation=tangent
            full[active]=state.reshape(full[active].shape)
            target=before+z[-1]*(final-before);res=float(np.linalg.norm(system.reduce(system.total_force(state,z[:-1])-target))/max(np.linalg.norm(target),1.))
            de=-np.einsum('eqai,ei->eqa',system.b,(z[:-1]-previous_z[:-1])[system.dofs]).reshape(-1,6)
            stress=(previous_state[:,:6]+state[:,:6])/2;dep=state[:,8:]-previous_state[:,8:];pw=np.einsum('qa,qa->q',stress,dep)
            iw=float(np.sum(np.einsum('qa,qa->q',stress,de)*system.weight.ravel()));pwtotal=float(np.sum(pw*system.weight.ravel()))
            ew=float((before+(previous_z[-1]+z[-1])*(final-before)/2)@(z[:-1]-previous_z[:-1]))
            du=system.reduce(z[:-1]-previous_z[:-1]);eig=np.linalg.eigvalsh((K+K.T)/2)
            row=dict(step=index,progress=float(z[-1]),previous_progress=float(previous_z[-1]),progress_increment=float(z[-1]-previous_z[-1]),
                displacement_m=z[:-1].reshape(-1,3).tolist(),material_state=full.tolist(),tangent=tangent.tolist(),input_only=True,step_size_m=ds,
                relative_force_residual=res,maximum_normalized_yield_drift=float(np.max(np.maximum(rate.terms(state,system.parameters)[0],0.))),
                contacts=contacts,audit=dict(internal_work_kn_m=iw,external_work_kn_m=ew,work_difference_kn_m=iw-ew,
                plastic_stress_work_kn_m=pwtotal,elastic_stress_work_kn_m=iw-pwtotal,minimum_point_plastic_stress_work_kpa=float(pw.min()),
                symmetric_eigen_min_kn_m=float(eig[0]),negative_symmetric_eigenvalues=int(np.sum(eig<0)),
                directional_incremental_work_kn_m=float(du@K@du),removed_material_frozen=bool(np.array_equal(full[~active],frozen)),
                active_contact_points=0,contact_history_frozen=True))
            frames.append(row);write(out/'samples.json',frames);write(out/'restart.json',row)
            previous_state=state.copy();previous_z=z.copy()
            if time.monotonic()-started>180:raise ValueError('independent arc wall time budget')
    except (ValueError,np.linalg.LinAlgError) as exc:stop=str(exc)
    write(out/'linear-systems.json',system.saved_systems);write(out/'segments.json',segments)
    report=dict(passed=False,history_converged=False,production_adopted=False,input_only=True,equilibrium_samples=len(frames),
        maximum_sampled_progress=max([seed['progress']]+[f['progress'] for f in frames]),
        positive_progress_increments=sum(f['progress_increment']>0 for f in frames),negative_progress_increments=sum(f['progress_increment']<0 for f in frames),
        maximum_bordered_true_linear_residual=system.maximum_bordered_residual,stop_reason=stop,elapsed_seconds=time.monotonic()-started,
        full_removal_reached=False,limitations=['continuum rate arclength trace, not a native map or production Run',
        'no yield projection; stop on existing drift/event gates','plastic stress work does not specify hardening free energy or physical stability'])
    write(out/'report.json',report);return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--reference',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--step-size',type=float,required=True);p.add_argument('--count',type=int,required=True);p.add_argument('--resume',type=Path);a=p.parse_args()
    assert 0<a.step_size<=1e-4 and 0<a.count<=256
    print(json.dumps(run(a.reference.resolve(),a.output.resolve(),a.step_size,a.count,a.resume.resolve() if a.resume else None),indent=2))
