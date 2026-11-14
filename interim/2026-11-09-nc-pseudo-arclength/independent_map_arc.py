"""Pseudo-arclength correctors with an input-only strain-driven MCC RK4 map.

This second independent check uses immutable step-start histories, numerical
material derivatives and the actual force corrector, rather than only rate
integration. It never reads native fields, progress, material calls or matrices.
"""
import argparse,json,sys,time
from pathlib import Path
import numpy as np
from independent_branches import StrainRate,ClosedStickSystem,physical,expected_inputs,read,write

def read_arc_seed(arc):
    """A missing/empty rate trace is failure, never an equilibrium origin."""
    if arc is None or not (arc/'samples.json').is_file():
        raise ValueError('independent arc seed unavailable: samples.json missing: '+str(arc))
    try:
        samples=read(arc/'samples.json')
        report=read(arc/'report.json')
    except (OSError,ValueError) as exc:
        raise ValueError('independent arc seed unreadable: '+str(exc)) from exc
    if not isinstance(report,dict):raise ValueError('independent arc report.json is not an object')
    count=report.get('equilibrium_samples',0)
    if not isinstance(samples,list) or not samples or type(count) is not int or count<=0:
        raise ValueError('independent arc seed unavailable: equilibrium_samples=0 or empty samples.json')
    if len(samples)!=report['equilibrium_samples']:
        raise ValueError('independent arc sample count differs from report.json')
    if report.get('diagnostic_status')=='failed':
        raise ValueError('independent arc seed unavailable: failed trace; retain partial samples: '+str(report.get('stop_reason')))
    return samples[0]


def seed_failure(out,exc):
    report=dict(passed=False,nc01_complete=False,history_converged=False,production_adopted=False,
                equilibrium_samples=0,full_removal_reached=False,static_limit_confirmed=False,
                diagnostic_status='failed',stop_reason=str(exc),failure_stage='independent_arc_seed')
    write(out/'report.json',report)
    write(out/'failure.json',report)
    (out/'failure.log').write_text(str(exc)+'\n',encoding='utf-8')
    return report


class Model:
    def __init__(self,reference,arc,count):
        self.mesh,self.ground,self.case=expected_inputs();self.seed=read_arc_seed(arc)
        self.active=np.array(read(reference/'reference.json')['last_accepted']['active_elements']);self.system=ClosedStickSystem(self.mesh,self.ground,self.active,self.case)
        self.oldfull=np.array(self.seed['material_state']);self.material=StrainRate(self.system.parameters);self.count=count;self.calls=0
        s=self.system;n=len(s.free);self.R=np.zeros((s.xyz.size,n))
        for j,g in enumerate(s.free):v=np.zeros(s.xyz.size);v[g]=1.;self.R[:,j]=s.lift(v)
        self.B=np.einsum('eqai,eij->eqaj',s.b,self.R[s.dofs]).reshape(-1,6,n);self.w=s.weight.ravel()
        unload=read(reference/'unload.json');self.before=s.total_force(np.array(unload['material_state'])[self.active].reshape(-1,14),np.array(unload['displacement_m']).ravel())
        final=physical(self.mesh,self.ground,self.active,0.,True,joint_system=s.joint)[0].reshape(-1);self.delta=final-self.before;self.L=self.R.T@self.delta
    def load(self,lam):return self.before+lam*self.delta
    def evaluate(self,du,full):
        de=-np.einsum('qai,i->qa',self.B,du)
        if np.abs(de).sum(axis=1).max()>.001:raise ValueError('independent map arc strain bound')
        old=full[self.active].reshape(-1,14);new=self.material.integrate(old,de,self.count);self.calls+=1
        D=np.zeros((len(old),6,6));h=1e-8
        for j in range(6):
            plus=de.copy();plus[:,j]+=h;minus=de.copy();minus[:,j]-=h
            D[:,:,j]=(self.material.integrate(old,plus,self.count)[:,:6]-self.material.integrate(old,minus,self.count)[:,:6])/(2*h);self.calls+=2
        f=-np.einsum('qai,qa,q->i',self.B,new[:,:6],self.w);K=np.einsum('qai,qab,qbj,q->ij',self.B,D,self.B,self.w,optimize=True)
        result=full.copy();result[self.active]=new.reshape(full[self.active].shape);return f,K,result
    def audit(self,old,new,du,lam0,lam,K):
        x=old[self.active].reshape(-1,14);y=new[self.active].reshape(-1,14);stress=(x[:,:6]+y[:,:6])/2;de=-np.einsum('qai,i->qa',self.B,du);dep=y[:,8:]-x[:,8:]
        pw=np.einsum('qa,qa->q',stress,dep);iw=float(np.sum(np.einsum('qa,qa->q',stress,de)*self.w));ew=float((self.load(lam0)+self.load(lam))/2@(self.R@du));eig=np.linalg.eigvalsh((K+K.T)/2)
        return dict(internal_work_kn_m=iw,external_work_kn_m=ew,work_difference_kn_m=iw-ew,
            plastic_stress_work_kn_m=float(np.sum(pw*self.w)),minimum_point_plastic_stress_work_kpa=float(pw.min()),
            symmetric_eigen_min_kn_m=float(eig[0]),negative_symmetric_eigenvalues=int(np.sum(eig<0)),directional_incremental_work_kn_m=float(du@K@du),
            removed_material_frozen=bool(np.array_equal(new[~self.active],self.oldfull[~self.active])),active_contact_points=0,contact_history_frozen=True)

def run(reference,arc,out,ds,count,material_steps):
    out.mkdir(parents=True,exist_ok=False)
    try:read_arc_seed(arc)
    except ValueError as exc:return seed_failure(out,exc)
    m=Model(reference,arc,material_steps);old=m.oldfull.copy();u=np.array(m.seed['displacement_m']);lam=m.seed['progress'];t=np.array(m.seed['tangent']);n=len(m.L)
    events=[];samples=[];systems=[];stop=None;started=time.monotonic()
    write(out/'controls.json',dict(input_only=True,native_fields_or_fractions_used=False,source_reference=str(reference),source_arc=str(arc),seed_sample=0,
        step_size_m=ds,count=count,material_rk4_steps=material_steps,material_fd_width=1e-8,arc_lambda_scale_m=1.,force_tolerance=1e-8,true_linear_tolerance=1e-10,arc_constraint_tolerance_m=1e-12,maximum_newton_iterations=20,production_adopted=False))
    write(out/'seed.json',m.seed)
    try:
        for index in range(count):
            du=ds*t[:-1];next_lam=lam+ds*t[-1]
            for it in range(21):
                f,K,new=m.evaluate(du,old);target=m.load(next_lam);rhs=m.R.T@target-f;rel=float(np.linalg.norm(rhs)/max(np.linalg.norm(target),1.));g=ds-t[:-1]@du-t[-1]*(next_lam-lam)
                events.append(dict(step=index,iteration=it,progress=float(next_lam),reduced_increment_m=du.tolist(),relative_force_residual=rel,arc_residual_m=float(g)))
                A=np.zeros((n+1,n+1));A[:n,:n]=K;A[:n,n]=-m.L;A[n,:n]=t[:-1]*1000;A[n,n]=t[-1]*1000
                if it>=1 and rel<1e-8 and abs(g)<1e-12:
                    b=np.r_[np.zeros(n),1000.];next_t=np.linalg.solve(A,b);lr=float(np.linalg.norm(A@next_t-b)/np.linalg.norm(b))
                    systems.append(dict(kind='tangent',matrix=A.tolist(),rhs=b.tolist(),solution=next_t.tolist(),true_relative_residual=lr))
                    if lr>1e-10:raise ValueError('independent map arc tangent linear gate')
                    next_t/=np.linalg.norm(next_t)
                    if t@next_t<0:next_t=-next_t
                    u+=(m.R@du).reshape(-1,3)
                    samples.append(dict(step=index,progress=float(next_lam),progress_increment=float(next_lam-lam),previous_progress=float(lam),relative_force_residual=rel,arc_residual_m=float(g),
                        displacement_m=u.tolist(),material_state=new.tolist(),tangent=next_t.tolist(),audit=m.audit(old,new,du,lam,next_lam,K)))
                    old=new;lam=next_lam;t=next_t;write(out/'samples.json',samples);break
                if it==20:raise ValueError('independent map arc Newton limit')
                b=np.r_[rhs,g*1000];dx=np.linalg.solve(A,b);lr=float(np.linalg.norm(A@dx-b)/max(np.linalg.norm(b),1e-30));systems.append(dict(kind='corrector',matrix=A.tolist(),rhs=b.tolist(),solution=dx.tolist(),true_relative_residual=lr))
                if lr>1e-10:raise ValueError('independent map arc corrector linear gate')
                du+=dx[:-1];next_lam+=dx[-1]
            if time.monotonic()-started>180:raise ValueError('independent map arc wall time budget')
    except (ValueError,np.linalg.LinAlgError) as exc:stop=str(exc)
    write(out/'events.json',events);write(out/'linear-systems.json',systems)
    report=dict(passed=False,history_converged=False,production_adopted=False,input_only=True,equilibrium_samples=len(samples),maximum_sampled_progress=max([m.seed['progress']]+[s['progress'] for s in samples]),
        positive_progress_increments=sum(s['progress_increment']>0 for s in samples),negative_progress_increments=sum(s['progress_increment']<0 for s in samples),
        material_map_evaluations=m.calls,material_rate_evaluations=m.material.calls,stop_reason=stop,elapsed_seconds=time.monotonic()-started,full_removal_reached=False)
    write(out/'report.json',report);return report
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--reference',type=Path,required=True);p.add_argument('--arc',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--step-size',type=float,required=True);p.add_argument('--count',type=int,required=True);p.add_argument('--material-steps',type=int,choices=[4,8],default=4);a=p.parse_args()
    assert 0<a.step_size<=1e-4 and 0<a.count<=256
    result=run(a.reference.resolve(),a.arc.resolve(),a.output.resolve(),a.step_size,a.count,a.material_steps)
    print(json.dumps(result,indent=2))
    sys.exit(1 if result.get('diagnostic_status')=='failed' or result['equilibrium_samples']==0 else 0)
