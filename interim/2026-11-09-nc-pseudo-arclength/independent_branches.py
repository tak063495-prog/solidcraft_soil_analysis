"""Independent strain-path integration and local same-load roots, no native data.

The local ODE reuses the existing continuum MCC equations and existing RK4
event/drift checks. No stress is projected onto yield. Material count 4 or 8
is a numerical quadrature, not a material parameter or acceptance threshold.
"""
import argparse,json,sys
from pathlib import Path
import numpy as np
from independent_arc import read,write,ClosedStickSystem,physical,expected_inputs
import mcc_rate_fem as rate
from mcc_rate_rk4 import advance_rk4

class StrainRate:
    constitutive_terms=staticmethod(rate.terms)
    select=rate.System.select
    def __init__(self,parameters):self.parameters=parameters;self.calls=0;self.linear_iterations=0;self.dense_fallbacks=0
    def rate(self,state,strain,plastic):
        yf,d,n,dn,den=rate.terms(state,self.parameters)
        if np.any(den[plastic]<=0):raise ValueError('independent strain consistency denominator')
        potential=np.einsum('qa,qa->q',dn,strain)/np.where(plastic,den,np.maximum(abs(den),1.))
        flow=(potential*plastic)[:,None]*n;dy=np.zeros_like(state)
        dy[:,:6]=np.einsum('qab,qb->qa',d,strain-flow)
        dy[:,6]=np.exp(state[:,7])/(self.parameters[:,1]-self.parameters[:,2])*flow[:,:3].sum(axis=1)
        dy[:,7]=-strain[:,:3].sum(axis=1);dy[:,8:]=flow;self.calls+=1
        return np.zeros(1),dy,potential*np.linalg.norm(n,axis=1),yf
    def integrate(self,state,strain,count):return advance_rk4(self,state,np.zeros(1),strain,count)[0]

def run(reference,arc,out,count):
    out.mkdir(parents=True,exist_ok=False);mesh,ground,case=expected_inputs();samples=read(arc/'samples.json')
    peak=max(range(len(samples)),key=lambda i:samples[i]['progress']);assert peak>=8
    anchor=samples[peak-8];active=np.array(read(reference/'reference.json')['last_accepted']['active_elements'])
    oldfull=np.array(anchor['material_state']);old=oldfull[active].reshape(-1,14);nd=len(anchor['tangent'])-1
    system=ClosedStickSystem(mesh,ground,active,case);material=StrainRate(system.parameters)
    R=np.zeros((system.xyz.size,nd))
    for j,g in enumerate(system.free):
        v=np.zeros(system.xyz.size);v[g]=1.;R[:,j]=system.lift(v)
    B=np.einsum('eqai,eij->eqaj',system.b,R[system.dofs]).reshape(-1,6,nd);weights=system.weight.ravel()
    unload=read(reference/'unload.json');before=system.total_force(np.array(unload['material_state'])[active].reshape(-1,14),np.array(unload['displacement_m']).ravel())
    final=physical(mesh,ground,active,0.,True,joint_system=system.joint)[0].reshape(-1);target=anchor['progress'];load=before+target*(final-before)
    direction=np.array(anchor['tangent'])[:-1];direction/=np.linalg.norm(direction);guesses=[-4e-4,-2e-4,0.,1e-4,2e-4,3e-4,4e-4]
    write(out/'controls.json',dict(input_only=True,native_fields_used=False,native_fractions_used=False,anchor_sample=peak-8,peak_sample=peak,
        target_progress=target,material_rk4_steps=count,material_fd_width=1e-8,maximum_total_strain_l1=.001,
        force_tolerance=1e-8,true_linear_tolerance=1e-10,root_correction_check_m=1e-11,production_adopted=False))
    calls=0;records=[];solutions=[]
    def evaluate(du,tangent=True):
        nonlocal calls
        de=-np.einsum('qai,i->qa',B,du)
        if np.abs(de).sum(axis=1).max()>.001:raise ValueError('independent root total strain bound')
        state=material.integrate(old,de,count);calls+=1
        f=-np.einsum('qai,qa,q->i',B,state[:,:6],weights)
        if not tangent:return f,None,state
        D=np.zeros((len(old),6,6));h=1e-8
        for j in range(6):
            plus=de.copy();plus[:,j]+=h;minus=de.copy();minus[:,j]-=h
            D[:,:,j]=(material.integrate(old,plus,count)[:,:6]-material.integrate(old,minus,count)[:,:6])/(2*h);calls+=2
        K=np.einsum('qai,qab,qbj,q->ij',B,D,B,weights,optimize=True)
        return f,K,state
    for amplitude in guesses:
        du=amplitude*direction;steps=[];stop=None;converged=False
        try:
            for it in range(21):
                f,K,new=evaluate(du);rhs=R.T@load-f;rel=float(np.linalg.norm(rhs)/max(np.linalg.norm(load),1.));dx=np.linalg.solve(K,rhs)
                lr=float(np.linalg.norm(K@dx-rhs)/max(np.linalg.norm(rhs),1e-30));steps.append(dict(iteration=it,relative_force_residual=rel,correction_norm_m=float(np.linalg.norm(dx)),true_relative_linear_residual=lr))
                if lr>1e-10:raise ValueError('independent root linear gate')
                if it>=1 and rel<1e-8 and np.linalg.norm(dx)<1e-11:
                    full=oldfull.copy();full[active]=new.reshape(full[active].shape);stress=(old[:,:6]+new[:,:6])/2;de=-np.einsum('qai,i->qa',B,du);dep=new[:,8:]-old[:,8:]
                    pw=np.einsum('qa,qa->q',stress,dep);iw=float(np.sum(np.einsum('qa,qa->q',stress,de)*weights));ew=float(load@(R@du))
                    eig=np.linalg.eigvalsh((K+K.T)/2)
                    solutions.append(dict(progress=target,initial_amplitude_m=amplitude,reduced_increment_m=du.tolist(),
                        displacement_m=(np.array(anchor['displacement_m'])+(R@du).reshape(-1,3)).tolist(),material_state=full.tolist(),
                        relative_force_residual=rel,audit=dict(plastic_stress_work_kn_m=float(np.sum(pw*weights)),minimum_point_plastic_stress_work_kpa=float(pw.min()),
                        internal_work_kn_m=iw,external_work_kn_m=ew,work_difference_kn_m=iw-ew,symmetric_eigen_min_kn_m=float(eig[0]),
                        negative_symmetric_eigenvalues=int(np.sum(eig<0)),directional_incremental_work_kn_m=float(du@K@du),
                        removed_material_frozen=bool(np.array_equal(full[~active],oldfull[~active])),active_contact_points=0)))
                    converged=True;break
                if it==20:raise ValueError('independent root Newton limit')
                good=False
                for alpha in [1.,.5,.25,.125,.0625,.03125,.015625]:
                    candidate=du+alpha*dx
                    if np.abs(np.einsum('qai,i->qa',B,candidate)).sum(axis=1).max()>.001:continue
                    ff,_,_=evaluate(candidate,False)
                    if np.linalg.norm(R.T@load-ff)<np.linalg.norm(rhs):du=candidate;good=True;break
                if not good:raise ValueError('independent root line search exhausted')
        except (ValueError,np.linalg.LinAlgError) as exc:stop=str(exc)
        records.append(dict(initial_amplitude_m=amplitude,converged=converged,stop_reason=stop,iterations=steps));write(out/'trials.json',records);write(out/'solutions.json',solutions)
    report=dict(passed=False,history_converged=False,production_adopted=False,input_only=True,trials=len(records),converged_trials=sum(t['converged'] for t in records),
        material_map_evaluations=calls,material_rate_evaluations=material.calls,target_progress=target,
        limitations=['bounded local same-history roots only; not a disconnected/global branch search','no hardening free energy, dynamics or physical stability certificate'])
    write(out/'report.json',report);return report
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--reference',type=Path,required=True);p.add_argument('--arc',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--material-steps',type=int,choices=[4,8],default=4);a=p.parse_args()
    print(json.dumps(run(a.reference.resolve(),a.arc.resolve(),a.output.resolve(),a.material_steps),indent=2))
