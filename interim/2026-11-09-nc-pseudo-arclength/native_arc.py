"""Bounded pseudo-arclength on the unchanged native material map, diagnostic only.

The evolving tangent hyperplane uses (u, a*lambda), a=1 m for numerical units.
Every trial starts from the last diagnostic state. No production state is committed.
"""
import argparse,hashlib,json,subprocess,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
from verify_high_order import prepare
from probe_nc_midpoint import inputs

def read(p):return json.loads(p.read_text())
def write(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')

class Native:
    def __init__(self,evidence,binary,out):
        self.out=out;out.mkdir(parents=True,exist_ok=False);self.binary=binary;self.call=0
        inp=evidence/'prior/fresh-head-midpoint/input'
        self.mesh,self.ground,self.case=[read(inp/(n+'.json')) for n in ('mesh','ground','stages')]
        assert (self.mesh,self.ground,self.case)==inputs('midpoint')
        self.q,self.offsets=prepare(self.mesh)
        self.seed=next(r['event'] for r in read(evidence/'field-records.json') if r['label']=='terminal-rejected:0')
        self.cp=self.seed['committed'];self.active=self.cp['active_elements']
        assert not any(p['active'] for p in self.seed['joint_report']['points']) and not np.any(self.seed['joint_internal'])
        nd=len(self.seed['reduced_residual']);self.R=np.zeros((len(self.mesh['nodes'])*3,nd))
        for g,row in enumerate(self.seed['dof_rows']):
            for i,c in row:self.R[g,i]=c
        self.rb=[];self.w=[];self.pi=[];self.par=[]
        params={m['id']:{**m['mcc'],'integration':self.case['settings']['mcc_integration']} for m in self.ground['materials']}
        for ei in np.flatnonzero(self.active):
            q=self.q[ei];gd=(q['nodes'][:,None]*3+np.arange(3)).ravel()
            for qi,w in enumerate(q['weight']):
                self.rb.append(q['b'][qi]@self.R[gd]);self.w.append(w);self.pi.append(int(self.offsets[ei]+qi));self.par.append(params[self.mesh['material_ids'][ei]])
        self.rb=np.array(self.rb);self.w=np.array(self.w);self.pi=np.array(self.pi)
        self.checkpoint=read(evidence/'prior/fresh-head-midpoint/job/recovery-checkpoint-0.json')
        assert self.checkpoint['run']['last_committed']==self.cp
        self.holding=np.array(self.checkpoint['run']['excavation_audits'][0]['initial_holding_load_kn']).ravel()
        self.load0=np.array(self.cp['applied_nodal_load_kn']).ravel();self.lam0=self.cp['stage_progress'];self.L=-self.R.T@self.holding
        self.k0=np.array(read(evidence/'assembled-matrices.json')['last-accepted:7:ad'])
        self.processes=[]
    def load(self,lam):return self.load0-(lam-self.lam0)*self.holding
    def evaluate(self,du,old,label):
        de=-np.einsum('qai,i->qa',self.rb,du)
        if np.abs(de).sum(axis=1).max()>.001:raise ValueError('native fixed partition strain bound')
        req=[dict(label=label,point=int(pi),parameters=p,state=old[pi],strain=e.tolist(),partitions=1,widths=[]) for pi,p,e in zip(self.pi,self.par,de)]
        ip=self.out/f'call-{self.call:05d}-request.json';op=self.out/f'call-{self.call:05d}-result.json';write(ip,req)
        t=time.monotonic();args=[str(self.binary),str(ip),str(op)]
        with (self.out/f'call-{self.call:05d}.log').open('x') as log:r=subprocess.run(args,stdout=log,stderr=subprocess.STDOUT,timeout=30)
        self.processes.append(dict(command=args,exit_code=r.returncode,elapsed_seconds=time.monotonic()-t));r.check_returncode()
        points=read(op)['points'];self.call+=1
        if any('error' in p for p in points):raise ValueError('material error '+str([p for p in points if 'error' in p][:1]))
        assert all(p['ad_primal_exact'] and p['substep_replay_exact'] for p in points)
        states=list(old)
        for pi,p in zip(self.pi,points):states[pi]=p['updated_state']
        stress=np.array([states[pi]['effective_stress_compression_positive_kpa'] for pi in self.pi])
        k=np.einsum('qai,qab,qbj,q->ij',self.rb,np.array([p['implicit_tangent'] for p in points]),self.rb,self.w,optimize=True)
        force=-np.einsum('qai,qa,q->i',self.rb,stress,self.w)
        return force,k,states
    def audit(self,old,new,du,lam_before,lam,k):
        def field(s,key):return np.array([s[pi][key] for pi in self.pi])
        stress=(field(old,'effective_stress_compression_positive_kpa')+field(new,'effective_stress_compression_positive_kpa'))/2
        de=-np.einsum('qai,i->qa',self.rb,du);dep=field(new,'plastic_strain_compression_positive')-field(old,'plastic_strain_compression_positive')
        pw=np.einsum('qa,qa->q',stress,dep)
        total=float(np.sum(np.einsum('qa,qa->q',stress,de)*self.w))
        plastic=float(np.sum(pw*self.w));external=float(((self.load(lam_before)+self.load(lam))/2)@(self.R@du))
        eig=np.linalg.eigvalsh((k+k.T)/2)
        return dict(internal_work_kn_m=total,external_work_kn_m=external,work_difference_kn_m=total-external,
                    plastic_stress_work_kn_m=plastic,minimum_point_plastic_stress_work_kpa=float(pw.min()),
                    elastic_stress_work_kn_m=total-plastic,symmetric_eigen_min_kn_m=float(eig[0]),
                    negative_symmetric_eigenvalues=int(np.sum(eig<0)),directional_incremental_work_kn_m=float(du@k@du),
                    removed_material_frozen=all(new[i]==self.cp['material_states'][i] for i in range(len(new)) if i not in self.pi),
                    active_contact_points=0,contact_history_frozen=True,
                    energy_limit='stress work split only; no specified hardening free energy or physical stability verdict')

def tangent(k,L,previous,a=1.):
    n=len(L);A=np.zeros((n+1,n+1));A[:n,:n]=k;A[:n,n]=-L;A[n,:n]=previous[:n];A[n,n]=a*a*previous[n]
    rhs=np.r_[np.zeros(n),1.];x=np.linalg.solve(A,rhs)
    linear=np.linalg.norm(A@x-rhs)/max(np.linalg.norm(rhs),1e-30)
    if linear>1e-10:raise ValueError('tangent augmented linear residual')
    x/=np.sqrt(x[:-1]@x[:-1]+(a*x[-1])**2)
    if previous[:-1]@x[:-1]+a*a*previous[-1]*x[-1]<0:x=-x
    return x

def run(evidence,binary,out,ds,count,resume=None):
    m=Native(evidence,binary,out);n=len(m.L);a=1.;events=[];samples=[];linear_records=[]
    if resume:
        saved=read(resume);old=saved['states'];u=np.array(saved['displacement_m']);lam=saved['progress'];t=np.array(saved['tangent']);start_index=saved['step']+1
        assert saved['diagnostic_only'] and saved['step_size_m']==ds and saved['seed_checkpoint_sha256']==hashlib.sha256(json.dumps(m.cp,sort_keys=True).encode()).hexdigest()
    else:
        old=m.cp['material_states'];u=np.array(m.cp['displacement_m']);lam=m.lam0;rate=np.linalg.solve(m.k0,m.L)
        t=np.r_[rate,1.];t/=np.linalg.norm(t);start_index=0
    write(out/'controls.json',dict(step_size_m=ds,steps=count,arc_lambda_scale_m=a,force_tolerance=1e-8,true_linear_tolerance=1e-10,
        arc_constraint_tolerance_m=1e-12,maximum_newton_iterations=20,minimum_correctors=1,native_committed_increments=0,
        input_only=False,resume_from=str(resume) if resume else None,production_adopted=False,passed=False))
    f,k,s=m.evaluate(np.zeros(n),old,'seed-zero');seedres=float(np.linalg.norm(m.R.T@m.load(lam)-f)/max(np.linalg.norm(m.load(lam)),1.))
    write(out/'seed-check.json',dict(relative_force_residual=seedres,state_exact=s==old))
    stop=None
    try:
        for index in range(start_index,start_index+count):
            du=ds*t[:-1];next_lam=lam+ds*t[-1]
            for it in range(21):
                f,k,new=m.evaluate(du,old,f'step-{index}:newton-{it}');target=m.load(next_lam);rhs=m.R.T@target-f
                rel=float(np.linalg.norm(rhs)/max(np.linalg.norm(target),1.));g=ds-t[:-1]@du-a*a*t[-1]*(next_lam-lam)
                events.append(dict(step=index,iteration=it,call=m.call-1,progress=float(next_lam),relative_force_residual=rel,arc_residual_m=float(g)))
                if it>=1 and rel<1e-8 and abs(g)<1e-12:
                    next_t=tangent(k,m.L,t,a);u+=(m.R@du).reshape(u.shape)
                    row=dict(step=index,progress=float(next_lam),previous_progress=float(lam),progress_increment=float(next_lam-lam),
                        displacement_m=u.tolist(),states=new,tangent=next_t.tolist(),relative_force_residual=rel,arc_residual_m=float(g),
                        audit=m.audit(old,new,du,lam,next_lam,k),diagnostic_only=True,step_size_m=ds,
                        seed_checkpoint_sha256=hashlib.sha256(json.dumps(m.cp,sort_keys=True).encode()).hexdigest())
                    samples.append(row);old=new;lam=next_lam;t=next_t
                    write(out/'samples.json',samples);write(out/'restart.json',row);break
                if it==20:raise ValueError('native pseudo-arclength Newton limit')
                A=np.zeros((n+1,n+1));A[:n,:n]=k;A[:n,n]=-m.L;A[n,:n]=t[:-1]*1000;A[n,n]=a*a*t[-1]*1000
                b=np.r_[rhs,g*1000];x=np.linalg.solve(A,b);lr=float(np.linalg.norm(A@x-b)/max(np.linalg.norm(b),1e-30))
                linear_records.append(dict(step=index,iteration=it,matrix=A.tolist(),rhs=b.tolist(),solution=x.tolist(),true_relative_residual=lr))
                if lr>1e-10:raise ValueError('native pseudo-arclength true linear residual')
                du+=x[:-1];next_lam+=x[-1]
    except (ValueError,np.linalg.LinAlgError,subprocess.SubprocessError) as exc:stop=str(exc)
    write(out/'events.json',events);write(out/'linear-systems.json',linear_records);write(out/'processes.json',m.processes)
    report=dict(passed=False,history_converged=False,production_adopted=False,native_committed_increments=0,
        equilibrium_samples=len(samples),maximum_sampled_progress=max([m.lam0]+[r['progress'] for r in samples]),
        negative_progress_increments=sum(r['progress_increment']<0 for r in samples),positive_progress_increments=sum(r['progress_increment']>0 for r in samples),
        stop_reason=stop,seed_relative_force_residual=seedres,full_removal_reached=False,material_calls=m.call)
    write(out/'report.json',report);return report

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--evidence',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--step-size',type=float,required=True);p.add_argument('--count',type=int,required=True);p.add_argument('--resume',type=Path);a=p.parse_args()
    assert 0<a.step_size<=1e-4 and 0<a.count<=256
    print(json.dumps(run(a.evidence.resolve(),a.binary.resolve(),a.output.resolve(),a.step_size,a.count,a.resume.resolve() if a.resume else None),indent=2))
