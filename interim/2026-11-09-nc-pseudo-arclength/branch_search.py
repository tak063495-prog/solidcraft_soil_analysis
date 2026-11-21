"""Bounded fixed-history multi-start roots; no proof of global existence/absence."""
import argparse,json
from pathlib import Path
import numpy as np
from native_arc import Native,read,write

def run(evidence,binary,out):
    m=Native(evidence,binary,out);old=m.cp['material_states'];direction=np.linalg.solve(m.k0,m.L);direction/=np.linalg.norm(direction)
    targets=[m.lam0,0.9514002799987793,1.];guesses=[-4e-4,-2e-4,0.,1e-4,2e-4,3e-4,4e-4]
    trials=[];solutions=[]
    for target in targets:
        for guess in guesses:
            du=guess*direction;failure=None;steps=[];converged=False
            try:
                for it in range(21):
                    f,k,new=m.evaluate(du,old,f'lambda-{target}:guess-{guess}:iteration-{it}')
                    rhs=m.R.T@m.load(target)-f;rel=float(np.linalg.norm(rhs)/max(np.linalg.norm(m.load(target)),1.))
                    dx=np.linalg.solve(k,rhs);lr=float(np.linalg.norm(k@dx-rhs)/max(np.linalg.norm(rhs),1e-30))
                    steps.append(dict(iteration=it,relative_force_residual=rel,correction_norm_m=float(np.linalg.norm(dx)),true_linear_residual=lr))
                    if lr>1e-10:raise ValueError('fixed-history root true linear gate')
                    if it>=1 and rel<1e-8 and np.linalg.norm(dx)<1e-11:
                        row=dict(target_progress=target,initial_amplitude_m=guess,reduced_increment_m=du.tolist(),
                            displacement_m=(np.array(m.cp['displacement_m'])+(m.R@du).reshape(-1,3)).tolist(),states=new,
                            relative_force_residual=rel,audit=m.audit(old,new,du,m.lam0,target,k))
                        solutions.append(row);converged=True;break
                    if it==20:raise ValueError('bounded root Newton limit')
                    accepted=False
                    for alpha in [1.,.5,.25,.125,.0625,.03125,.015625]:
                        candidate=du+alpha*dx
                        if np.abs(np.einsum('qai,i->qa',m.rb,candidate)).sum(axis=1).max()>.001:continue
                        ff,_,_=m.evaluate(candidate,old,'root-line-search')
                        if np.linalg.norm(m.R.T@m.load(target)-ff)<np.linalg.norm(rhs):du=candidate;accepted=True;break
                    if not accepted:raise ValueError('bounded root line search no decreasing admissible trial')
            except (ValueError,np.linalg.LinAlgError) as exc:failure=str(exc)
            trials.append(dict(target_progress=target,initial_amplitude_m=guess,converged=converged,stop_reason=failure,iterations=steps))
            write(out/'trials.json',trials);write(out/'solutions.json',solutions)
    write(out/'processes.json',m.processes)
    report=dict(passed=False,production_adopted=False,targets=targets,initial_amplitudes_m=guesses,
        trials=len(trials),converged_trials=sum(t['converged'] for t in trials),native_committed_increments=0,
        limitations=['same native starting material history, bounded total strain <=0.001',
        'a second incremental root is not a disconnected branch proof or an admissible monotone jump',
        'no found root is not global nonexistence; no physical energy or stability acceptance'])
    write(out/'report.json',report);return report
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--evidence',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    print(json.dumps(run(a.evidence.resolve(),a.binary.resolve(),a.output.resolve()),indent=2))
