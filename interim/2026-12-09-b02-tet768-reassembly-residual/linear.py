"""Linear-system observation and a discarded refinement experiment, no gate change."""
import numpy as np
from scipy.linalg import lu_solve
from precision import LD, norm

def diagnose(k, rhs, answer, factor):
    r64=k@answer-rhs
    rld=k.astype(LD)@answer.astype(LD)-rhs.astype(LD)
    correction=lu_solve(factor,np.asarray(-rld,float))
    trial=answer+correction
    refined=k.astype(LD)@trial.astype(LD)-rhs.astype(LD)
    scale=max(norm(rhs),1.)
    return dict(true_f64_relative=norm(r64)/scale,extended_relative=norm(rld)/scale,
                backward_error_inf=float(np.max(np.abs(rld))/(np.linalg.norm(k,np.inf)*np.max(np.abs(answer))+np.max(np.abs(rhs)))) if rhs.size else 0.,
                discarded_one_refinement_relative=norm(refined)/scale,discarded_correction_maximum=float(np.max(np.abs(correction),initial=0.)),
                rhs_norm=norm(rhs),solution_norm=norm(answer),matrix_inf_norm=float(np.linalg.norm(k,np.inf)),
                solution_updated=False,mantissa_bits=int(np.finfo(LD).nmant+1))
