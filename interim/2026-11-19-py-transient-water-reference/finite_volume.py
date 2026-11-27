"""Independent cell-centred conservative FD/FV in the physical head variable.

Richards uses nonlinear mixed-form storage and arithmetic face K, not U's
analytical linearized PDE. Backward Euler + damped Newton; no core imports.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_banded
from scipy.sparse import lil_matrix

from reference import ALPHA, KS, SS_TERZAGHI, boundary, primary, storage


def faces(name, h, t):
    n = len(h)
    dx = 1/n
    left, right = boundary(name, t)
    rich = name == "richards"
    gravity = 1. if rich else 0.
    a = ALPHA if rich else 0.
    conductivity = SS_TERZAGHI if name == "terzaghi" else 1.
    kk = KS*np.exp(a*h) if rich else np.full(n, conductivity)
    kb = KS*np.exp(a*left) if rich else conductivity
    q, dl, dr = np.zeros(n+1), np.zeros(n+1), np.zeros(n+1)
    k = (kk[:-1]+kk[1:])/2
    grad = np.diff(h)/dx+gravity
    q[1:n] = -k*grad
    dl[1:n] = k/dx-.5*a*kk[:-1]*grad
    dr[1:n] = -k/dx-.5*a*kk[1:]*grad
    k0 = (kb+kk[0])/2
    g0 = (h[0]-left)/(dx/2)+gravity
    q[0] = -k0*g0
    dr[0] = -k0/(dx/2)-.5*a*kk[0]*g0
    if right[0] == "flux":
        q[-1] = right[1]
    else:
        q[-1] = -conductivity*(right[1]-h[-1])/(dx/2)
        dl[-1] = conductivity/(dx/2)
    return q, dl, dr


def flows(q):
    return max(q[0], 0)+max(-q[-1], 0), max(-q[0], 0)+max(q[-1], 0)


def residual(name, h, old_storage, t, dt):
    s, cap = storage(name, h)
    q, dl, dr = faces(name, h, t)
    dx = 1/len(h)
    r = dx*(s-old_storage)/dt+np.diff(q)
    band = np.zeros((3, len(h)))
    band[0, 1:] = dr[1:-1]
    band[1] = dx*cap/dt+dl[1:]-dr[:-1]
    band[2, :-1] = -dl[1:-1]
    return r, band, q


def backward_euler(name, n, start, end, steps, initial=None):
    if n < 2 or steps < 1 or end <= start:
        raise ValueError("invalid grid or time interval")
    x = (np.arange(n)+.5)/n
    h = primary(name, x, start) if initial is None else np.array(initial, float)
    if h.shape != x.shape or not np.all(np.isfinite(h)):
        raise ValueError("invalid initial state")
    dt = (end-start)/steps
    s0 = storage(name, h)[0].sum()/n
    inflow = outflow = max_step_defect = 0.
    total_iterations = 0
    for step in range(1, steps+1):
        t = start+step*dt
        old = storage(name, h)[0]
        trial = h.copy()
        for iteration in range(30):
            r, band, q = residual(name, trial, old, t, dt)
            norm = np.max(np.abs(r))
            if norm < 2e-13:
                break
            delta = solve_banded((1, 1), band, -r)
            factor = 1.
            for _ in range(24):
                candidate = trial+factor*delta
                if name != "richards" or np.all(candidate <= 1e-12):
                    rr = residual(name, candidate, old, t, dt)[0]
                    if np.all(np.isfinite(rr)) and np.max(np.abs(rr)) < norm:
                        trial = candidate
                        break
                factor /= 2
            else:
                raise RuntimeError("Newton line search failed")
        else:
            raise RuntimeError("Newton iteration limit; not accepted")
        total_iterations += iteration
        h = trial
        fi, fo = flows(q)
        inflow += dt*fi
        outflow += dt*fo
        change = (storage(name, h)[0]-old).sum()/n
        max_step_defect = max(max_step_defect, abs(change-dt*(q[0]-q[-1])))
    delta_storage = storage(name, h)[0].sum()/n-s0
    defect = delta_storage-(inflow-outflow)
    scale = max(abs(delta_storage), inflow+outflow, 1e-30)
    return h, {"n": n, "steps": steps, "dx": 1/n, "dt": dt,
               "start": start, "end": end, "inflow": inflow, "outflow": outflow,
               "delta_storage": delta_storage, "mass_defect": defect,
               "mass_relative": abs(defect)/scale, "max_step_defect": max_step_defect,
               "newton_iterations": total_iterations, "completed": True}


def semidiscrete(name, n, start, end, rtol=2e-10, atol=2e-12):
    """High accuracy method-of-lines time control on exactly the same FV grid.

    Integrates inlet/outlet separately as additional states. No BE solution is
    used as truth; repeat at tighter tolerances to bound temporal contamination.
    """
    x = (np.arange(n)+.5)/n
    h0 = primary(name, x, start)
    def rhs(t, y):
        q = faces(name, y[:n], t)[0]
        cap = storage(name, y[:n])[1]
        return np.r_[-n*np.diff(q)/cap, flows(q)]
    pattern = lil_matrix((n+2, n+2))
    for i in range(n):
        pattern[i, max(0, i-1):min(n, i+2)] = 1
    pattern[n:, 0] = 1
    pattern[n:, n-1] = 1
    sol = solve_ivp(rhs, (start, end), np.r_[h0, 0., 0.], method="Radau",
                    rtol=rtol, atol=atol, jac_sparsity=pattern.tocsr())
    if not sol.success or sol.t[-1] != end or not np.all(np.isfinite(sol.y)):
        raise RuntimeError("method-of-lines reference incomplete")
    h = sol.y[:n, -1]
    change = (storage(name, h)[0]-storage(name, h0)[0]).sum()/n
    fi, fo = sol.y[n:, -1]
    return h, {"inflow": float(fi), "outflow": float(fo),
               "delta_storage": float(change), "mass_defect": float(change-fi+fo),
               "nfev": sol.nfev, "completed": True}
