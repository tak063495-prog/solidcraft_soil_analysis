"""Independent prescribed-jump replay. No equilibrium solver or native imports."""
import numpy as np


def dot3(a, b):
    # Explicit x,y,z order, no BLAS/parallel reduction.
    return (a[..., 0]*b[..., 0] + a[..., 1]*b[..., 1]) + a[..., 2]*b[..., 2]


def energy(q, normal, kn, kt):
    qn = dot3(q, normal)
    tau = q - qn[:, None]*normal
    return qn*qn/(2*kn) + dot3(tau, tau)/(2*kt)


def integrate(start, traction, target, normal, kn, kt, mu, counts):
    """All GP, explicit fixed per-GP subdivisions, exact planar opening roots.

    Returns traction, plastic increment, D, endpoint work, numerical loss.
    It never calls a global equilibrium solver.
    """
    start, target = np.asarray(start), np.asarray(target)
    jump, q = start.copy(), np.asarray(traction).copy()
    normal, kn, kt, mu = map(np.asarray, (normal, kn, kt, mu))
    counts = np.broadcast_to(np.asarray(counts, dtype=int), (len(jump),))
    if np.any(counts < 1) or np.any(counts > 50000):
        raise ValueError('invalid diagnostic count')
    dp = np.zeros_like(jump)
    diss = np.zeros(len(jump)); work = diss.copy(); numerical = diss.copy()
    moving = np.any(target != start, axis=1)

    def closed(ids, end):
        if not len(ids):
            return
        n = normal[ids]; k = kt[ids]
        delta = end-jump[ids]
        oldq = q[ids]
        trial = oldq-dot3(oldq,n)[:,None]*n+k[:,None]*(delta-dot3(delta,n)[:,None]*n)
        qn = kn[ids]*np.minimum(dot3(end,n), 0.)
        cap = -mu[ids]*qn
        length = np.sqrt(dot3(trial,trial))
        scale = np.ones(len(ids))
        positive = length > 0
        scale[positive] = np.minimum(1.,cap[positive]/length[positive])
        tau = scale[:,None]*trial
        plastic = (trial-tau)/k[:,None]
        nextq = qn[:,None]*n+tau
        de = energy(nextq,n,kn[ids],k)-energy(oldq,n,kn[ids],k)
        dw = dot3(nextq,delta)
        dd = np.maximum(dot3(tau,plastic),0.)
        dp[ids] += plastic; diss[ids] += dd; work[ids] += dw
        numerical[ids] += np.maximum(dw-de-dd,0.)
        jump[ids] = end; q[ids] = nextq

    for i in range(1,int(max(counts))+1):
        ids = np.flatnonzero((counts >= i)&moving)
        end = start[ids]+(target[ids]-start[ids])*(i/counts[ids])[:,None]
        final = counts[ids] == i
        end[final] = target[ids[final]]
        n = normal[ids]
        g0, g1 = dot3(jump[ids],n),dot3(end,n)
        both = (g0 <= 0)&(g1 <= 0)
        closed(ids[both],end[both])
        crossing = (g0 > 0) != (g1 > 0)
        cross_ids = ids[crossing]
        if len(cross_ids):
            before = jump[cross_ids].copy(); dest = end[crossing]
            a = g0[crossing]/(g0[crossing]-g1[crossing])
            plane = before+(dest-before)*a[:,None]
            cn = normal[cross_ids]
            plane -= dot3(plane,cn)[:,None]*cn
            opening = g0[crossing] <= 0
            closed(cross_ids[opening],plane[opening])
            close_ids = cross_ids[~opening]
            jump[close_ids] = plane[~opening]; q[close_ids] = 0.
            closed(close_ids,dest[~opening])
        opened = g1 > 0
        jump[ids[opened]] = end[opened]; q[ids[opened]] = 0.
    return q, dp, diss, work, numerical
