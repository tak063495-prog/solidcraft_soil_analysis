"""Input-only scalar Coulomb sweeping solution, no solver/probe imports."""
import math


def advance_exact(state, target, kn, kt, mu):
    """Integrate affine (x, gap) exactly; open intervals have no shear memory."""
    x0, g0 = state['jump']
    x1, g1 = target
    out = dict(state)
    out['jump'] = list(target)
    if g0 > 0 and g1 > 0:
        out['tau'] = 0.0
        return out
    if (g0 > 0) != (g1 > 0):
        a = g0 / (g0 - g1)
        xc = x0 + a * (x1 - x0)
        if g0 > 0:
            base = dict(state, jump=[xc, 0.0], tau=0.0)
            return advance_exact(base, target, kn, kt, mu)
        closed = advance_exact(state, [xc, 0.0], kn, kt, mu)
        return dict(closed, jump=list(target), tau=0.0)
    cap0, cap1 = -mu * kn * min(g0, 0.0), -mu * kn * min(g1, 0.0)
    dx = x1 - x0
    tau0 = state['tau']
    trial = tau0 + kt * dx
    tau1 = max(-cap1, min(cap1, trial))
    ds = dx - (tau1 - tau0) / kt
    out['tau'] = tau1
    out['slip'] += ds
    out['accumulated'] += abs(ds)
    if trial > cap1 or trial < -cap1:
        sign = 1.0 if trial > cap1 else -1.0
        rate = sign * kt * dx - (cap1 - cap0)
        onset = (cap0 - sign * tau0) / rate
        # Affine admissible bounds imply the first root lies in [0,1].
        if not -1e-10 <= onset <= 1 + 1e-10:
            raise ArithmeticError(('yield root outside segment', onset))
        onset = max(0.0, min(1.0, onset))
        cap_onset = cap0 + onset * (cap1 - cap0)
        out['dissipation'] += 0.5 * (cap_onset + cap1) * abs(ds)
    return out


def initial():
    return dict(jump=[0.0, 0.0], tau=0.0, slip=0.0, accumulated=0.0,
                dissipation=0.0)


def refine(vertices, n):
    result = []
    start = [0.0, 0.0]
    for segment, end in enumerate(vertices):
        for i in range(1, n + 1):
            q = list(end) if i == n else [a + (b-a)*i/n for a,b in zip(start,end)]
            result.append((segment, i, q))
        start = end
    return result
