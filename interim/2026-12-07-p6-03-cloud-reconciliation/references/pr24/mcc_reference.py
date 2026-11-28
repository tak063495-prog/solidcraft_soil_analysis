"""Independent axisymmetric MCC rates. Imports only stdlib, NumPy and SciPy.
Definitions and derivations: README.ja.md, source IDs S1--S4. No solver imports.
Compression-positive effective p,q; ev=log(v0/v), ea=ev/3+es.
"""
from dataclasses import dataclass
import math
import numpy as np
from scipy.optimize import brentq


@dataclass(frozen=True)
class Material:
    lam: float = .2
    kap: float = .04
    M: float = 1.2
    nu: float = .3
    N: float = 3.5
    pref: float = 1.

    def __post_init__(self):
        if not all(math.isfinite(x) for x in vars(self).values()):
            raise ValueError('nonfinite material')
        if not (self.lam > self.kap > 0 and self.M > 0 and -1 < self.nu < .5 and self.pref > 0):
            raise ValueError('invalid MCC parameters')

    @property
    def g_over_k(self):
        return 3 * (1 - 2 * self.nu) / (2 * (1 + self.nu))

    def volume(self, p, pc):
        if p <= 0 or pc < p:
            raise ValueError('initial p<=0 or pc<p')
        v = self.N - self.lam * math.log(pc / self.pref) + self.kap * math.log(pc / p)
        if v <= 1:
            raise ValueError('initial e<=0')
        return v


def initial(mat=Material(), p=100., ocr=1.):
    pc = p * ocr
    return np.array([p, 0., math.log(pc), mat.volume(p, pc), 0.])


def invariants(y):
    p, q, logpc, v, ev = map(float, y)
    return dict(p=p, q=q, pc=math.exp(logpc), e=v - 1, eps_v=ev)


def elastic(ea, y0, mode, mat):
    """Exact elastic axisymmetric path; y0 is an isotropic initial state."""
    p0, q0, lpc, v0, ev0 = y0
    if q0 != 0 or ev0 != 0:
        raise ValueError('elastic origin must be isotropic, zero total strain')
    if mode == 'undrained':
        ev = 0.
    elif mode == 'drained':
        ev = ea / (1 / 3 + 1 / mat.g_over_k)
    elif mode == 'k0':
        ev = ea
    else:
        raise ValueError('unknown mode')
    v = v0 * math.exp(-ev)
    p = p0 * math.exp((v0 - v) / mat.kap)
    if mode == 'undrained':
        q = 3 * mat.g_over_k * v0 * p0 / mat.kap * ea
    elif mode == 'drained':
        q = 3 * (p - p0)
    else:
        q = 2 * mat.g_over_k * (p - p0)
    return np.array([p, q, lpc, v, ev])


def yield_value(y, mat):
    p, q, lpc, _, _ = y
    pc = math.exp(lpc)
    return (q / pc)**2 + mat.M**2 * (p / pc) * (p / pc - 1)


def yield_start(y0, mode, end, mat):
    if abs(yield_value(y0, mat)) < 1e-14:
        return 0.
    if yield_value(elastic(end, y0, mode, mat), mat) <= 0:
        return None
    return brentq(lambda x: yield_value(elastic(x, y0, mode, mat), mat), 0, end,
                  xtol=1e-15, rtol=1e-14)


def rate(y, mode, mat):
    p, q, lpc, v, _ = y
    pc = math.exp(lpc)
    if p <= 0 or v <= 1 or not np.isfinite(y).all():
        raise ValueError('invalid rate state; no floors/projection')
    K = v * p / mat.kap
    G3 = 3 * mat.g_over_k * K
    A, B = mat.M**2 * (2 * p - pc), 2 * q
    z = np.array([K * A, G3 * B])
    hard = mat.M**2 * p * pc * v / (mat.lam - mat.kap) * A
    denom = K * A * A + G3 * B * B + hard
    if denom <= 0:
        raise ValueError('nonpositive consistency denominator')
    D = np.diag([K, G3]) - np.outer(z, z) / denom
    if mode == 'k0':
        ev, es = 1., 2 / 3
    elif mode == 'undrained':
        ev, es = 0., 1.
    elif mode == 'drained':
        base = np.array([1., 2 / 3])
        radial = np.array([2., -2 / 3])
        stress_radial = np.array([1., -1 / 3])
        r = -(stress_radial @ D @ base) / (stress_radial @ D @ radial)
        ev, es = base + r * radial
    else:
        raise ValueError('unknown mode')
    dmult = float(z @ np.array([ev, es]) / denom)
    if dmult < -1e-12:
        raise ValueError('unloading plastic branch unsupported')
    dp, dq = D @ np.array([ev, es])
    return np.array([dp, dq, v / (mat.lam - mat.kap) * A * dmult, -v * ev, ev])


def integrate(mode, n, end=.08, ocr=1., mat=Material(), samples=16):
    """Event-aligned explicit midpoint; n plastic substeps per sample interval.

    All levels have identical output abscissae. Exact elastic event location is
    shared, not an adaptive material error hidden inside a spatial refinement.
    """
    if not isinstance(n, int) or n < 1 or end <= 0 or samples < 1:
        raise ValueError('invalid refinement')
    y0 = initial(mat, ocr=ocr)
    event = yield_start(y0, mode, end, mat)
    times = np.linspace(0., end, samples + 1)
    rows = [dict(step=0, axial_strain=0., **invariants(y0))]
    y = y0.copy()
    previous = 0.
    for i, t in enumerate(times[1:], 1):
        if event is None or t <= event:
            y = elastic(float(t), y0, mode, mat)
        else:
            start = max(previous, event)
            if previous <= event:
                y = elastic(event, y0, mode, mat)
            h = (t - start) / n
            for j in range(n):
                k1 = rate(y, mode, mat)
                mid = y + .5 * h * k1
                if mode == 'k0':
                    mid[4] = start + (j+.5)*h
                elif mode == 'undrained':
                    mid[4] = 0.
                mid[3] = y0[3] * math.exp(-mid[4])
                y = y + h * rate(mid, mode, mat)
                if mode == 'k0':
                    y[4] = start + (j+1)*h
                elif mode == 'undrained':
                    y[4] = 0.
                y[3] = y0[3] * math.exp(-y[4])
        rows.append(dict(step=i, axial_strain=float(t), **invariants(y)))
        previous = t
    return rows, event


def isotropic(pressures, mat=Material(), p0=100., pc0=100.):
    """Exact hydrostatic loading/unloading, v-lnp lines, ev=ln(v0/v)."""
    v0 = mat.volume(p0, pc0)
    p, pc, v = p0, pc0, v0
    rows = []
    for target in pressures:
        if not math.isfinite(target) or target <= 0:
            raise ValueError('nonpositive pressure')
        if target > pc:
            v -= mat.kap * math.log(pc / p) + mat.lam * math.log(target / pc)
            pc = target
        else:
            v -= mat.kap * math.log(target / p)
        p = target
        if v <= 1:
            raise ValueError('e<=0')
        rows.append(dict(p=p, q=0., pc=pc, e=v-1, eps_v=math.log(v0/v)))
    return rows


def richardson(values, ratio=2., min_order=.5, max_order=8., max_contraction=.7):
    """Reference-free vector differences. Guards from S5; no pass threshold.

    Order is null for zero/roundoff differences or non-asymptotic sequences.
    No order is inferred from a known reference error.
    """
    a, b, c = (np.asarray(x, dtype=float) for x in values)
    if a.shape != b.shape or b.shape != c.shape or ratio <= 1:
        raise ValueError('incompatible refinement series')
    if not all(np.isfinite(x).all() for x in [a, b, c]):
        raise ValueError('nonfinite refinement series')
    if not a.size or not math.isfinite(ratio):
        raise ValueError('empty/nonfinite refinement series')
    d1, d2 = np.linalg.norm(b-a), np.linalg.norm(c-b)
    result = dict(observed_order=None, extrapolated=None, fine_error_estimate=None,
                  difference_coarse=float(d1), difference_fine=float(d2), raw_observed_order=None,
                  difference_direction_cosine=None, accepted=None,
                  tolerance='未定義')
    scale = max(np.linalg.norm(a), np.linalg.norm(b), np.linalg.norm(c), 1.)
    if min(d1, d2) <= 128 * np.finfo(float).eps * scale:
        result['reason'] = 'roundoff_or_exact; order 未定義'
        return result
    result['raw_observed_order'] = math.log(d1/d2) / math.log(ratio)
    result['difference_direction_cosine'] = float(np.vdot(b-a,c-b)/(d1*d2))
    if np.vdot(b-a, c-b) <= 0:
        result['reason'] = 'direction_reversal; order 未確認'
        return result
    order = math.log(d1/d2) / math.log(ratio)
    if not min_order <= order <= max_order or d2 / d1 >= max_contraction:
        result['reason'] = 'non_asymptotic; order 未確認'
        return result
    estimate = c + (c-b) / (ratio**order-1)
    result.update(observed_order=order, extrapolated=estimate.tolist(),
                  fine_error_estimate=float(np.linalg.norm(estimate-c)), reason='diagnostic_only')
    return result


def mc_column(z, stage):
    """S6 exact original MC column, not MCC. z=height above fixed base, m."""
    z = np.asarray(z, dtype=float)
    H, E, nu, gamma, P, cohesion, phi = .2, 1000., .25, 20., 12., 2., math.radians(10)
    if np.any((z < 0) | (z > H)):
        raise ValueError('outside column')
    K = E / (3 * (1 - 2 * nu))
    C = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))
    A = 3 * (1 + math.sin(phi)) / (3 - math.sin(phi))
    B = 4 * cohesion * math.cos(phi) / (3 - math.sin(phi))
    weight = gamma * (H - z)
    if stage == 'gravity':
        u = -gamma * (H*z-z*z/2) / C
        sv = weight
        sh = nu / (1-nu) * sv
    elif stage in ('load', 'unload'):
        u = -((P+gamma*H-B)*z-gamma*z*z/2) / (A*K)
        sv = P + weight
        mean = (sv-B)/A
        sh = (3*mean-sv)/2
        if stage == 'unload':
            u += P*z/C
            sv -= P
            sh -= nu*P/(1-nu)
    else:
        raise ValueError('unknown stage')
    return dict(z=z.tolist(), uz=u.tolist(), sigma_v=sv.tolist(), sigma_h=sh.tolist())
