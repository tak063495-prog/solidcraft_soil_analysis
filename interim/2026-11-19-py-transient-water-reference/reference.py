"""Closed/semi-closed 1-D references. No application imports or solver calls.

Units: m, h, m/h. Coordinates, equations and sources: README.ja.md.
"""
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.optimize import brentq


@dataclass(frozen=True)
class Case:
    name: str
    start: float
    end: float
    head_kind: str
    coordinate: str


CASES = {
    "terzaghi": Case("terzaghi", .02, .2, "excess_pressure_head", "down_from_drain"),
    "moving_level": Case("moving_level", 0., .2, "total_head", "horizontal_from_left"),
    "richards": Case("richards", .5, 5., "total_head", "up_from_water_table"),
}
KS, ALPHA, THETA_R, DELTA_THETA = .01, 1., .2, .25
RA, RB = .1, .9
D_RICHARDS = KS / (ALPHA * DELTA_THETA)
OMEGA, AMPLITUDE = 2 * np.pi, .1
MV, GAMMA_W = .001, 9.81  # 1/kPa, kPa/m; cv=1 m²/h
SS_TERZAGHI = MV*GAMMA_W  # 1/m; k=cv*Ss=0.00981 m/h


def validate(name, x, t):
    if name not in CASES:
        raise ValueError("unknown case")
    x = np.asarray(x, dtype=float)
    if not np.all(np.isfinite(x)) or not np.isfinite(t):
        raise ValueError("non-finite coordinate/time")
    if np.any((x < 0) | (x > 1)) or t < 0:
        raise ValueError("reference domain is 0<=x<=1, t>=0")
    return x


@lru_cache(maxsize=8)
def richards_modes(terms=128):
    """Robin eigenvalues and independently derived projection coefficients.

    k cos(k) + beta sin(k)=0; integral of 2*(RA-RB)*sinh(beta*z)*sin(k*z).
    Original paper's printed series is NOT transcribed; see derivation in README.
    """
    beta = ALPHA / 2
    k = np.array([brentq(lambda v: v*np.cos(v)+beta*np.sin(v),
                        (j+.5)*np.pi, (j+1)*np.pi, xtol=5e-15)
                  for j in range(terms)])
    norm = .5 - np.sin(2*k)/(4*k)
    b = 2*(RA-RB)*beta*np.exp(beta)*np.sin(k)/((beta**2+k**2)*norm)
    k.setflags(write=False)
    b.setflags(write=False)
    return k, b


def richards_u(x, t, terms=128):
    x = validate("richards", x, t)
    if t == 0:
        return RA+(1-RA)*np.exp(-ALPHA*x)
    k, b = richards_modes(terms)
    modes = np.sin(np.outer(x.ravel(), k))
    decay = b*np.exp(-D_RICHARDS*(k*k+ALPHA**2/4)*t)
    return RB+(1-RB)*np.exp(-ALPHA*x) + (
        np.exp(-ALPHA*x/2)*(modes@decay).reshape(x.shape))


def fields(name, x, t, terms=128):
    """Return comparison head and instantaneous flux in increasing x direction."""
    x = validate(name, x, t)
    if name == "terzaghi":
        if t == 0:
            raise ValueError("instantaneous drainage flux is singular at t=0")
        mu = (np.arange(terms)+.5)*np.pi
        decay = np.exp(-mu*mu*t)
        head = (np.sin(np.outer(x.ravel(), mu))@(2*decay/mu)).reshape(x.shape)
        flux = -SS_TERZAGHI*(np.cos(np.outer(x.ravel(), mu))@(2*decay)).reshape(x.shape)
    elif name == "moving_level":
        k = np.sqrt(1j*OMEGA)
        c = -1j*AMPLITUDE*np.exp(1j*OMEGA*t)/np.sinh(k)
        head = 1+np.real(c*np.sinh(k*(1-x)))
        flux = np.real(c*k*np.cosh(k*(1-x)))
    else:
        u = richards_u(x, t, terms)
        if np.any(u <= 0) or np.any(u > 1+1e-12):
            raise ValueError("exponential unsaturated model outside 0<U<=1")
        head = np.log(u)/ALPHA+x
        if t == 0:
            flux = np.full_like(x, -KS*RA)
        else:
            k, b = richards_modes(terms)
            decay = b*np.exp(-D_RICHARDS*(k*k+ALPHA**2/4)*t)
            xx = np.outer(x.ravel(), k)
            flux = -KS*RB-KS/ALPHA*np.exp(-ALPHA*x/2)*(
                (k*np.cos(xx)+ALPHA/2*np.sin(xx))@decay).reshape(x.shape)
    return head, flux


def primary(name, x, t):
    """FV unknown: pressure head psi for Richards, comparison head otherwise."""
    h, _ = fields(name, x, t)
    return h-x if name == "richards" else h


def storage(name, h):
    h = np.asarray(h)
    if name == "richards":
        u = np.exp(ALPHA*h)
        return THETA_R+DELTA_THETA*u, ALPHA*DELTA_THETA*u
    if name == "terzaghi":
        return SS_TERZAGHI*h, np.full_like(h, SS_TERZAGHI)
    return h, np.ones_like(h)


def boundary(name, t):
    """Left Dirichlet; right (Dirichlet or flux) in increasing coordinate."""
    if name == "terzaghi":
        return 0., ("flux", 0.)
    if name == "moving_level":
        return 1+AMPLITUDE*np.sin(OMEGA*t), ("head", 1.)
    return 0., ("flux", -KS*RB)


def metadata(name):
    c = CASES[name]
    return {"schema": "transient-water-1", "case": name,
            "benchmark": "fixed-v1", "length_unit": "m", "time_unit": "h",
            "head_unit": "m", "flux_unit": "m/h", "flux_kind": "instantaneous",
            "head_kind": c.head_kind, "coordinate": c.coordinate,
            "flux_sign": "positive_increasing_coordinate",
            "complete": True, "run_status": "completed"}
