"""Input-only small exact references; these are not solver acceptance results.

Contracts: interim/2026-11-16-issue-triage/todo.ja.md A-01/A-03/A-04/B-01/B-02;
interim/2026-11-13-joint-ext-combos/README.ja.md (pressure, work, equilibrium);
interim/2026-11-14-plate6-curved-offset/README.ja.md (geometry and strip formulae).
"""
import numpy as np


def tri6(r, s):
    a = 1-r-s
    n = np.array([a*(2*a-1), r*(2*r-1), s*(2*s-1), 4*a*r, 4*r*s, 4*s*a])
    dr = np.array([1-4*a, 4*r-1, 0., 4*(a-r), 4*s, -4*s])
    ds = np.array([1-4*a, 0., 4*s-1, -4*r, 4*r, 4*(a-s)])
    return n, dr, ds


def triangle_rule(order=5):
    """Independent tensor Gauss rule on 0<=r,s,r+s<=1."""
    x, w = np.polynomial.legendre.leggauss(order)
    x, w = (x+1)/2, w/2
    for i, r in enumerate(x):
        for j, v in enumerate(x):
            yield r, (1-r)*v, w[i]*w[j]*(1-r)


def hydrostatic_triangle(points, level, gamma):
    """All six consistent nodal forces, fully submerged quadratic straight face."""
    points = np.asarray(points)
    loads = np.zeros((6, 3))
    for r, s, w in triangle_rule():
        n, dr, ds = tri6(r, s)
        area_vector = np.cross(dr@points, ds@points)
        pressure = gamma*(level-(n@points)[2])
        assert pressure >= 0, 'this reference is for fully submerged faces'
        loads += w*n[:, None]*pressure*area_vector
    return loads


def darcy(z, bottom_head, top_head, height, conductivity, gamma):
    """Saturated 1D: total head linear; upward flux=-K dH/dz."""
    head = bottom_head+(top_head-bottom_head)*np.asarray(z)/height
    return head, gamma*(head-z), -conductivity*(top_head-bottom_head)/height


def diffusion(z, time, height, conductivity, storage, amplitude):
    """Saturated Richards perturbation about hydrostatic H0.

    S H_t=K H_zz, zero perturbation at both ends. Requires positive pressure
    throughout; H0 is added by the caller. No claim about unsaturated VG flow.
    Returns perturbation, time derivative, upward flux perturbation.
    """
    wave = np.pi/height
    decay = np.exp(-conductivity/storage*wave**2*time)
    h = amplitude*np.sin(wave*np.asarray(z))*decay
    return h, -conductivity/storage*wave**2*h, -conductivity*wave*amplitude*np.cos(wave*np.asarray(z))*decay


def strip(x, young, thickness, offset, axial_force, length, support):
    """nu=0, unit width: README(b) kappa=eN/D, eps_s=N/(Et)+e*kappa."""
    d = young*thickness**3/12
    curvature = offset*axial_force/d
    strain = axial_force/(young*thickness)+offset*curvature
    x = np.asarray(x)
    if support == 'cantilever':
        w, theta = -curvature*x*x/2, curvature*x
    elif support == 'simple':
        w, theta = curvature*x*(length-x)/2, curvature*(x-length/2)
    else:
        raise ValueError(support)
    return strain*x, w, theta, curvature
