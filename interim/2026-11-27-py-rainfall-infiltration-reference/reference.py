"""Independent conservative 1D Richards diagnostics; NumPy/SciPy only.

Depth x and Darcy flux q positive downward. See README for the derivation.
No production solver imports, subprocesses, or production-output reads.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import quad
from scipy.special import erfc, erfcx
from scipy.sparse import diags, eye
from scipy.sparse.linalg import splu
from scipy.linalg import solve_banded


@dataclass(frozen=True)
class Gardner:
    ks: float = 1e-4
    alpha: float = 1.0
    theta_r: float = .08
    theta_s: float = .45
    initial_u: float = .2
    length: float = 4.0

    @property
    def delta(self):
        return self.theta_s - self.theta_r

    @property
    def v(self):
        return self.ks / self.delta

    @property
    def diffusion(self):
        return self.v / self.alpha


def step_u(x, t, jump, p=Gardner(), eps=2e-11):
    """Semi-infinite flux-step solution by an independently derived heat kernel.

    Substitute tau=y*y to remove the integrable surface t=0 singularity.
    quad's absolute/relative tolerances concern quadrature, not acceptance.
    """
    if t <= 0:
        return 0.0
    d, v = p.diffusion, p.v
    c = v / (2*np.sqrt(d))

    def integrand(y):
        if y == 0:
            return 2/np.sqrt(np.pi) if x == 0 else 0.0
        w = (x+v*y*y)/(2*np.sqrt(d)*y)
        return np.exp(-(x-v*y*y)**2/(4*d*y*y)) * (
            2/np.sqrt(np.pi) - 2*y*c*erfcx(w))

    return jump/(p.delta*np.sqrt(d))*quad(
        integrand, 0, np.sqrt(t), epsabs=eps, epsrel=eps, limit=150)[0]


def step_q(x, t, jump, p=Gardner()):
    if t <= 0:
        return 0.0
    d, v = p.diffusion, p.v
    a = (x-v*t)/(2*np.sqrt(d*t))
    b = (x+v*t)/(2*np.sqrt(d*t))
    return jump*.5*(erfc(a)+np.exp(-a*a)*erfcx(b))


def analytical(x, t, p=Gardner(), jumps=((0., 2e-5), (300., -2e-5))):
    x = np.asarray(x, dtype=float)
    u = np.full(x.shape, p.initial_u)
    q = np.full(x.shape, p.ks*p.initial_u)
    for at, jump in jumps:
        u += np.array([step_u(float(a), t-at, jump, p) for a in x.flat]).reshape(x.shape)
        q += np.array([step_q(float(a), t-at, jump, p) for a in x.flat]).reshape(x.shape)
    return u, q


def gardner_fv(n, dt, end=600., p=Gardner(), method='BE',
               jumps=((0., 2e-5), (300., -2e-5))):
    """Cell-centred FV in u=exp(alpha*psi), exact Gardner storage.

    Central Darcy faces, prescribed rainfall flux, deep initial-u boundary.
    Each step stops at a pulse edge; trapezoidal face flux for CN.
    """
    if n < 8 or n % 8 or dt <= 0 or end <= 0 or method not in ('BE', 'CN'):
        raise ValueError('invalid resolution, time or method')
    dx = p.length/n
    x = (np.arange(n)+.5)*dx
    u = np.full(n, p.initial_u)
    # q_i = a*u_left + b*u_right at internal downward faces.
    a, b = p.ks*(.5+1/(p.alpha*dx)), p.ks*(.5-1/(p.alpha*dx))
    lower = np.full(n-1, a/(p.delta*dx))
    upper = np.full(n-1, -b/(p.delta*dx))
    diagonal = np.full(n, (b-a)/(p.delta*dx))
    diagonal[0] = -a/(p.delta*dx)
    # Bottom boundary u=initial_u at half-cell distance.
    bot_a, bot_b = p.ks*(.5+2/(p.alpha*dx)), p.ks*(.5-2/(p.alpha*dx))
    diagonal[-1] = (b-bot_a)/(p.delta*dx)
    matrix = diags([lower, diagonal, upper], [-1, 0, 1], format='csc')
    weight = 1. if method == 'BE' else .5
    t, cumulative_in, cumulative_out, maximum_mass = 0., 0., 0., 0.
    recharge_excess = 0.
    recharge_face = n//8  # x=L/8=.5 m for all study meshes (multiples of 8).
    history, cache = [], {}
    while t < end-1e-10:
        next_edge = min([at for at, _ in jumps if at > t+1e-9]+[end])
        step = min(dt, end-t, next_edge-t)
        rain = p.ks*p.initial_u + sum(j for at, j in jumps if at <= t+1e-9)
        forcing = np.zeros(n)
        forcing[0] = rain/(p.delta*dx)
        forcing[-1] -= bot_b*p.initial_u/(p.delta*dx)
        if step not in cache:
            cache[step] = splu(eye(n, format='csc')-weight*step*matrix)
        old = u.copy()
        u = cache[step].solve(old+(1-weight)*step*(matrix@old)+step*forcing)
        if np.any(u <= 0) or np.any(u >= 1):
            raise ValueError('Gardner unsaturated branch exceeded')
        q = np.r_[rain, a*u[:-1]+b*u[1:], bot_a*u[-1]+bot_b*p.initial_u]
        old_bottom = bot_a*old[-1]+bot_b*p.initial_u
        old_recharge = a*old[recharge_face-1]+b*old[recharge_face]
        recharge_excess += step*(weight*q[recharge_face]+(1-weight)*old_recharge-p.ks*p.initial_u)
        cumulative_in += step*rain
        cumulative_out += step*(weight*q[-1]+(1-weight)*old_bottom)
        mass = p.delta*dx*np.sum(u-p.initial_u)-cumulative_in+cumulative_out
        maximum_mass = max(maximum_mass, abs(mass))
        t += step
        history.append([t, float(q[-1]), cumulative_in, cumulative_out, float(mass), recharge_excess])
    return dict(x=x, u=u, psi=np.log(u)/p.alpha, q=q, history=np.asarray(history),
                max_mass_m=maximum_mass, dt=dt, n=n, end=end, method=method)


def vg(psi):
    """Original rainfall VG/Mualem law; theta and specific-storage offset.

    van Genuchten (1980), m=1-1/n, Mualem pore exponent 1/2.
    """
    psi = np.asarray(psi)
    m, n, alpha, ks = 1-1/1.6, 1.6, 1.2, 1e-6
    s = (1+(alpha*np.maximum(-psi, 0))**n)**(-m)
    theta = .08+.37*s
    k = ks*np.sqrt(s)*(1-(1-s**(1/m))**m)**2
    storage = theta-.45+1e-5*np.maximum(psi, 0)
    return theta, k, storage


def vg_column(n=24, dt=100., end=10800.):
    """Uniform-top 1D analogue, NOT the original finite-square 3D problem.

    Own damped Newton, numerical tridiagonal Jacobian, actual storage change.
    Atmospheric complementarity uses zero-ponding capacity; excess is rejected.
    """
    if n < 4 or dt <= 0 or not 0 < end <= 10800:
        raise ValueError('invalid column configuration')
    dx = 1/n
    x = (np.arange(n)+.5)*dx
    psi = x-1
    initial = vg(psi)[2]
    net, maximum_mass, rejected, frames = 0., 0., 0., []
    t, max_residual = 0., 0.

    def fields(h, rain):
        theta, k, storage = vg(h)
        harmonic = 2*k[:-1]*k[1:]/(k[:-1]+k[1:])
        q = np.empty(n+1)
        q[1:-1] = harmonic*(1-np.diff(h)/dx)
        capacity = 2*k[0]*1e-6/(k[0]+1e-6)*(1-h[0]/(dx/2))
        q[0] = min(rain, capacity)
        q[-1] = 2*k[-1]*1e-6/(k[-1]+1e-6)*(1+h[-1]/(dx/2))
        return theta, storage, q

    while t < end-1e-8:
        edge = min([s for s in (3600., 7200., 10800.) if s > t+1e-8]+[end])
        step = min(dt, end-t, edge-t)
        rain = 3e-6 if t < 3600.-1e-8 or t >= 7200.-1e-8 else 0.
        old_storage = vg(psi)[2]

        def residual(h):
            _, storage, q = fields(h, rain)
            return (storage-old_storage)*dx/step+np.diff(q)

        h = psi.copy()
        for iteration in range(35):
            r = residual(h)
            if np.max(np.abs(r)) < 2e-13:
                break
            jac = np.zeros((3, n))
            epsilon = 1e-6
            for j in range(n):
                hp, hm = h.copy(), h.copy()
                hp[j] += epsilon
                hm[j] -= epsilon
                column = (residual(hp)-residual(hm))/(2*epsilon)
                for i in range(max(0, j-1), min(n, j+2)):
                    jac[1+i-j, j] = column[i]
            correction = solve_banded((1, 1), jac, -r)
            fraction = 1.
            while fraction >= 2**-16:
                trial = h+fraction*correction
                if np.linalg.norm(residual(trial)) < np.linalg.norm(r):
                    h = trial
                    break
                fraction /= 2
            else:
                raise RuntimeError('VG Newton line search failed')
        else:
            raise RuntimeError('VG Newton iteration limit; not accepted')
        max_residual = max(max_residual, float(np.max(np.abs(residual(h)))))
        psi = h
        theta, storage, q = fields(h, rain)
        net += step*(q[0]-q[-1])
        rejected += step*(rain-q[0])
        maximum_mass = max(maximum_mass, abs(float(dx*np.sum(storage-initial)-net)))
        t += step
        if any(abs(t-s) < 1e-7 for s in (3600., 7200., 10800., end)):
            frames.append(dict(time_s=t, head_m=(psi+1-x).tolist(),
                               theta=theta.tolist(), storage_offset_m=float(dx*sum(storage)),
                               net_infiltration_m=net, rejected_rain_m=rejected,
                               top_flux_down_m_s=float(q[0]), bottom_flux_down_m_s=float(q[-1])))
    return dict(n=n, dt_s=dt, x_m=x.tolist(), frames=frames,
                max_mass_m=maximum_mass, max_residual_m_s=max_residual,
                accepted=None, original_3d_gate_value=None)


def mechanics(z, water_level, gamma=9.81, total=50., modulus=10000.):
    """Fixed total compressive stress, Terzaghi alpha_B=1, small strain column.

    Saturated receptor only: negative pore pressure is not extrapolated here.
    """
    pressure = gamma*(water_level-np.asarray(z))
    if np.any(pressure < 0):
        raise ValueError('mechanical receptor must remain submerged')
    effective = total-pressure
    return effective, effective/modulus


def mechanical_fv(n, water_level, gamma=9.81, total=50., modulus=10000.):
    """Own 1D linear FE equilibrium; unit area, bottom fixed, top total load.

    Receptor elevation z=-y, y in [0,1]; downward displacement positive.
    Hydrostatic cell pressure is integrated exactly, constant strain elements.
    """
    dy = 1/n
    y = np.arange(n+1)*dy
    pressure = gamma*(water_level+(y[:-1]+y[1:])/2)
    stiffness = np.zeros((n+1, n+1))
    pore_force = np.zeros(n+1)
    for i in range(n):
        stiffness[i:i+2, i:i+2] += modulus/dy*np.array([[1., -1.], [-1., 1.]])
        pore_force[i:i+2] += pressure[i]*np.array([1., -1.])
    external = np.zeros(n+1)
    external[0] = total
    displacement = np.r_[np.linalg.solve(stiffness[:-1, :-1], (external-pore_force)[:-1]), 0.]
    stress = -modulus*np.diff(displacement)/dy
    exact = ((total-gamma*water_level)*(1-y)-gamma*(1-y*y)/2)/modulus
    residual = stiffness@displacement+pore_force-external
    return dict(displacement=displacement, effective_stress=stress,
                exact_displacement=exact,
                max_displacement_error=float(np.max(np.abs(displacement-exact))),
                max_stress_error=float(np.max(np.abs(stress-(total-pressure)))),
                equilibrium_residual=float(np.max(np.abs(residual[:-1]))))
