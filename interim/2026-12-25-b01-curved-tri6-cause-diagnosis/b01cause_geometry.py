"""Observation of quadratic geometry, never a curved shell constitutive element."""
import math
import numpy as np

GRAD = np.array([[-1., -1.], [1., 0.], [0., 1.]])
EDGES = ((0, 1), (1, 2), (2, 0))
PARAM = np.array([[0., 0.], [1., 0.], [0., 1.], [.5, 0.], [.5, .5], [0., .5]])
SAMPLES = np.vstack((PARAM, [[1/3, 1/3], [.091576213509771, .091576213509771],
                            [.816847572980459, .091576213509771],
                            [.091576213509771, .816847572980459],
                            [.445948490915965, .445948490915965]]))


def dot(a, b):
    a, b = np.asarray(a, float).ravel(order='C'), np.asarray(b, float).ravel(order='C')
    if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('finite equally sized components required')
    return math.fsum(float(x)*float(y) for x, y in zip(a, b))


def norm(a):
    return math.sqrt(dot(a, a))


def matvec(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return np.array([dot(row, b) for row in a])


def metric(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.shape != b.shape:
        raise ValueError('shape mismatch')
    absolute, denominator = norm(a-b), norm(b)
    return dict(absolute_l2=absolute, reference_l2=denominator,
                relative_l2=absolute/max(denominator, 1e-12),
                maximum_absolute=float(np.max(np.abs(a-b), initial=0.)))


def shape(r, s):
    l = np.array([1-r-s, r, s])
    n = list(l*(2*l-1))
    d = [(4*l[i]-1)*GRAD[i] for i in range(3)]
    dd = [4*np.outer(GRAD[i], GRAD[i]) for i in range(3)]
    for i, j in EDGES:
        n.append(4*l[i]*l[j])
        d.append(4*(l[i]*GRAD[j]+l[j]*GRAD[i]))
        dd.append(4*(np.outer(GRAD[i], GRAD[j])+np.outer(GRAD[j], GRAD[i])))
    return np.array(n), np.array(d), np.array(dd)


def frame(points, r, s):
    x = np.asarray(points, float)
    n, d, dd = shape(r, s)
    tangent = d.T@x
    second = np.einsum('iab,ic->abc', dd, x)
    c = np.cross(tangent[0], tangent[1])
    jac = norm(c)
    if not jac > 0:
        raise ValueError('rank deficient sampled geometry')
    normal = c/jac
    dc = [np.cross(second[0, a], tangent[1])+np.cross(tangent[0], second[1, a]) for a in range(2)]
    dn = np.array([(v-normal*dot(normal, v))/jac for v in dc])
    return dict(position=n@x, tangent=tangent, metric=tangent@tangent.T,
                jacobian=jac, normal=normal, normal_derivative=dn,
                curvature=np.einsum('abc,c->ab', second, normal))


def geometry(points):
    x = np.asarray(points, float)
    edge = x[1]-x[0]
    normal = np.cross(edge, x[2]-x[0]); normal /= norm(normal)
    scale = max(norm(x[i]-x[j]) for i, j in EDGES)
    cutoff = 64*np.finfo(float).eps*scale
    delta = np.array([x[3+e]-(x[i]+x[j])*.5 for e, (i, j) in enumerate(EDGES)])
    frames = [frame(x, *p) for p in SAMPLES]
    # An exact rigid displacement has sym(a_alpha dot u_,beta)=0 on this geometry.
    rigid_strains = []
    for omega in np.eye(3):
        u = np.cross(np.tile(omega, (6, 1)), x)
        for p, f in zip(SAMPLES, frames):
            _, d, _ = shape(*p); du = d.T@u
            e = f['tangent']@du.T
            rigid_strains.extend(((e+e.T)*.5).ravel())
    return dict(scale_m=scale, guard_cutoff_m=cutoff, midpoint_delta_m=delta.tolist(),
                midpoint_distances_m=[norm(v) for v in delta],
                out_of_corner_plane_m=[dot(v, normal) for v in delta],
                sampled_area_jacobian_min=min(f['jacobian'] for f in frames),
                sampled_area_jacobian_max=max(f['jacobian'] for f in frames),
                maximum_normal_departure= max(norm(f['normal']-normal) for f in frames),
                maximum_normal_derivative=max(norm(f['normal_derivative']) for f in frames),
                maximum_second_form=max(norm(f['curvature']) for f in frames),
                exact_quadratic_rigid_covariant_membrane_max=max(map(abs, rigid_strains)),
                sample_parameters=SAMPLES.tolist(),
                samples=[{k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in f.items()} for f in frames],
                sampled_only=True, global_jacobian_positivity_certified=False)


def rigid(points):
    x = np.asarray(points, float)
    center = np.array([math.fsum(x[:, a])/len(x) for a in range(3)])
    modes = []
    for axis in np.eye(3):
        modes.append(np.tile(np.r_[axis, [0., 0., 0.]], (len(x), 1)).ravel())
    for axis in np.eye(3):
        modes.append(np.column_stack((np.cross(np.tile(axis, (len(x), 1)), x-center),
                                      np.tile(axis, (len(x), 1)))).ravel())
    return np.column_stack(modes)


def surface(kind, uv):
    u, v = np.asarray(uv).T
    if kind == 'plane':
        return np.column_stack((u, v, 0*u))
    if kind == 'cylinder':
        return np.column_stack((u, 2*np.sin(v/2), 2*(1-np.cos(v/2))))
    if kind == 'sphere':
        return np.column_stack((2*np.sin(u/2)*np.cos(v/2), 2*np.sin(v/2),
                                2*np.cos(u/2)*np.cos(v/2)-2))
    raise ValueError(kind)


def mesh(kind, n):
    """n=1/2 only in this study: 2/8 Tri6, shared quadratic edge nodes."""
    if n not in (1, 2):
        raise ValueError('this study fixes n=1/2; review resources before new meshes')
    uv = np.array([(i/(2*n), j/(2*n)) for i in range(2*n+1) for j in range(2*n+1)])
    def node(i, j):
        return i*(2*n+1)+j
    elements = []
    for i in range(n):
        for j in range(n):
            for tri in [((2*i, 2*j), (2*i+2, 2*j), (2*i+2, 2*j+2)),
                        ((2*i, 2*j), (2*i+2, 2*j+2), (2*i, 2*j+2))]:
                mids = [((tri[a][0]+tri[b][0])//2, (tri[a][1]+tri[b][1])//2) for a, b in EDGES]
                elements.append([node(*p) for p in [*tri, *mids]])
    boundary = [i for i, p in enumerate(uv) if any(v in (0., 1.) for v in p)]
    return surface(kind, uv), elements, boundary


def document_text(path):
    # universal newline is deliberate: suffix contracts must survive autocrlf.
    with open(path, encoding='utf-8', newline=None) as f:
        return f.read()
