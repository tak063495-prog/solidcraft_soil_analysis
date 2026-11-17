"""Opt-in arbitrary precision observation of the existing Coulomb branches.

No solver, state adoption, regularization, or production import. The three
derivatives are forward analytic derivatives, including crossing and substeps.
mpmath is imported only by this explicitly run diagnostic, never by CI tests.
"""
from mpmath import mp


class D:
    def __init__(self, value, derivative=None):
        self.v = mp.mpf(value)
        self.d = mp.zeros(1, 3) if derivative is None else derivative

    def __add__(self, other):
        other = dual(other)
        return D(self.v + other.v, self.d + other.d)

    __radd__ = __add__

    def __neg__(self):
        return D(-self.v, -self.d)

    def __sub__(self, other):
        return self + -dual(other)

    def __rsub__(self, other):
        return dual(other) + -self

    def __mul__(self, other):
        other = dual(other)
        return D(self.v * other.v, self.d * other.v + other.d * self.v)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = dual(other)
        return D(self.v / other.v,
                 (self.d * other.v - other.d * self.v) / other.v**2)

    def sqrt(self):
        value = mp.sqrt(self.v)
        return D(value, self.d / (2 * value))


def dual(value):
    return value if isinstance(value, D) else D(value)


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), D(0))


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def scale(a, scalar):
    return [x * scalar for x in a]


def point(normal, kn, kt, mu, before, oldq, rest, target, count):
    n = [dual(x) for x in normal]
    before, oldq, rest = ([dual(x) for x in v] for v in (before, oldq, rest))
    start = before
    end = [D(x, mp.eye(3)[i, :]) for i, x in enumerate(target)]
    branches, boundaries = [], []
    if all(x.v == y.v for x, y in zip(start, end)):
        raise ValueError('exact unchanged-target shortcut: derivative is not inferred')

    def closed(old, qold, after):
        delta = sub(after, old)
        tauold = sub(qold, scale(n, dot(qold, n)))
        dt = sub(delta, scale(n, dot(delta, n)))
        trial = add(tauold, scale(dt, kt))
        gap = dot(sub(after, rest), n)
        normalq = gap * kn if gap.v <= 0 else D(0)
        limit = -normalq * mu
        length2 = dot(trial, trial)
        length = length2.sqrt() if length2.v else D(0)
        if gap.v == 0:
            boundaries.append('gap=0')
        if length.v == limit.v:
            boundaries.append('friction=cap')
        if length.v > limit.v and length.v > 0:
            tau, branch = scale(trial, limit / length), 'slip'
        else:
            tau, branch = trial, 'stick'
        return add(scale(n, normalq), tau), branch

    for k in range(1, count + 1):
        after = end if k == count else add(start, scale(sub(end, start), mp.mpf(k) / count))
        g0, g1 = dot(sub(before, rest), n), dot(sub(after, rest), n)
        if g0.v <= 0 and g1.v <= 0:
            q, branch = closed(before, oldq, after)
        elif g0.v > 0 and g1.v > 0:
            q, branch = [D(0) for _ in range(3)], 'open'
        else:
            crossing = add(before, scale(sub(after, before), g0 / (g0 - g1)))
            crossing = sub(crossing, scale(n, dot(sub(crossing, rest), n)))
            if g0.v <= 0:
                q, branch = [D(0) for _ in range(3)], 'open'
            else:
                q, branch = closed(crossing, [D(0) for _ in range(3)], after)
        branches.append(branch)
        before, oldq = after, q
    return (mp.matrix([x.v for x in q]),
            mp.matrix([[x.d[j] for j in range(3)] for x in q]),
            dict(branch=branches[-1], path=branches, boundaries=boundaries))
