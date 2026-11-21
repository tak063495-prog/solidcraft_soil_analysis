"""Deterministic observation-only transfer; no constitutive state is updated."""
import math
import numpy as np


def dot(a, b):
    return math.fsum(float(x) * float(y) for x, y in zip(a, b, strict=True))


def solve(a, b):
    """Small dense solve with fixed pivot tie-breaking and sum order."""
    a = [list(map(float, row)) for row in a]
    b = [list(map(float, row)) for row in b]
    n = len(a)
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: (abs(a[i][k]), -i))
        if abs(a[pivot][k]) < 1e-15:
            raise ValueError('rank deficient projection')
        a[k], a[pivot] = a[pivot], a[k]
        b[k], b[pivot] = b[pivot], b[k]
        for i in range(k + 1, n):
            factor = a[i][k] / a[k][k]
            for j in range(k + 1, n):
                a[i][j] -= factor * a[k][j]
            for j in range(len(b[i])):
                b[i][j] -= factor * b[k][j]
    result = [[0.] * len(b[0]) for _ in range(n)]
    for i in reversed(range(n)):
        for j in range(len(b[0])):
            result[i][j] = (b[i][j] - math.fsum(a[i][k]*result[k][j] for k in range(i+1,n)))/a[i][i]
    return np.asarray(result)


def p2(r):
    x, y = r
    return [1., x, y, x*x, x*y, y*y]


def transfer(rule_rs, rule_weights, targets_rs, method):
    """Return target x GP observation matrix on a single parent/child face."""
    rs = np.asarray(rule_rs, dtype=float)
    targets = np.asarray(targets_rs, dtype=float)
    if method == 'nearest':
        return np.asarray([[float(j == min(range(len(rs)), key=lambda k: (dot(rs[k]-r, rs[k]-r), k))) for j in range(len(rs))] for r in targets])
    if method == 'p2_l2':
        basis = [p2(r) for r in rs]
        gram = [[math.fsum(w*b[i]*b[j] for w,b in zip(rule_weights,basis,strict=True)) for j in range(6)] for i in range(6)]
        rhs = [[w*b[i] for w,b in zip(rule_weights,basis,strict=True)] for i in range(6)]
        coefficients = solve(gram,rhs)
        return np.asarray([[dot(p2(r),column) for column in coefficients.T] for r in targets])
    if method != 'cubic':
        raise ValueError(method)
    uv = [(float(x), float(y)/(1-float(x))) for x,y in rs]
    def clustered(values):
        nodes=[]
        for v in sorted(values):
            if not nodes or abs(v-nodes[-1])>1e-12:
                nodes.append(v)
        return nodes
    # Geometry roundoff only: the prescribed tensor abscissae must form 4x4.
    us=clustered(x for x,y in uv)
    vs=clustered(y for x,y in uv)
    if len(us)*len(vs) != len(rs):
        raise ValueError('cubic requires a tensor collapsed rule')
    def lagrange(nodes, x, j):
        return math.prod((x-t)/(nodes[j]-t) for k,t in enumerate(nodes) if j != k)
    output=[]
    for r in targets:
        x,y=map(float,r); v=y/(1-x) if x != 1 else 0.
        output.append([lagrange(us,x,min(range(len(us)),key=lambda j:abs(us[j]-u)))*lagrange(vs,v,min(range(len(vs)),key=lambda j:abs(vs[j]-s))) for u,s in uv])
    return np.asarray(output)


def apply(matrix, values):
    values=np.asarray(values,dtype=float)
    if not np.isfinite(values).all():
        raise ValueError('nonfinite history')
    shape=values.shape
    flat=values.reshape((shape[0],-1))
    matrix=np.asarray(matrix,dtype=float)
    result=np.zeros((len(matrix),flat.shape[1]));correction=np.zeros_like(result)
    # GP index order is fixed; elementwise operations avoid BLAS reduction order.
    for i in range(len(flat)):
        term=matrix[:,i,None]*flat[i][None,:]
        updated=result+term
        correction+=np.where(abs(result)>=abs(term),(result-updated)+term,(term-updated)+result)
        result=updated
    return (result+correction).reshape((len(matrix),)+shape[1:])


def norm(values, weights=None):
    a=np.asarray(values,dtype=float)
    if not np.isfinite(a).all():
        raise ValueError('nonfinite norm')
    if weights is None:
        return math.sqrt(math.fsum(float(x)*float(x) for x in a.ravel(order='C')))
    w=np.broadcast_to(np.asarray(weights,dtype=float),a.shape)
    if not np.isfinite(w).all() or (w<0).any():
        raise ValueError('invalid norm weights')
    return math.sqrt(math.fsum(float(t)*float(x)*float(x) for t,x in zip(w.ravel(order='C'),a.ravel(order='C'),strict=True)))
