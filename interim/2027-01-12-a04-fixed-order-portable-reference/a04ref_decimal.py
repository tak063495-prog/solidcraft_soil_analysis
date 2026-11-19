"""Input-only Decimal RT0 state and flux arithmetic for the small A-04 soil fixture."""
from decimal import Decimal, localcontext
from functools import lru_cache

import numpy as np


PRECISION = 40
D = Decimal
ZERO = D(0)
ONE = D(1)


def dec(x):
    return D.from_float(float(x))


def norm(x):
    return sum((v*v for v in x), ZERO).sqrt()


@lru_cache(maxsize=32768)
def retention(psi, alpha, n, delta_theta, ss):
    if psi >= 0:
        return ss*psi, ONE
    m = ONE-ONE/n
    t = (alpha*(-psi))**n
    se = (ONE+t)**(-m)
    # 1 - Se^(1/m) = t/(1+t), including the saturation limit.
    kr = se.sqrt()*(ONE-(t/(ONE+t))**m)**2
    if kr <= 0:
        raise ArithmeticError('invalid retention; no conductivity floor')
    return delta_theta*(se-ONE), kr


def inverse_mass(a):
    """Four unpivoted positive-pivot eliminations, as in the original local solve."""
    n = len(a)
    b = [[D(i == j) for j in range(n)] for i in range(n)]
    for k in range(n):
        if a[k][k] <= 0:
            raise ArithmeticError('RT0 nonpositive pivot')
        for i in range(k+1, n):
            f = a[i][k]/a[k][k]
            for j in range(k+1, n):
                a[i][j] -= f*a[k][j]
            for j in range(n):
                b[i][j] -= f*b[k][j]
    for k in reversed(range(n)):
        for j in range(n):
            b[k][j] = (b[k][j]-sum(a[k][i]*b[i][j] for i in range(k+1,n)))/a[k][k]
    return np.array(b, dtype=object)


class Hydraulic:
    def __init__(self, geometry, ledger, log):
        self.g = geometry
        self.ledger = ledger
        self.log = log
        if geometry.datum != 0 or geometry.nc > 12:
            raise ValueError('limited to zero datum and at most 12 affine cells')
        self.laws = []
        self.weights = []
        self.z = []
        self.grams = []
        for e, material in enumerate(geometry.materials):
            law = material['retention']
            self.laws.append((dec(law['alpha_per_m']), dec(law['n']),
                dec(law['saturated_water_content'])-dec(law['residual_water_content']),
                dec(law['saturated_specific_storage_per_m'])))
            self.z.append([dec(v) for v in geometry.position[e,:,2]])
            w = [dec(v) for v in geometry.weight[e]]
            self.weights.append(w)
            # Same binary64 input geometry and Jacobi points as the host reference.
            basis = geometry.basis[e]
            gram = (basis @ geometry.inverse_k[e]) @ basis.transpose(0,2,1)
            self.grams.append(np.array([[[dec(x)*w[q] for x in row] for row in block]
                for q,block in enumerate(gram)], dtype=object).reshape(len(w),16).T)
        self.cell = lru_cache(maxsize=512)(self._cell)

    def _cell(self, e, h):
        props = [retention(h-z,*self.laws[e]) for z in self.z[e]]
        storage = sum(w*p[0] for w,p in zip(self.weights[e],props))
        inv_kr = np.array([ONE/p[1] for p in props], dtype=object)
        mass = (self.grams[e] @ inv_kr).reshape(4,4)
        return storage, inverse_mass(mass.tolist())

    def materials(self, h):
        values = [self.cell(e,v) for e,v in enumerate(h)]
        return np.array([v[0] for v in values]), np.stack([v[1] for v in values])

    def commit(self,h,trace,storage):
        for e,hv,sv in zip(self.g.elements,h,storage):
            self.ledger['cell'][int(e)] = (hv,sv)
        self.ledger['face'].update(zip(self.g.keys,trace))

    def initial(self, flow):
        g = self.g
        h = np.array([dec(x) for x in g.seed])
        a = flow['initial_head']
        value = a['head_m']+(g.face_centres-a.get('reference_point_m',[0.,0.,0.])) @ np.array(a.get('gradient',[0.,0.,0.]))
        trace = np.array([dec(v) for v in value])
        storage,inv = self.materials(h)
        flux = (inv @ (h[:,None]-trace[g.cf])[:,:,None])[:,:,0]
        self.commit(h,trace,storage)
        return {k:np.asarray(v,float) for k,v in dict(h=h,trace=trace,storage=storage,flux=flux).items()}

    def step(self, h_view, trace_view, storage_view, period, dt):
        g = self.g
        h,old = map(np.array, zip(*(self.ledger['cell'][int(e)] for e in g.elements)))
        trace = np.array([self.ledger['face'][k] for k in g.keys])
        # The inherited excavation controller selects cells/faces by identity.
        for raw,view in ((h,h_view),(trace,trace_view),(old,storage_view)):
            np.testing.assert_array_equal(raw.astype(float),view)
        heads, prescribed = g.boundaries(period)
        for f,v in heads.items():
            trace[f] = dec(v)
        free = [f for f in range(len(trace)) if f not in heads]
        source = np.array([dec(v) for v in g.source])
        boundary = np.array([dec(prescribed.get(f,0.)) for f in free])
        x = np.concatenate((h,trace[free]))
        interval = dec(dt)

        def evaluate(candidate):
            hh = candidate[:g.nc]
            ll = trace.copy()
            ll[free] = candidate[g.nc:]
            storage,inv = self.materials(hh)
            q = (inv @ (hh[:,None]-ll[g.cf])[:,:,None])[:,:,0]
            rate = (storage-old)/interval
            balance = q.sum(axis=1)+rate-source
            face = np.array([sum(q[e,j] for e,j in g.owners[f]) for f in free])-boundary
            cell_scale = abs(q).sum(axis=1)+abs(rate)+abs(source)
            face_scale = np.array([sum(abs(q[e,j]) for e,j in g.owners[f]) for f in free])+abs(boundary)
            return np.r_[balance,face], np.r_[cell_scale,face_scale], ll,storage,q

        max_linear = 0.
        for iteration in range(60):
            r,scale,ll,storage,q = evaluate(x)
            merit = max(abs(r)/(dec(1e-16)+dec(1e-8)*scale))
            old_norm = norm(r)
            row = dict(period=period['name'],dt=dt,iteration=iteration,merit=str(merit),
                       norm=str(old_norm),trials=[])
            self.log.append(row)
            if merit < dec(.01):
                break
            jac = np.empty((len(x),len(x)),dtype=float)
            for col in range(len(x)):
                delta = dec(1e-6)*(ONE+abs(x[col]))
                xp,xm = x.copy(),x.copy()
                xp[col] += delta
                xm[col] -= delta
                jac[:,col] = (evaluate(xp)[0]-evaluate(xm)[0])/(2*delta)
            dx = np.linalg.solve(jac,-r.astype(float))
            error = [dec(v)+rr for v,rr in zip(jac@dx,r)]
            linear = float(norm(error)/old_norm)
            row['linear'] = linear
            if linear > 1e-10:
                raise ArithmeticError(('true hydraulic linear residual',linear))
            max_linear = max(max_linear,linear)
            direction = np.array([dec(v) for v in dx])
            for cut in range(20):
                trial = x+direction/(2**cut)
                measure = norm(evaluate(trial)[0])
                row['trials'].append(dict(cut=cut,norm=str(measure)))
                if measure < old_norm:
                    x = trial
                    break
            else:
                raise ArithmeticError(('hydraulic line search',str(merit)))
        else:
            raise ArithmeticError('hydraulic Newton budget')
        h = x[:g.nc]
        self.commit(h,ll,storage)
        out = sum(q[e,j] for owners in g.owners if len(owners)==1 for e,j in owners)
        return dict(cell_head_m=h.astype(float),face_head_m=ll.astype(float),
                    outflow_m3_s=q.astype(float),storage_offset_m3=storage.astype(float),
                    boundary_outflow_m3_s=float(out),source_m3_s=float(sum(source)),
                    storage_rate_m3_s=float(sum(storage-old)/interval), iterations=iteration,
                    maximum_absolute_cell_residual_m3_s=float(max(abs(r[:g.nc]))),
                    maximum_scaled_residual=float(merit),maximum_true_linear_residual=max_linear)


def model_type(base,log):
    class Portable(base):
        def __init__(self,*args):
            super().__init__(*args)
            ledger = dict(cell={},face={})
            for model in (self.full,self.cut):
                solver = Hydraulic(model,ledger,log)
                model.step = solver.step
                if model is self.full:
                    self.initial = solver.initial

        def initial_hydraulic(self):
            return self.initial(self.flow)

        def solve(self,sink):
            with localcontext() as ctx:
                ctx.prec = PRECISION
                return super().solve(sink)
    return Portable
