"""Independent C-01 predictors, active-branch Jacobians and spherical continuation."""

import bisect
import json
import math

import numpy as np

STAGES = ("geostatic", "shear", "open", "recontact", "reverse", "reload", "return")


class Seed:
    def __init__(self, path, nodes):
        self.times = [[0.0] for _ in STAGES]
        self.fields = [[] for _ in STAGES]
        self.fields[0].append(np.zeros(nodes * 3))
        with path.open(encoding="utf-8") as stream:
            for line in stream:
                stage, time, u = json.loads(line)
                field = np.asarray(u, dtype=float)
                if not 0 <= stage < 7 or field.shape != (nodes, 3) or not np.all(np.isfinite(field)):
                    raise ValueError("invalid source displacement")
                if not self.times[stage][-1] < time <= 1:
                    raise ValueError("unordered source time")
                if not self.fields[stage]:
                    if self.times[stage - 1][-1] != 1:
                        raise ValueError("missing previous stage endpoint")
                    self.fields[stage].append(self.fields[stage - 1][-1])
                self.times[stage].append(time)
                self.fields[stage].append(field.ravel())
        if any(t[-1] != 1 for t in self.times):
            raise ValueError("incomplete source path")

    def interpolate(self, name, progress):
        s = STAGES.index(name)
        times = self.times[s]
        fields = self.fields[s]
        i = bisect.bisect_left(times, progress)
        if times[i] == progress:
            return fields[i]
        frac = (progress - times[i - 1]) / (times[i] - times[i - 1])
        return fields[i - 1] + frac * (fields[i] - fields[i - 1])


class Sphere:
    def __init__(self, scale):
        self.scale = scale
        self.pending = None
        self.calls = 0
        self.corrections = 0

    def retry(self, name, start, end, base, old, new, path_rejected):
        if name == "geostatic" or path_rejected:
            self.pending = None
            return None
        self.calls += 1
        length = math.hypot(
            end - start, float(np.linalg.norm(new - old)) / (math.sqrt(new.size) * self.scale)
        )
        if self.pending is None or self.pending[0] != (name, start, end):
            if length <= base:
                return None
            lo, low, hi, high = start, 0.0, end, length
        else:
            _, lo, low, hi, high = self.pending
            if length == base:
                self.pending = None
                return None
            if length < base:
                lo, low = end, length
            else:
                hi, high = end, length
        candidate = lo + (hi - lo) * (base - low) / (high - low)
        if length <= base and (candidate <= end or math.nextafter(lo, hi) >= hi):
            self.pending = None
            return None
        if math.nextafter(lo, hi) >= hi:
            candidate = lo
        else:
            candidate = min(math.nextafter(hi, lo), max(math.nextafter(lo, hi), candidate))
        self.pending = ((name, start, candidate), lo, low, hi, high)
        self.corrections += 1
        return candidate - start


def active_tangent(system, old, traction, target, counts):
    """Differentiate the independent disk projection at each local partition."""
    n = system.normal
    normal = np.einsum("qi,qj->qij", n, n)
    plane = np.eye(3) - normal
    result = np.zeros_like(plane)
    for count in np.unique(counts):
        mask = counts == count
        nn = n[mask]
        P = plane[mask]
        kn = system.kn[mask]
        kt = system.kt[mask]
        mu = system.mu[mask]
        delta = (target[mask] - old[mask]) / count
        dg = (delta * nn).sum(axis=1)
        dt = delta - dg[:, None] * nn
        tau = np.einsum("qij,qj->qi", P, traction[mask])
        jac = np.zeros_like(P)
        gap = (old[mask] * nn).sum(axis=1)
        for j in range(1, int(count) + 1):
            next_gap = gap + dg
            closed = next_gap <= 0
            recontact = (gap > 0) & (next_gap == 0)
            crossing = (gap > 0) & (next_gap < 0)
            tau[gap > 0] = 0
            jac[gap > 0] = 0
            frac = np.ones(len(gap))
            df = np.zeros_like(nn)
            frac[recontact] = 0
            frac[crossing] = next_gap[crossing] / dg[crossing]
            df[crossing] = (
                ((j * dg[crossing] - next_gap[crossing]) / (count * dg[crossing] ** 2))[:, None]
                * nn[crossing]
            )
            trial = tau + kt[:, None] * dt * frac[:, None]
            dj = jac + kt[:, None, None] * (
                frac[:, None, None] * P / count + dt[:, :, None] * df[:, None, :]
            )
            cap = mu * kn * np.maximum(-next_gap, 0)
            norm = np.linalg.norm(trial, axis=1)
            sliding = closed & (norm > cap)
            tau = trial.copy()
            jac = dj.copy()
            direction = trial[sliding] / norm[sliding, None]
            projector = np.eye(3) - direction[:, :, None] * direction[:, None, :]
            jac[sliding] = (cap[sliding] / norm[sliding])[:, None, None] * (projector @ dj[sliding]) - (
                mu[sliding] * kn[sliding] * j / count
            )[:, None, None] * direction[:, :, None] * nn[sliding, None, :]
            tau[sliding] = cap[sliding, None] * direction
            if recontact.any():
                # Take the closed-side derivative of the trial shear capped by Coulomb.
                increment = kt[recontact, None] * dt[recontact]
                length = np.linalg.norm(increment, axis=1)
                cap_slope = mu[recontact] * kn[recontact]
                slips = length > cap_slope * gap[recontact]
                divisor = np.where(slips, length, gap[recontact])
                slope = (increment / divisor[:, None]) * np.where(slips, cap_slope, 1)[:, None]
                jac[recontact] = -(j / count) * slope[:, :, None] * nn[recontact, None, :]
            tau[~closed] = 0
            jac[~closed] = 0
            gap = next_gap
        jac[mu == 0] = 0
        result[mask] = jac + (kn * (gap <= 0))[:, None, None] * normal[mask]
    return result


class Methods:
    def __init__(self, sink, seed, active, arc, history=None):
        self.history = history
        self.sink = sink
        self.seed = seed
        self.active = active
        self.arc = arc
        self.guesses = 0
        self.tangents = 0

    def tangent(self, system, old, traction, target, counts):
        self.tangents += 1
        return active_tangent(system, old, traction, target, counts)

    def initial(self, name, time, old):
        if self.seed is None:
            return old
        self.guesses += 1
        return self.seed.interpolate(name, time)

    def commit(self, name, progress, u, dissipation):
        if self.history is not None:
            print(json.dumps([STAGES.index(name), progress, dissipation.tolist()], separators=(",", ":")), file=self.history)
        if self.sink:
            print(
                json.dumps(
                    [STAGES.index(name), progress, u.reshape(-1, 3).tolist()], separators=(",", ":")
                ),
                file=self.sink,
            )
