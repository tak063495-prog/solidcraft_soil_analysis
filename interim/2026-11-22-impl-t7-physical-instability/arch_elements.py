"""Planar corotational Euler-Bernoulli beam discretization of the published deep arch."""

from pathlib import Path
import argparse
import json

import numpy as np
from continuation import Controls, Evaluation, serial, simple_point, trace
from references import A, EI, R, force_unit, span, theta0


class Arch:
    def __init__(self, n):
        # Keep a node at the concentrated load even on the odd grid.
        self.apex = n // 2
        t = np.r_[
            np.linspace(0, 0.5, self.apex + 1),
            np.linspace(0.5, 1, n - self.apex + 1)[1:],
        ]
        angle = theta0 + span * t
        xy = np.stack((np.cos(angle), np.sin(angle)), axis=1)
        self.d0 = np.diff(xy, axis=0)
        self.h = np.linalg.norm(self.d0, axis=1)
        self.free = np.setdiff1d(np.arange(3 * (n + 1)), [0, 1, 2, 3 * n, 3 * n + 1])
        self.load = np.zeros(3 * (n + 1))
        self.load[3 * self.apex + 1] = -1
        self.ids = 3 * np.arange(n)[:, None] + np.arange(6)
        self.L = np.array([[-1.0, 0, 0, 1, 0, 0], [0, -1.0, 0, 0, 1, 0]])
        self.D = np.zeros((n, 3, 3))
        self.D[:, 0, 0] = 1e6 * A * R**2 / EI / self.h
        self.D[:, 1:, 1:] = np.array([[4.0, 2.0], [2.0, 4.0]]) / self.h[:, None, None]

    def __call__(self, parent, z):
        u = np.zeros_like(self.load)
        u[self.free] = z[:-1]
        e = u[self.ids]
        d = self.d0 + e @ self.L.T
        ell = np.linalg.norm(d, axis=1)
        t = d / ell[:, None]
        n = np.stack((-t[:, 1], t[:, 0]), axis=1)
        rot = np.arctan2(
            self.d0[:, 0] * d[:, 1] - self.d0[:, 1] * d[:, 0],
            np.sum(self.d0 * d, axis=1),
        )
        basic = np.column_stack((ell - self.h, e[:, 2] - rot, e[:, 5] - rot))
        f = (self.D @ basic[:, :, None])[:, :, 0]
        g = t @ self.L
        h = n @ self.L
        b = np.stack((g, -h / ell[:, None], -h / ell[:, None]), axis=1)
        b[:, 1, 2] += 1
        b[:, 2, 5] += 1
        ke = b.transpose(0, 2, 1) @ self.D @ b
        ke += np.einsum("e,ei,ej->eij", f[:, 0] / ell, h, h)
        m = (f[:, 1] + f[:, 2]) / ell**2
        ke += np.einsum("e,ei,ej->eij", m, h, g) + np.einsum("e,ei,ej->eij", m, g, h)
        internal = np.zeros_like(u)
        k = np.zeros((len(u), len(u)))
        np.add.at(internal, self.ids, np.einsum("eai,ea->ei", b, f))
        np.add.at(k, (self.ids[:, :, None], self.ids[:, None, :]), ke)
        q = -self.load[self.free]
        res = internal[self.free] + q * z[-1]
        return Evaluation(
            res,
            k[np.ix_(self.free, self.free)],
            q,
            {},
            max(1.0, abs(z[-1])),
            details=dict(crown=-u[3 * self.apex + 1] * R),
        )


def study(n):
    model = Arch(n)
    z = np.zeros(len(model.free) + 1)
    forward = z.copy()
    forward[-1] = 1

    def stop(points):
        loads = [p["z"][-1] for p in points]
        return loads[-1] < max(loads) - 0.45

    path = trace(
        model, {}, z, forward, Controls(0.18, 1500, 1e-9, 1e-10, 20), stop=stop
    )
    points = [p["z"] for p in path["accepted"]]
    curve = [
        dict(load=z[-1] * force_unit, downward_crown=model({}, z).details["crown"])
        for z in points
    ]
    peak = int(np.argmax([c["load"] for c in curve]))
    fold = (
        simple_point(model, {}, points[peak - 1], points[peak + 1], 1e-9)
        if 0 < peak < len(points) - 1
        else dict(kind="D")
    )
    return dict(
        elements=n,
        termination=path["termination_reason"],
        trials=path["trials"],
        curve=curve,
        peak_sample=curve[peak],
        fold=fold,
        fold_load=fold["z"][-1] * force_unit if fold["kind"] == "simple fold" else None,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    for n in (15, 30, 60):
        with (args.output / f"beam-{n}.json").open("x") as f:
            json.dump(study(n), f, default=serial, allow_nan=False)
