"""Independent principal-stress return for the uncapped nonassociated C-03 material."""

import argparse
import gzip
import json
from pathlib import Path

import numpy as np


POINTS = (108, 143, 152, 182, 185, 273)
LABELS = ("elastic", "shear_face", "compression_edge", "extension_edge", "shear_apex")
LEVELS = (32, 64, 128)
WORK_LIMIT = 0.01
REFERENCE_LIMIT = 0.001
EYE = np.eye(3)


def tensor(voigt, engineering=False):
    x, y, z, xy, yz, zx = np.moveaxis(np.asarray(voigt), -1, 0)
    if engineering:
        xy, yz, zx = xy / 2, yz / 2, zx / 2
    return np.stack([x, xy, zx, xy, y, yz, zx, yz, z], axis=-1).reshape(
        x.shape + (3, 3)
    )


def work(stress, dep):
    """Direct contraction, then the independently evaluated pressure and q terms."""
    p = np.trace(stress, axis1=-2, axis2=-1) / 3
    ev = np.trace(dep, axis1=-2, axis2=-1)
    dev_s = stress - p[..., None, None] * EYE
    dev_e = dep - (ev / 3)[..., None, None] * EYE
    shear = np.linalg.norm(dev_s, axis=(-1, -2)) * np.linalg.norm(dev_e, axis=(-1, -2))
    return np.stack([(stress * dep).sum(axis=(-1, -2)), p * ev, shear], axis=-1)


class PointReturn:
    """One point per row. Pool adjacent ordered principal stresses at edges."""

    def __init__(self, material):
        if (
            material["flow"] != "non_associated"
            or material["tension_cutoff_kpa"] is not None
        ):
            raise ValueError("requires uncapped nonassociated Mohr-Coulomb")
        self.c = material["cohesion_kpa"]
        young, nu = material["young_kpa"], material["poisson"]
        self.g2 = young / (1 + nu)
        self.lame = self.g2 * nu / (1 - 2 * nu)
        phi, psi = np.radians([material["friction_deg"], material["dilation_deg"]])
        self.k = 2 * self.c * np.cos(phi)
        self.normal = np.array([-1 - np.sin(phi), 0.0, 1 - np.sin(phi)])
        flow = np.array([-1 - np.sin(psi), 0.0, 1 - np.sin(psi)])
        self.shift = self.g2 * flow + self.lame * flow.sum()
        self.projections = []
        for groups in (
            ((0,), (1,), (2,)),
            ((0, 1), (2,)),
            ((0,), (1, 2)),
            ((0, 1, 2),),
        ):
            P = np.zeros((3, 3))
            for group in groups:
                P[np.ix_(group, group)] = 1 / len(group)
            self.projections.append((P, self.normal @ P @ self.shift))

    def step(self, old, strain):
        trial = (
            old
            + self.g2 * strain
            + self.lame * np.trace(strain, axis1=1, axis2=2)[:, None, None] * EYE
        )
        eig, axes = np.linalg.eigh(trial)
        width = 1e-10 * np.maximum.reduce(
            [np.abs(eig).max(axis=1), np.full(len(eig), self.c), np.ones(len(eig))]
        )
        pending = eig @ self.normal - self.k > width
        new = trial.copy()
        dep = np.zeros_like(trial)
        labels = np.zeros(len(old), dtype=int)
        # For fixed multiplier the sorted return is the isotonic projection of
        # trial - multiplier*C*dg. Each pooling has a linear consistency equation.
        for P, slope in self.projections:
            ids = np.flatnonzero(pending)
            if len(ids) == 0:
                break
            if slope <= 0:
                continue
            e = eig[ids]
            multiplier = (e @ P @ self.normal - self.k) / slope
            raw = e - multiplier[:, None] * self.shift
            corrected = raw @ P
            tol = width[ids, None]
            gaps = np.diff(corrected, axis=1)
            dual = np.cumsum(raw - corrected, axis=1)[:, :2]
            valid = (
                (multiplier >= 0)
                & (gaps >= -tol).all(axis=1)
                & (dual >= -tol).all(axis=1)
                & (np.abs(dual * gaps) <= tol * np.maximum(1.0, np.abs(gaps))).all(
                    axis=1
                )
            )
            take = ids[valid]
            if not len(take):
                continue
            final = corrected[valid]
            delta = eig[take] - final
            plastic = (
                delta
                - self.lame
                * delta.sum(axis=1, keepdims=True)
                / (self.g2 + 3 * self.lame)
            ) / self.g2
            basis = axes[take]
            new[take] = (basis * final[:, None, :]) @ basis.transpose(0, 2, 1)
            dep[take] = (basis * plastic[:, None, :]) @ basis.transpose(0, 2, 1)
            ties = np.abs(gaps[valid]) <= width[take, None]
            labels[take] = 1 + ties[:, 0].astype(int) + 2 * ties[:, 1].astype(int)
            pending[take] = False
        if np.any(pending):
            raise ValueError("no admissible ordered principal return")
        return new, dep, labels


def integrate(data, count):
    law = PointReturn(data["material"])
    stress = tensor(data["initial_stress_kpa"])
    cumulative = np.zeros((len(stress), 3))
    history, labels = [], []
    for row in data["frames"]:
        de = tensor(row["strain"], engineering=True) / count
        for _ in range(count):
            stress, plastic, code = law.step(stress, de)
            cumulative += work(stress, plastic)
        history.append(cumulative.copy())
        labels.append(code)
    return np.array(history), np.array(labels)


def discrepancies(candidate, reference):
    scale = np.max(np.abs(reference), axis=0)
    error = np.max(np.abs(candidate - reference), axis=0)
    # An identically zero reference work requires exact zero for that component.
    return np.divide(
        error, scale, out=np.where(error == 0.0, 0.0, np.inf), where=scale != 0.0
    )


def compare(data):
    frames = data["frames"]
    if (
        data["points"] != list(POINTS)
        or len(frames) != 224
        or [(r["stage"], r["progress"]) for r in frames]
        != [(s, i / 32) for s in range(7) for i in range(1, 33)]
    ):
        raise ValueError("incomplete C-03 point history")
    levels = [integrate(data, count) for count in LEVELS]
    counts = list(LEVELS)
    for count in (256, 512, 1024):
        previous, current = levels[-2:]
        if np.all(
            discrepancies(previous[0], current[0]) <= REFERENCE_LIMIT
        ) and np.array_equal(previous[1], current[1]):
            break
        counts.append(count)
        levels.append(integrate(data, count))
    fine, labels = levels[-1]
    fem = np.cumsum([r["fem_work_increment_kpa"] for r in frames], axis=0)
    fem_labels = np.array(
        [[LABELS.index(label) for label in r["fem_label"]] for r in frames]
    )
    refinement = [discrepancies(a[0], b[0]) for a, b in zip(levels, levels[1:])]
    err = discrepancies(fem, fine)
    points = []
    for j, gp in enumerate(POINTS):
        mismatch = np.flatnonzero(labels[:, j] != fem_labels[:, j])
        ref_labels_ok = np.array_equal(levels[-2][1][:, j], labels[:, j])
        points.append(
            dict(
                gp=gp,
                reference_levels=counts,
                matched=224 - len(mismatch),
                reference_labels_stable=ref_labels_ok,
                classification_passed=ref_labels_ok and not len(mismatch),
                mismatches=[
                    dict(
                        stage=frames[i]["stage"],
                        progress=frames[i]["progress"],
                        fem=LABELS[fem_labels[i, j]],
                        reference=LABELS[labels[i, j]],
                    )
                    for i in mismatch
                ],
                fem_final_kpa=fem[-1, j].tolist(),
                reference_final_kpa=fine[-1, j].tolist(),
                relative_max=err[j].tolist(),
                refinement=[r[j].tolist() for r in refinement],
                work_passed=(
                    (err[j] <= WORK_LIMIT) & (refinement[-1][j] <= REFERENCE_LIMIT)
                ).tolist(),
            )
        )
    return points, dict(
        reference_running_work_kpa=fine.tolist(), reference_label=labels.tolist()
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ("inputs", "histories", "summary"):
        p.add_argument(arg, type=Path)
    a = p.parse_args()
    with gzip.open(a.inputs, "rt", encoding="utf8") as f:
        routes = json.load(f)
    if set(routes) != {"m0c1", "m1c0", "m1c1"}:
        raise ValueError("expected three native strain histories")
    summary, histories = {}, {}
    for route, data in routes.items():
        summary[route], histories[route] = compare(data)
    with gzip.GzipFile(str(a.histories), "wb", mtime=0) as out:
        out.write(
            json.dumps(histories, separators=(",", ":"), allow_nan=False).encode()
        )
    a.summary.write_text(
        json.dumps(
            dict(
                levels=LEVELS,
                work_limit=WORK_LIMIT,
                reference_limit=REFERENCE_LIMIT,
                labels=LABELS,
                routes=summary,
            ),
            indent=2,
            allow_nan=False,
        )
        + "\n",
        encoding="utf8",
    )
