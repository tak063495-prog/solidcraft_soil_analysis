"""Search the accepted histories, including restart boundaries, for simple branch points."""

import argparse
from pathlib import Path
import gzip
import json

import numpy as np
from scipy.optimize import brentq, root

from continuation import TrialRejected, bifurcation, serial, spectrum


def branch_in_interval(op, parent, accepted, event, tol, mode):
    """Recheck the locally single-valued parent and test for a simple branch point."""
    left, right = mode
    ordered = sorted(accepted, key=lambda z: z[-1])
    segments = list(zip(ordered, ordered[1:]))
    lo, hi = ordered[0][-1], ordered[-1][-1]
    if hi == lo:
        return dict(kind="D", reason="no monotone load interval")

    def parent_at(load):
        if not lo <= load <= hi:
            raise TrialRejected("outside the accepted parent interval")
        a, b = next((a, b) for a, b in segments if load <= b[-1])
        seed = a[:-1] + (load - a[-1]) / (b[-1] - a[-1]) * (b[:-1] - a[:-1])

        def ev(x):
            return op(parent, np.r_[x, load])

        fit = root(lambda x: ev(x).residual, seed, jac=lambda x: ev(x).k, tol=tol)
        e = ev(fit.x)
        if (
            not e.eligible
            or tuple(e.event) != event
            or np.linalg.norm(e.residual) > tol * e.denominator
        ):
            raise TrialRejected("parent equilibrium or active side unavailable")
        return fit.x

    def bordered_test(load):
        x = parent_at(load)
        e = op(parent, np.r_[x, load])
        bordered = np.block([[e.k, left[:, None]], [right[None, :], np.zeros((1, 1))]])
        return np.linalg.solve(bordered, np.r_[np.zeros(len(x)), 1.0])[-1]

    try:
        if len(ordered) == 3:
            # The middle equilibrium exhibited a simple zero. Recheck all
            # three points from the same checkpoint before seeding its modes.
            parent_at(lo)
            parent_at(hi)
            load = ordered[1][-1]
        else:
            load = brentq(bordered_test, lo, hi, xtol=np.finfo(float).eps)
        z = np.r_[parent_at(load), load]
        # Use the shorter leg when a saved critical point splits the interval.
        span = min(np.linalg.norm(b - a) for a, b in segments)
        eps = span / np.array([8.0, 16.0, 32.0])
        result = bifurcation(op, parent, z, eps, parent_at, tol)
        return dict(z=z, offsets=eps, **result)
    except (TrialRejected, np.linalg.LinAlgError, ValueError) as exc:
        return dict(kind="D", reason=str(exc))


def scan(op, points, tol, prepare):
    rows, observations, candidates = [], [], []
    for i, point in enumerate(points):
        z = np.asarray(point["z"])
        prior = points[max(i - 1, 0)]
        state = prepare(prior["state"])
        ev = op(state, z)
        if not ev.eligible or np.linalg.norm(ev.residual) > tol * ev.denominator:
            raise ValueError("saved equilibrium failed original acceptance")
        s = spectrum(ev.k, ev.q, ev.residual)
        # Fix the sign of each singular pair together by the first maximum
        # physical component. Subspaces are compared independently of sign.
        left, right = s["left"][:, -1], s["right"][:, -1]
        if right[np.argmax(abs(right))] < 0:
            left, right = -left, -right
        event = tuple(ev.event)
        row = dict(
            progress=z[-1],
            nullity=s["nullity"],
            sigma_min=s["singular_values"][-1],
            load_projection=s["null_load_norm"],
            residual_projection=s["null_residual_norm"],
            inertia=s.get("negative_eigen"),
            ldlt_inertia=s.get("negative_ldlt"),
            parent_index=max(0, i - 1),
        )
        row["physical_tangent_scale"] = "per scaled displacement"
        row["weak_load_projection"] = float(left @ ev.q)
        row["weak_residual_projection"] = float(left @ ev.residual)
        row["work_kn_m"] = (ev.details or {}).get("work_kn_m")
        if "eigenvalues" in s:
            row["eigenvalues"] = s["eigenvalues"]
        crossed = False
        if observations and observations[-1]["event"] == event:
            pl, pr = observations[-1]["mode"]
            row["weak_alignment"] = [float(pl @ left), float(pr @ right)]
            # A change in the oriented pair triggers localization.
            crossed = (pl @ left) * (pr @ right) < 0
        rows.append(row)
        observations.append(
            dict(z=z, event=event, mode=(left, right), nullity=s["nullity"], crossed=crossed)
        )

    i = 1
    while i < len(observations):
        start, stop = i - 1, i
        first, second = observations[start : stop + 1]
        if first["nullity"] == 0 and second["nullity"] == 1 and i + 1 < len(observations):
            third = observations[i + 1]
            a, b, c = (obs["z"][-1] for obs in (first, second, third))
            if third["nullity"] == 0 and (a < b < c or a > b > c):
                # Keep the saved zero inside one interval for the load probes.
                stop += 1
        interval = observations[start : stop + 1]
        if second["crossed"] or any(obs["nullity"] == 1 for obs in interval):
            if any(obs["event"] != first["event"] for obs in interval[1:]):
                result = dict(kind="D", reason="active-set change")
            else:
                result = branch_in_interval(
                    op,
                    prepare(points[start]["state"]),
                    [obs["z"] for obs in interval],
                    first["event"],
                    tol,
                    first["mode"],
                )
            candidates.append(
                dict(accepted_indices=list(range(start, stop + 1)), **result)
            )
        i = stop + 1
    return dict(accepted_points=rows, branch_candidates=candidates)


def study(folder, bridge):
    from compare import accepted, checkpoint
    from parallel import Model, Native, independent

    original = json.loads((folder / "run.json").read_bytes())
    model = Model(
        *(
            json.loads((folder / (k + ".json")).read_bytes())
            for k in ("mesh", "ground", "case")
        )
    )
    with gzip.open(folder / "independent-last.json.gz") as f:
        path = json.load(f)["stage_path"]
    path["holding"] = np.asarray(path["holding"])
    out = {}
    for series in json.loads((folder / "summary.json").read_bytes())["series"]:
        label = series["label"]
        engine, _, parts = label.split("-")
        points = accepted(folder / (label + ".json.gz"))
        segments = sorted(
            folder.glob(label + "-segment-*.json.gz"),
            key=lambda p: int(p.name.split("-segment-")[-1].split(".")[0]),
        )
        for segment in segments:
            batch = accepted(segment)
            if (
                batch[0]["z"] != points[-1]["z"]
                or batch[0]["state"] != points[-1]["state"]
            ):
                raise ValueError("restart does not join the last accepted point")
            points.extend(batch[1:])
        native = engine == "native"
        op = (
            Native(bridge, folder, model, original, int(parts))
            if native
            else independent(model, path, int(parts))
        )
        try:
            tol = (
                model.case["settings"]["relative_force_tolerance"] if native else 1e-10
            )
            out[label] = scan(op, points, tol, lambda state: checkpoint(state, engine))
        finally:
            if native:
                op.close()
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("bridge", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    data = {
        str(n): study(args.root / f"mesh-{n}", args.bridge.resolve()) for n in (18, 144)
    }
    args.output.write_text(
        json.dumps(data, default=serial, separators=(",", ":"), allow_nan=False) + "\n",
        encoding="utf-8",
    )
