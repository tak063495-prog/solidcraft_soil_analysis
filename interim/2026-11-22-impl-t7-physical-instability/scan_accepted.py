"""Search the accepted histories, including restart boundaries, for simple branch points."""

import argparse
from pathlib import Path
import gzip
import json

import numpy as np
from scipy.optimize import brentq, root

from continuation import TrialRejected, bifurcation, serial, spectrum


def branch_in_interval(op, parent, za, zb, event, tol, left, right):
    """Follow the locally single-valued parent and locate a zero of its bordered test."""
    a, b = map(np.asarray, (za, zb))
    lo, hi = sorted((a[-1], b[-1]))
    if hi == lo:
        return dict(kind="D", reason="no monotone load interval")

    def parent_at(load):
        if not lo <= load <= hi:
            raise TrialRejected("outside the accepted parent interval")
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
        load = brentq(bordered_test, lo, hi, xtol=np.finfo(float).eps)
        z = np.r_[parent_at(load), load]
        # The span is a predictor scale. A successful result must equilibrate
        # on both sides at three smaller offsets with the original data.
        eps = np.linalg.norm(b - a) / np.array([8.0, 16.0, 32.0])
        result = bifurcation(op, parent, z, eps, parent_at, tol)
        return dict(z=z, offsets=eps, **result)
    except (TrialRejected, np.linalg.LinAlgError, ValueError) as exc:
        return dict(kind="D", reason=str(exc))


def scan(op, points, tol, prepare):
    rows, candidates = [], []
    previous = None
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
        if previous is not None and previous["event"] == event:
            pl, pr = previous["left"], previous["right"]
            row["weak_alignment"] = [float(pl @ left), float(pr @ right)]
            # A sign change in the *oriented pair* or an observed simple
            # null triggers localization; neither is a physical verdict.
            crossed = (pl @ left) * (pr @ right) < 0
            if crossed or s["nullity"] == 1 or previous["nullity"] == 1:
                result = branch_in_interval(
                    op, state, prior["z"], z, event, tol, pl, pr
                )
                candidates.append(dict(accepted_indices=[i - 1, i], **result))
        rows.append(row)
        previous = dict(event=event, left=left, right=right, nullity=s["nullity"])
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
