"""Re-equilibrate selected coordinates within accepted brackets and compare physical fields."""

import argparse
import gzip
import json
from pathlib import Path

import ijson
import numpy as np
from scipy.linalg import subspace_angles

from continuation import Controls, direction, identity, serial, spectrum, trace
from parallel import Model, Native, U, independent, save


def accepted(path):
    with gzip.open(path) as f:
        return next(ijson.items(f, "accepted", use_float=True))


def checkpoint(state, engine):
    if engine == "independent":
        state["u"] = np.array(state["u"])
        state["joint"] = {k: np.array(v) for k, v in state["joint"].items()}
    return state


def at_progress(op, points, target, tol, engine):
    for a, b in zip(points, points[1:]):
        za, zb = np.array(a["z"]), np.array(b["z"])
        if za[-1] == target:
            return (
                za,
                checkpoint(a["state"], engine),
                dict(
                    reason="original accepted coordinate", parent=identity(a["state"])
                ),
            )
        if not za[-1] < target <= zb[-1]:
            continue
        if a["event"] != b["event"]:
            return None, None, dict(reason="active-set event bracket")
        parent = checkpoint(a["state"], engine)
        chord = zb - za
        try:
            tangent = direction(op(parent, za), chord)
        except ArithmeticError as error:
            return None, None, dict(reason=str(error))
        if tangent[-1] <= 0:
            return None, None, dict(reason="no monotone load control")
        plane = np.zeros(len(chord))
        plane[-1] = 1.0
        history = trace(
            op,
            parent,
            za,
            chord,
            Controls(
                (target - za[-1]) / tangent[-1], 1, tol, 1e-10, 18, event_bisections=0
            ),
            local_control=plane,
        )
        end = history["accepted"][-1]
        same = len(history["accepted"]) == 2 and tuple(end["event"]) == tuple(
            b["event"]
        )
        report = dict(
            reason=history["termination_reason"],
            parent=identity(parent),
            bracket=[za[-1], zb[-1]],
            matched_side=same,
            actual=end["z"][-1],
            trials=history["trials"],
        )
        if same:
            report["arrival"] = op(parent, end["z"])
            return end["z"], end["state"], report
        return None, None, report
    return None, None, dict(reason="no accepted bracket")


def fields(op, z, state, engine, model, arrival=None):
    ev = arrival if arrival is not None else op(state, z)
    sys = model.system
    if engine == "native":
        pts = [state["joints"]["points"][i] for i in op.order]
        u = np.array(state["displacement_m"]).ravel()
        q = np.array([p["traction_kpa"] for p in pts])
        jump = np.array([p["jump_m"] for p in pts])
        slip = np.array([p["plastic_slip_m"] for p in pts])
        friction = np.array([p["plastic_dissipation_kpa_m"] for p in pts])
        stress = np.array(
            [
                p["effective_stress_compression_positive_kpa"]
                for p in state["material_states"]
            ]
        )
        event = tuple(p["status"] for p in pts)
        generations = [v["generation"] for v in state["joints"]["lifecycle"]]
        basis = op.t
    else:
        joint = state["joint"]
        u = state["u"]
        q = joint["traction"]
        jump = joint["jump"]
        slip = joint["plastic"]
        friction = joint["pd"]
        stress = np.asarray(ev.candidate["fields"]["stress"])
        event = tuple(("open", "stick", "slip")[i] for i in joint["status"])
        generations = state["generation"]
        basis = model.spaces[False][0]
    diag = spectrum(ev.k / U, ev.q, ev.residual)
    surfaces = {}
    for seam in ("LowerSeam", "UpperSeam"):
        mask = np.array([p["joint_id"] == seam for p in sys.points])
        w = sys.weight[mask]
        surfaces[seam] = dict(
            area_m2=sum(w),
            normal_jump_m=w
            @ np.einsum("qa,qa->q", jump[mask], sys.normal[mask])
            / sum(w),
            traction_force_kn=w @ q[mask],
            friction_work_kn_m=w @ friction[mask],
        )
    return dict(
        z=z,
        stage="open",
        generations=generations,
        u=u,
        stress=stress,
        traction=q,
        slip=slip,
        jump=jump,
        surfaces=surfaces,
        event=event,
        relative_force=np.linalg.norm(ev.residual) / ev.denominator,
        spectral=diag,
        tangent_evaluation="arrival" if arrival else "checkpoint replay",
        right=basis @ diag["right"],
        left=basis @ diag["left"],
    )


def difference(a, b):
    out = {}
    for name in ("u", "stress", "traction", "slip", "jump"):
        x, y = map(np.asarray, (a[name], b[name]))
        absolute = np.linalg.norm(x - y)
        den = np.linalg.norm(y)
        out[name] = dict(
            absolute_l2=absolute,
            reference_l2=den,
            relative_l2=absolute / den if den else None,
        )
    for side in ("left", "right"):
        out[side + "_principal_angles_rad"] = subspace_angles(a[side], b[side])
    return out


def compare(folder, bridge):
    m = Model(
        *[
            json.loads((folder / (k + ".json")).read_bytes())
            for k in ("mesh", "ground", "case")
        ]
    )
    original = json.loads((folder / "run.json").read_bytes())
    with gzip.open(folder / "independent-last.json.gz") as f:
        path = json.load(f)["stage_path"]
    path["holding"] = np.asarray(path["holding"])
    report = {}
    for engine in ("native", "independent"):
        if engine == "independent" and path["name"] != "open":
            report[engine] = dict(reason="different stage and generation history")
            continue
        label = engine + "-0.02-1"
        files = [
            folder / (label + ".json.gz"),
            *sorted(
                folder.glob(label + "-segment-*.json.gz"),
                key=lambda p: int(p.name.split("-segment-")[-1].split(".")[0]),
            ),
        ]
        points = accepted(files[0])
        for f in files[1:]:
            points.extend(accepted(f)[1:])
        op = (
            Native(bridge, folder, m, original, 1)
            if engine == "native"
            else independent(m, path, 1)
        )
        results = {}
        try:
            for target in (0.0, 0.02, 0.125):
                z, state, rec = at_progress(
                    op, points, target, 1e-8 if engine == "native" else 1e-10, engine
                )
                if state is not None:
                    rec["fields"] = fields(
                        op, z, state, engine, m, rec.pop("arrival", None)
                    )
                results[str(target)] = rec
        finally:
            if engine == "native":
                op.close()
        report[engine] = results
    comparisons = {}
    for key in ("0.0", "0.02", "0.125"):
        a = report["native"].get(key, {}).get("fields")
        b = report["independent"].get(key, {}).get("fields")
        if (
            a
            and b
            and a["generations"] == b["generations"]
            and a["event"] == b["event"]
        ):
            comparisons[key] = difference(a, b)
        else:
            comparisons[key] = dict(
                reason="same-stage, same-side field pairing unavailable"
            )
    return dict(engines=report, common_progress=comparisons)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("bridge", type=Path)
    args = parser.parse_args()
    results = {
        str(n): compare(args.root / f"mesh-{n}", args.bridge.resolve())
        for n in (18, 144)
    }
    surfaces = {}
    for key in ("0.0", "0.02", "0.125"):
        a = results["18"]["engines"]["native"][key].get("fields")
        b = results["144"]["engines"]["native"][key].get("fields")
        if a and b:
            surfaces[key] = {
                s: {
                    k: np.asarray(a["surfaces"][s][k]) - b["surfaces"][s][k]
                    for k in a["surfaces"][s]
                }
                for s in a["surfaces"]
            }
    save(
        args.root / "comparison.json.gz",
        dict(grids=results, surface_18_minus_144=surfaces),
    )
