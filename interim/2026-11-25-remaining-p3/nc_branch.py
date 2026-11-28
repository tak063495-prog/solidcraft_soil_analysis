"""Original MCC section from the stored simple fold, with adaptive pseudo-arclength."""
import gzip
import json
from pathlib import Path
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
P2 = HERE.parent / "2026-11-25-remaining5-p2"
sys.path.insert(0, str(P2))
from nc_continue import Section


def correct(m, states, lam, tangent, ds, tol):
    dz = ds * tangent
    counts = np.ones(len(states), dtype=int)
    iterations = []
    for it in range(m.settings["maximum_newton_iterations"] + 1):
        try:
            evaluated = m.evaluate(states, dz[:-1], counts)
            if evaluated is None:
                iterations.append(dict(repartition=int(counts.max())))
                dz = ds * tangent
                continue
            f, k, new, plastic = evaluated
            load = lam + dz[-1]
            ext = m.load(load)
            rhs = ext[m.free] - f
            rel = np.linalg.norm(rhs) / max(1., np.linalg.norm(ext))
            arc = ds - tangent @ dz
            a = np.vstack((np.column_stack((k, m.h)), tangent))
            row = dict(load=load, relative_force=rel, arc_residual_m=arc)
            iterations.append(row)
            if any("linear" in item for item in iterations) and rel <= tol and abs(arc) <= 1e-12:
                left, s, vh = np.linalg.svd(k)
                w = left[:, -1] * (np.sign(left[:, -1] @ m.h) or 1)
                v = vh[-1] * (np.sign(vh[-1] @ tangent[:-1]) or 1)
                border = np.linalg.svd(a, compute_uv=False)[-1]
                next_t = np.linalg.solve(a, np.r_[np.zeros(len(k)), 1.])
                next_t /= np.linalg.norm(next_t)
                diag = dict(relative_force=rel, lambda_tangent=next_t[-1], sigma_min=s[-1],
                            sigma_next=s[-2], bordered_min=border, projection_kn=w @ m.h,
                            signed_stiffness=w @ k @ v, plastic_points=plastic,
                            partitions_max=int(counts.max()))
                return (dz, new, diag), iterations, None
            b = np.r_[rhs, arc]
            dx = np.linalg.solve(a, b)
            row["linear"] = np.linalg.norm(a @ dx - b) / max(1e-30, np.linalg.norm(b))
            row["physical_linear"] = np.linalg.norm(k @ dx[:-1] + m.h * dx[-1] - rhs) / max(1e-30, np.linalg.norm(rhs))
            if max(row["linear"], row["physical_linear"]) > m.settings["relative_linear_tolerance"]:
                raise ArithmeticError("true linear residual")
            dz += dx
        except (ArithmeticError, np.linalg.LinAlgError) as error:
            return None, iterations, str(error)
    return None, iterations, "Newton iteration limit"


def run(inputs, binary, out, max_ds):
    out.mkdir(parents=True)
    start = time.monotonic()
    m = Section(inputs, binary)
    fold = json.loads((P2 / "nc/fold.json").read_bytes())
    with gzip.open(P2 / "nc/0.0002/accepted.jsonl.gz", "rt") as f:
        for line in f:
            parent = json.loads(line)
            if parent["step"] == fold["parent_step"]:
                break
            previous = parent
    def point(row):
        return np.r_[np.ravel(row["displacement_m"])[m.free], row["load"]]
    t = point(parent) - point(previous)
    t /= np.linalg.norm(t)
    cert = fold["certificate"]
    root, log, failure = correct(m, parent["active_material_states"], parent["load"], t, cert["local_arc_m"], 1e-13)
    assert failure is None, failure
    dz, old, check = root
    u = np.ravel(parent["displacement_m"]).copy()
    u[m.free] += dz[:-1]
    lam = parent["load"] + dz[-1]
    assert abs(lam - cert["load"]) < 1e-13 and check["sigma_min"] < 1e-7
    report = dict(start_lambda=lam, source_sha256=m.run_sha, fold_replay=check,
                  fold_iterations=log, max_step_m=max_ds, force_gate=m.settings["relative_force_tolerance"],
                  linear_gate=m.settings["relative_linear_tolerance"], mesh_elements=96,
                  budget=16520, original_attempt_cap=2000)
    t = dz / np.linalg.norm(dz)
    ds, arc, accepted = min(0.0008, max_ds), 0., 0
    stop = "original attempt cap"
    def checkpoint(name):
        with gzip.open(out / (name + ".json.gz"), "wt") as f:
            json.dump(dict(load=lam, arc_m=arc, displacement_m=u.reshape(-1,3).tolist(),
                           active_material_states=old), f, allow_nan=False)
    checkpoint("fold")
    last_dl = None
    last_slope = None
    turn_parent = None
    with (out / "trials.jsonl").open("x", buffering=1) as stream:
        for attempt in range(1, 2000):
            result, log, reason = correct(m, old, lam, t, ds, min(m.settings["relative_force_tolerance"], 1e-12))
            row = dict(attempt=attempt, parent_lambda=lam, parent_arc_m=arc, ds_m=ds,
                       iterations=log, accepted=result is not None)
            if result is None:
                row["reason"] = reason
                ds /= 2
                if ds < m.settings["minimum_increment"]:
                    stop = "failed corrector at minimum step"
            else:
                dz, new, diag = result
                if last_slope is not None and last_slope < 0 <= diag["lambda_tangent"]:
                    turn_parent = dict(load=lam, arc_m=arc, displacement_m=u.reshape(-1,3).tolist(),
                                       active_material_states=old, tangent=t.tolist(), ds_m=ds)
                    with gzip.open(out / "turn-parent.json.gz", "wt") as f:
                        json.dump(turn_parent, f, allow_nan=False)
                last_slope = diag["lambda_tangent"]
                old = new
                u[m.free] += dz[:-1]
                lam += dz[-1]
                arc += ds
                accepted += 1
                t = dz / np.linalg.norm(dz)
                row.update(diag, load=lam, arc_m=arc, delta_lambda=dz[-1])
                if last_dl is not None and last_dl < 0 < dz[-1]:
                    stop = "next load reversal bracketed"
                elif lam >= 1:
                    stop = "lambda=1 bracketed"
                elif lam <= 0:
                    stop = "lower boundary of physical release interval crossed"
                last_dl = dz[-1]
                ds = min(max_ds, ds * (1.3 if len(log) <= 8 else 1.))
            stream.write(json.dumps(row, allow_nan=False) + "\n")
            if stop != "original attempt cap":
                break
    checkpoint("last")
    report.update(stop=stop, accepted=accepted, attempts=attempt, last_lambda=lam, arc_m=arc,
                  seconds=time.monotonic()-start)
    (out / "summary.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    run(Path(sys.argv[1]), Path(sys.argv[2]).resolve(), Path(sys.argv[3]), float(sys.argv[4]))
