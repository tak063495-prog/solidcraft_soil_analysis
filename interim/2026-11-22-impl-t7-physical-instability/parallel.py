"""Re-evaluate the two accepted histories without exchanging checkpoints."""

import argparse
import copy
import gzip
import json
from pathlib import Path
import subprocess
import sys
import time
import hashlib
import platform
import scipy

import numpy as np

from work_balance import interval, native_state
from continuation import (
    assess_path,
    direction,
    Controls,
    Evaluation,
    TrialRejected,
    identity,
    serial,
    spectrum,
    trace,
)

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "interim/2026-12-28-a01-joint-combo-cloud-followup"))
from a01follow_binding import fixtures
from a01follow_contact import evaluate as differential
from a01follow_reference import Model, run
from verify_joint_coupling import contact
from wet_joint_cycle_reference import contact_work
from types import SimpleNamespace

U = 0.01  # One scaled translation = 0.01 m; lambda = stage progress.


def refine(mesh):
    """Red refinement of each Tet4. Keep coincident nodes on opposite sides separate."""
    fine = copy.deepcopy(mesh)
    xyz = fine["nodes"]
    mids = {}

    def mid(a, b):
        edge = tuple(sorted((a, b)))
        if edge not in mids:
            mids[edge] = len(xyz)
            xyz.append(((np.array(xyz[a]) + xyz[b]) / 2).tolist())
        return mids[edge]

    cells = []
    for a, b, c, d in mesh["tetrahedra"]:
        ab, ac, ad, bc, bd, cd = [
            mid(i, j) for i, j in ((a, b), (a, c), (a, d), (b, c), (b, d), (c, d))
        ]
        children = [
            [a, ab, ac, ad],
            [ab, b, bc, bd],
            [ac, bc, c, cd],
            [ad, bd, cd, d],
            [ab, ac, ad, cd],
            [ab, ad, bd, cd],
            [ab, bd, bc, cd],
            [ab, bc, ac, cd],
        ]
        for tet in children:
            pts = np.asarray([xyz[i] for i in tet])
            if np.linalg.det(pts[1:] - pts[0]) < 0:
                tet[1], tet[2] = tet[2], tet[1]
        cells.extend(children)
    fine["tetrahedra"] = cells
    for key in ("region_ids", "material_ids"):
        fine[key] = np.repeat(mesh[key], 8).tolist()
    for name, ids in mesh["node_sets"].items():
        subset = set(ids)
        fine["node_sets"][name] = ids + [
            m for (a, b), m in mids.items() if a in subset and b in subset
        ]
    for name, faces in mesh["face_sets"].items():
        parts = []
        for a, b, c in faces:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            parts.extend([[a, ab, ca], [ab, b, bc], [ca, bc, c], [ab, bc, ca]])
        fine["face_sets"][name] = parts
    return fine


def array_joint(fields):
    keys = {
        "jump": "jump",
        "traction": "traction",
        "rest": "elastic_offset",
        "accumulated": "accumulated",
        "plastic": "plastic_slip",
        "pd": "plastic_dissipation",
        "nd": "numerical_dissipation",
        "work": "endpoint_work",
        "initial": "initial_traction",
        "reference": "reference_jump",
        "initial_energy": "initial_energy",
        "status": "contact_status",
    }
    return {k: np.asarray(fields[v]).copy() for k, v in keys.items()}


def independent_start(model):
    events = {}
    accepted = {}

    def observe(e):
        if e["kind"] == "birth_preparation":
            events[e["stage"]] = e

    def sink(row, fields):
        accepted.update(row=row, fields=copy.deepcopy(fields))

    try:
        run(model, 8, sink, observer=observe)
    except (ArithmeticError, np.linalg.LinAlgError) as error:
        reason = str(error)
    else:
        raise ValueError("independent load control did not stop")
    stage = next(reversed(events))
    prep = events[stage]
    t, offset = model.spaces[False]
    row = accepted["row"]
    f = accepted["fields"]
    step = next(v for v in model.case["stages"] if v["name"] == stage)
    if step["kind"] == "boundary":
        origin = next(
            (
                v["target_factor"]
                for v in reversed(
                    model.case["stages"][: model.case["stages"].index(step)]
                )
                if v["kind"] == "boundary"
            ),
            0.0,
        )
        rate = step["target_factor"] - origin
        joint = array_joint(f)
        u = f["u"].ravel()
        progress = row["progress"] if row["stage"] == stage else 0.0
        generations = f["generation"].tolist()
    elif stage == "birth":
        origin = rate = 0.0
        generations = [1] * len(model.case["joint_system"]["joints"])
        if row["stage"] == "birth":
            joint = array_joint(f)
            u = f["u"].ravel()
            progress = row["progress"]
        else:
            if row["stage"] != "geostatic" or row["progress"] != 1.0:
                raise ValueError("missing pre-birth equilibrium")
            joint = prep["joint"]
            u = prep["u"]
            progress = 0.0
    else:
        raise ValueError(f"unhandled independent stage {stage}")
    parent = dict(u=u, joint=joint, generation=generations, progress=progress)
    path = dict(name=stage, origin=origin, rate=rate, holding=prep["holding"])
    return (
        parent,
        np.r_[t.T @ (u - (origin + progress * rate) * offset) / U, progress],
        path,
        reason,
        accepted,
    )


def history(load, reaction, counts, branches, u):
    return dict(
        load_resultant_kn=np.sum(np.reshape(load, (-1, 3)), axis=0),
        reaction_resultant_kn=np.sum(np.reshape(reaction, (-1, 3)), axis=0),
        contact_substeps=counts,
        control_displacement_norm_m=np.linalg.norm(u),
        contact_counts={
            name: sum(b["branch"] == name for b in branches)
            for name in ("open", "stick", "slip")
        },
    )


def derivatives(model, joint, target, counts, t, offset_rate, load_rate):
    _, deriv, branches = differential(model.system, joint, target, counts)
    cr = model.system.c @ t
    terms = np.einsum("qai,qab,q->iqb", cr, deriv, model.system.weight, optimize=True)
    full = t.T @ model.k + terms.reshape(len(t.T), -1) @ model.system.c.reshape(
        -1, model.nd
    )
    return full @ t, full @ offset_rate - t.T @ load_rate, branches


def independent(model, path, multiple):
    t, offset = model.spaces[False]
    system = model.system
    holding = path["holding"]
    rate = path["rate"]
    origin = path["origin"]

    def applied(progress):
        return (
            model.base
            + (origin + progress * rate) * model.boundary
            + (1 - progress) * holding
        )

    def evaluate(parent, z):
        if not 0 <= z[-1] <= 1:
            raise TrialRejected("stage progress outside [0,1]")
        old = parent["joint"]
        u = t @ (U * z[:-1]) + (origin + z[-1] * rate) * offset
        target = system.c @ u
        count = (
            np.maximum(
                1,
                np.ceil(
                    np.linalg.norm(target - old["jump"], axis=1)
                    / model.case["joint_system"]["maximum_increment_m"]
                ).astype(int),
            )
            * multiple
        )
        if np.max(count) > 50000:
            raise TrialRejected("contact partition budget")
        k, q, branches = derivatives(
            model, old, target, count, t, rate * offset, rate * model.boundary - holding
        )
        joint = copy.deepcopy(old)
        for n in np.unique(count):
            ix = count == n
            sub = SimpleNamespace(
                **{
                    key: getattr(system, key)[ix]
                    for key in ("normal", "kn", "kt", "mu")
                }
            )
            before = old["jump"][ix] - old["rest"][ix]
            after = target[ix] - old["rest"][ix]
            traction, plastic, slip = contact(
                sub, before, old["traction"][ix], after, int(n)
            )
            pd, nd, work = contact_work(sub, before, old["traction"][ix], after, int(n))
            joint["traction"][ix] = traction[:, 0]
            joint["plastic"][ix] += plastic[:, 0]
            joint["accumulated"][ix] += slip[:, 0]
            joint["pd"][ix] += pd
            joint["nd"][ix] += nd
            joint["work"][ix] += work
        joint["jump"] = target
        event = tuple(b["branch"] for b in branches)
        joint["status"] = np.array([("open", "stick", "slip").index(b) for b in event])
        load = applied(z[-1])
        internal = (
            model.k @ u
            + model.seed
            + np.einsum("qag,qa,q->g", system.c, joint["traction"], system.weight)
        )
        residual = t.T @ (internal - load)
        fields = model.fields(
            u,
            joint,
            True,
            [1, 1],
            [0, 0],
            load,
            internal,
            None,
            None,
            (1 - z[-1]) * holding,
            {},
        )
        state = dict(
            u=u,
            joint=joint,
            generation=parent["generation"],
            progress=z[-1],
            fields=fields,
        )
        work = interval(
            model,
            t,
            dict(u=parent["u"], **old),
            dict(u=u, **joint),
            applied(parent["progress"]),
            load,
        )
        return Evaluation(
            residual,
            k * U,
            q,
            state,
            max(1.0, np.linalg.norm(load)),
            eligible=0.0 <= z[-1] <= 1.0,
            event=event,
            details=dict(
                physical_k=k,
                q=q,
                counts=count,
                load=load,
                internal=internal,
                branches=branches,
                work_kn_m=work,
                history=history(load, internal - load, count, branches, u),
            ),
        )

    return evaluate


class Native:
    def __init__(self, exe, folder, model, run_data, multiple):
        self.process = subprocess.Popen(
            [
                str(exe),
                *[str(folder / (k + ".json")) for k in ("mesh", "ground", "case")],
            ],
            text=True,
            encoding="utf-8",
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        self.m = model
        self.multiple = multiple
        self.parent = run_data["last_committed"]
        self.stage = self.parent["stage_index"] + (self.parent["stage_progress"] == 1)
        self.item = model.case["stages"][self.stage - 1]
        if self.item["kind"] != "boundary":
            raise ValueError("no same-generation boundary continuation")
        self.origin = model.case["stages"][self.stage - 2]["target_factor"]
        self.rate = self.item["target_factor"] - self.origin
        self.system = model.system
        # Rust and independent quadrature may order points differently. Require
        # the same face, owner, physical location, normal and weight bijectively.
        layout = run_data["joint_integration_points"]
        self.order = []
        for p in self.system.points:
            match = [
                i
                for i, s in enumerate(layout)
                if all(
                    s[k] == p[k]
                    for k in ("joint_id", "face_index", "owner_a", "owner_b")
                )
                and np.linalg.norm(np.array(s["position_m"]) - p["position_m"]) < 1e-12
                and np.linalg.norm(np.array(s["normal"]) - p["normal"]) < 1e-12
                and abs(s["weight_m2"] - p["weight_m2"]) < 1e-14
            ]
            if len(match) != 1:
                raise ValueError("ambiguous native/physical quadrature correspondence")
            self.order.append(match[0])
        if len(set(self.order)) != len(layout):
            raise ValueError("incomplete quadrature mapping")
        p = (
            self.parent["stage_progress"]
            if self.parent["stage_index"] == self.stage
            else 0.0
        )
        first = self.query(self.parent, None, p)
        self.start = np.r_[np.asarray(first["free_m"]) / U, p]
        self.t = np.zeros((model.nd, len(self.start) - 1))
        for i, terms in enumerate(first["dof_rows"]):
            for j, c in terms:
                self.t[i, j] = c
        self.offset = model.spaces[False][1] * self.rate

    def query(self, parent, free, p):
        msg = dict(
            parent=parent,
            stage=self.stage,
            origin_factor=self.origin,
            progress=p,
            free_m=free,
            contact_subdivision=self.multiple,
        )
        self.process.stdin.write(
            json.dumps(msg, default=serial, allow_nan=False) + "\n"
        )
        self.process.stdin.flush()
        result = json.loads(self.process.stdout.readline())
        if "error" in result:
            if result["error"].startswith(
                ("material update failed:", "numerical solve failed:")
            ):
                raise TrialRejected(result["error"])
            raise ValueError(result["error"])
        return result

    def close(self):
        self.process.stdin.close()
        if self.process.wait() != 0:
            raise RuntimeError("native evaluation process failed")

    def __call__(self, parent, z):
        if not 0 <= z[-1] <= 1:
            raise TrialRejected("stage progress outside [0,1]")
        r = self.query(parent, z[:-1] * U, z[-1])
        state = r["candidate"]
        old = [parent["joints"]["points"][i] for i in self.order]
        joint = {
            k: np.array([p[field] for p in old])
            for k, field in [("jump", "jump_m"), ("traction", "traction_kpa")]
        }
        joint["rest"] = np.array([p["birth"]["elastic_offset_m"] for p in old])
        new = [state["joints"]["points"][i] for i in self.order]
        target = np.array([p["jump_m"] for p in new])
        counts = np.array(r["partitions"][-len(old) :])[self.order]
        k, q, branches = derivatives(
            self.m,
            joint,
            target,
            counts,
            self.t,
            self.offset,
            self.m.boundary * self.rate,
        )
        # Analytic bulk reconstruction is valid only on the elastic material branch.
        elastic = all(
            p.get("last_return") == "elastic" and p["equivalent_plastic_strain"] == 0.0
            for p in state["material_states"]
        )
        solver = np.array(r["tangent_columns"]).T
        work = (
            interval(
                self.m,
                self.t,
                native_state(parent, self.order),
                native_state(state, self.order),
                np.asarray(parent["applied_nodal_load_kn"]),
                np.asarray(state["applied_nodal_load_kn"]),
            )
            if elastic
            else None
        )
        return Evaluation(
            np.asarray(r["residual_kn"]),
            k * U,
            q,
            state,
            max(1.0, np.linalg.norm(state["applied_nodal_load_kn"])),
            0.0 <= z[-1] <= 1.0
            and elastic
            and r["supported"]
            and r["constraint_error_m"] <= 1e-10,
            tuple(b["branch"] for b in branches),
            solver * U,
            dict(
                physical_k=k,
                q=q,
                partitions=r["partitions"],
                branches=branches,
                work_kn_m=work,
                reaction=r["reaction_kn"],
                supported=r["supported"],
                dof_rows=r["dof_rows"],
                history=history(
                    state["applied_nodal_load_kn"],
                    r["reaction_kn"],
                    counts,
                    branches,
                    state["displacement_m"],
                ),
            ),
        )


def save(path, value):
    # Keep the small restart header first. Encode each critical trial in C's
    # JSON encoder separately so neither startup nor sealing needs a full
    # second copy of the history in memory.
    def encode(v):
        return json.dumps(
            v, default=serial, allow_nan=False, separators=(",", ":")
        ).encode()

    with (
        path.open("xb") as raw,
        gzip.GzipFile(
            filename="", fileobj=raw, mode="wb", compresslevel=1, mtime=0
        ) as f,
    ):
        if "critical" not in value:
            f.write(encode(value))
            return
        header = {
            k: v
            for k, v in value.items()
            if k not in ("critical", "accepted", "trials")
        }
        header["last_point"] = value["accepted"][-1]
        f.write(encode(header)[:-1])
        for key in ("accepted", "trials"):
            f.write(b"," + encode(key) + b":[")
            for i, item in enumerate(value[key]):
                if i:
                    f.write(b",")
                f.write(encode(item))
            f.write(b"]")
        f.write(b',"critical":{')
        for i, (key, item) in enumerate(value["critical"].items()):
            if i:
                f.write(b",")
            f.write(encode(str(key)) + b":" + encode(item))
        f.write(b"}}")


def series(evaluator, parent, z, ds, steps, tol):
    direction = np.zeros(len(z))
    direction[-1] = 1
    return trace(
        evaluator,
        parent,
        z,
        direction,
        Controls(ds, steps, tol, 1e-10, 18, event_bisections=4),
        stop=lambda points: points[-1]["z"][-1] >= 0.125,
    )


def opening_check(op, points, model, t, offset, tol):
    # Physical area average on LowerSeam, normal directed from a to b.
    sys = model.system
    mask = np.array([p["joint_id"] == "LowerSeam" for p in sys.points])
    row = np.einsum(
        "q,qa,qag->g", sys.weight[mask], sys.normal[mask], sys.c[mask]
    ) / sum(sys.weight[mask])
    plane = np.r_[row @ t * U, row @ offset] / U
    for a, b, c in zip(points, points[1:], points[2:]):
        za, zb, zc = [p["z"] for p in (a, b, c)]
        if not (
            a["event"] == b["event"] == c["event"]
            and (plane @ (zb - za)) * (plane @ (zc - zb)) > 0
        ):
            continue
        parent = b["state"]
        incoming = zc - zb
        try:
            tangent = direction(op(parent, zb), incoming)
        except ArithmeticError:
            continue
        if abs(plane @ tangent) <= np.finfo(float).eps * np.linalg.norm(plane):
            continue
        ds = (plane @ incoming) / (plane @ tangent)
        if ds <= 0:
            continue
        check = trace(
            op,
            parent,
            zb,
            incoming,
            Controls(ds, 1, tol, 1e-10, 18, event_bisections=0),
            local_control=plane,
        )
        ok = (
            len(check["accepted"]) == 2 and check["accepted"][-1]["event"] == c["event"]
        )
        return dict(
            interval=[zb, zc],
            direction=plane,
            trace=check,
            matched_side=ok,
            difference_scaled=check["accepted"][-1]["z"] - zc if ok else None,
        )
    return dict(
        matched_side=False,
        reason="no three same-side equilibria with monotone LowerSeam opening",
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    parser.add_argument("cli", type=Path)
    parser.add_argument("bridge", type=Path)
    parser.add_argument("--elements", type=int, choices=(18, 144), required=True)
    args = parser.parse_args()
    d = args.output
    d.mkdir(parents=True, exist_ok=False)
    mesh, ground, case = fixtures.stacked(steps=8)
    if args.elements == 144:
        mesh = refine(mesh)
    for name, v in zip(("mesh", "ground", "case"), (mesh, ground, case)):
        (d / (name + ".json")).write_text(
            json.dumps(v, default=serial), encoding="utf-8"
        )
    m = Model(mesh, ground, case)
    t0 = time.monotonic()
    parent, z, path, termination, accepted = independent_start(m)
    save(
        d / "independent-last.json.gz",
        dict(checkpoint=accepted, stage_origin=parent, stage_path=path),
    )
    with (d / "native.log").open("w") as log:
        command = [
            str(args.cli.resolve()),
            "solve-stages",
            *[str(d / (n + ".json")) for n in ("mesh", "ground", "case", "run")],
        ]
        rc = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT).returncode
    native_run = json.loads((d / "run.json").read_bytes())
    summary = dict(
        elements=args.elements,
        independent_load_stop=termination,
        native_load_stop=native_run["failure_reason"],
        native_exit=rc,
        prehistory_seconds=time.monotonic() - t0,
        series=[],
    )
    schedule = (
        ((0.04, 1), (0.02, 1), (0.01, 1), (0.02, 2), (0.02, 4))
        if args.elements == 18
        else ((0.02, 1),)
    )
    for engine in ("independent", "native"):
        for ds, multiple in schedule:
            label = f"{engine}-{ds}-{multiple}"
            started = time.monotonic()
            if (
                engine == "native"
                and case["stages"][
                    native_run["last_committed"]["stage_index"]
                    + (native_run["last_committed"]["stage_progress"] == 1)
                    - 1
                ]["kind"]
                != "boundary"
            ):
                summary["series"].append(
                    dict(
                        label=label,
                        reason="last checkpoint at a generation transition",
                        physical="D",
                    )
                )
                continue
            op = (
                independent(m, path, multiple)
                if engine == "independent"
                else Native(args.bridge.resolve(), d, m, native_run, multiple)
            )
            start, z0 = (
                (parent, z) if engine == "independent" else (op.parent, op.start)
            )
            initial = op(start, z0)
            try:
                result = series(
                    op,
                    start,
                    z0,
                    ds,
                    64,
                    1e-10
                    if engine == "independent"
                    else case["settings"]["relative_force_tolerance"],
                )
                tol = (
                    1e-10
                    if engine == "independent"
                    else case["settings"]["relative_force_tolerance"]
                )
                basis, rate = (
                    (m.spaces[False][0], m.spaces[False][1] * path["rate"])
                    if engine == "independent"
                    else (op.t, op.offset)
                )
                result["fold_candidates"] = assess_path(op, result["accepted"], tol)
                result["opening_control"] = opening_check(
                    op, result["accepted"], m, basis, rate, tol
                )
            finally:
                if engine == "native":
                    op.close()
            result.update(
                stage=path["name"] if engine == "independent" else op.item["name"],
                subdivision=multiple,
                input_sha={
                    name: identity(v)
                    for name, v in zip(("mesh", "ground", "case"), (mesh, ground, case))
                },
                initial_spectrum=spectrum(initial.k / U, initial.q, initial.residual),
            )
            for record in result["critical"].values():
                if "iterations" not in record:
                    continue
                ev = next(
                    (
                        i["evaluation"]
                        for i in reversed(record["iterations"])
                        if "evaluation" in i
                    ),
                    None,
                )
                if ev is None:
                    continue
                ev["spectrum"] = spectrum(
                    ev["details"]["physical_k"], ev["q"], ev["residual"]
                )
            result["environment"] = dict(
                python=platform.python_version(),
                numpy=np.__version__,
                scipy=scipy.__version__,
                cli_sha256=hashlib.sha256(args.cli.read_bytes()).hexdigest(),
                bridge_sha256=hashlib.sha256(args.bridge.read_bytes()).hexdigest(),
                sources={
                    p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in (
                        Path(__file__),
                        Path(__file__).with_name("continuation.py"),
                        Path(__file__).with_name("work_balance.py"),
                    )
                },
                displacement_scale_m=U,
                load_coefficient_scale=1,
                force_unit="kN",
                length_unit="m",
            )
            save(d / (label + ".json.gz"), result)
            summary["series"].append(
                dict(
                    label=label,
                    stop=result["termination_reason"],
                    physical=result["physical_classification"],
                    fold_candidates=[c["kind"] for c in result["fold_candidates"]],
                    monotone_opening_match=result["opening_control"]["matched_side"],
                    opening_difference_max_scaled=float(
                        np.max(abs(result["opening_control"]["difference_scaled"]))
                    )
                    if result["opening_control"]["matched_side"]
                    else None,
                    attempts=len(result["trials"]),
                    accepted=len(result["accepted"]) - 1,
                    last_progress=result["accepted"][-1]["z"][-1],
                    relative_force=max(
                        (
                            r["relative_force"]
                            for r in result["trials"]
                            if r["accepted"]
                        ),
                        default=None,
                    ),
                    seconds=time.monotonic() - started,
                )
            )
            (d / "summary.json").write_text(
                json.dumps(summary, default=serial, indent=2), encoding="utf-8"
            )
    (d / "summary.json").write_text(
        json.dumps(summary, default=serial, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
