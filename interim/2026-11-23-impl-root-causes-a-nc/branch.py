"""Opt-in dry A-01 continuation from installed-stick BIRTH through shear into open."""
import argparse
import copy
import gzip
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import tempfile
import sys

import numpy as np
import scipy

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE.parent / "2026-11-22-impl-t7-physical-instability"))
from continuation import Controls, TrialRejected, direction, identity, serial, trace
from parallel import Model, Native, U, array_joint, fixtures, independent, refine, run
from compare import checkpoint

LIMIT = 64  # Original T7 open interval, across all resumed segments.


def read(p):
    return json.loads(p.read_bytes())


def sha(p):
    with p.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def write(p, data):
    with p.open("xb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as stream:
        stream.write(json.dumps(data, default=serial, allow_nan=False,
                                separators=(",", ":")).encode())


def environment(folder, bridge):
    sources = {Path(__file__).resolve(), ROOT / "tools/joint_tip_fixture.py",
               ROOT / "tools/split_joint_mesh.py"}
    for mod in list(sys.modules.values()):
        name = getattr(mod, "__file__", None)
        if name and name.endswith(".py") and Path(name).resolve().is_relative_to(ROOT):
            sources.add(Path(name).resolve())
    return dict(inputs={k: sha(folder / (k + ".json")) for k in ("mesh", "ground", "case")},
                sources={p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(sources)},
                bridge=sha(bridge), python=platform.python_version(), numpy=np.__version__,
                scipy=scipy.__version__)


class OpenOrigin(Exception):
    """Stop the prehistory before the first open increment."""


def native_history(cli, inputs, report):
    result = []
    offset = sum(a["accepted"] for a in report["geostatic"]["attempts"])
    attempts = [a for a in report["attempts"] if a["accepted"] and a["stage"] in ("birth", "shear")]
    with tempfile.TemporaryDirectory() as tmp:
        for count, attempt in enumerate(attempts, 1):
            path = Path(tmp) / f"{count}.json"
            subprocess.run([cli, "checkpoint-stages", *inputs, str(offset + count), path],
                           check=True, stdout=subprocess.DEVNULL)
            run_data = read(path)["run"]
            state = run_data["last_committed"]
            assert state["stage_progress"] == attempt["progress_to"]
            assert run_data["attempts"][-1]["stage"] == attempt["stage"]
            result.append(dict(stage=attempt["stage"], state=state))
    return result


def prepare(args):
    d = args.folder
    d.mkdir(parents=True, exist_ok=False)
    mesh, ground, case = fixtures.stacked(steps=8)
    if args.elements == 144:
        mesh = refine(mesh)
    files = [d / (k + ".json") for k in ("mesh", "ground", "case")]
    for p, data in zip(files, (mesh, ground, case), strict=True):
        p.write_text(json.dumps(data, default=serial), encoding="utf-8")
    native = d / "native.json"
    proc = subprocess.run([args.cli.resolve(), "solve-stages", *files, native], capture_output=True)
    (d / "native.log").write_bytes(proc.stdout + proc.stderr)
    report = read(native)
    if not report["failure_reason"] and not report["complete"]:
        raise ValueError("missing native stop record")
    # Keep the actual native shear checkpoint. Later native failures stay in native.json.
    origin, = [f["state"] for f in report["frames"]
               if f["name"] == "shear" and f["state"]["stage_progress"] == 1]
    history = []
    path = {}
    model = Model(mesh, ground, case)

    def observe(event):
        if event["kind"] == "birth_preparation" and event["stage"] == "open":
            path["holding"] = event["holding"]
            raise OpenOrigin

    def sink(row, fields):
        history.append(dict(row=row, fields=copy.deepcopy(fields)))

    try:
        run(model, 8, sink, observer=observe, birth_predictor="installed_stick")
    except OpenOrigin:
        pass
    end = history[-1]
    assert (end["row"]["stage"], end["row"]["progress"]) == ("shear", 1.0)
    path.update(name="open", origin=end["row"]["factor"])
    opened = next(s for s in case["stages"] if s["name"] == "open")
    path["rate"] = opened["target_factor"] - path["origin"]
    f = end["fields"]
    parent = dict(u=f["u"].ravel(), joint=array_joint(f), generation=f["generation"], progress=0.)
    t, offset = model.spaces[False]
    z = np.r_[t.T @ (parent["u"] - path["origin"] * offset) / U, 0.]
    # Obtain the admitted contact sides on the published short interval
    # from the independent input-only solve. Pass labels only to the verifier.
    short_case = copy.deepcopy(case)
    short_case["stages"] = short_case["stages"][:3]
    short_case["stages"][-1]["target_factor"] = path["origin"] + 0.02 * path["rate"]
    short = []

    def short_point(row, fields):
        if row["stage"] == "open":
            short.append(dict(progress=0.02 * row["progress"], factor=row["factor"],
                              generation=fields["generation"].tolist(),
                              event=[("open", "stick", "slip")[i] for i in fields["contact_status"]]))

    run(Model(mesh, ground, short_case), 8, short_point, birth_predictor="installed_stick")
    assert len(short) == 8 and short[-1]["progress"] == 0.02
    native_points = native_history(args.cli.resolve(), files, report)
    assert native_points[-1]["state"] == origin
    write(d / "prehistory.gz", dict(independent=history, native=native_points, native_cli_sha256=sha(args.cli),
                                     native_exit=proc.returncode, native_stop=report["failure_reason"]))
    write(d / "origins.gz", dict(native=origin, independent=parent, independent_z=z, path=path, short=short,
                                 prehistory=sha(d / "prehistory.gz"), native_run=sha(native),
                                 environment=environment(d, args.bridge)))


def unpack(p):
    with gzip.open(p) as stream:
        return json.load(stream)


def advance(args):
    d = args.folder
    origin = unpack(d / "origins.gz")
    binding = environment(d, args.bridge)
    if (origin["environment"] != binding or sha(d / "native.json") != origin["native_run"]
            or sha(d / "prehistory.gz") != origin["prehistory"]):
        raise ValueError("prehistory, executable, input or source bytes changed")
    model = Model(*(read(d / (k + ".json")) for k in ("mesh", "ground", "case")))
    engine = args.engine
    if engine == "native":
        native = read(d / "native.json")
        native["last_committed"] = origin["native"]
        op = Native(args.bridge.resolve(), d, model, native, 1)
        parent, z = op.parent, op.start
        tol = model.case["settings"]["relative_force_tolerance"]
    else:
        path = origin["path"]
        path["holding"] = np.asarray(path["holding"])
        op = independent(model, path, 1)
        parent, z = checkpoint(origin[engine], engine), np.array(origin["independent_z"])
        tol = 1e-10
    try:
        short = origin["short"]
        fixed = tuple(short[-1]["event"])
        restart = None
        used = 0
        chain = sha(d / "origins.gz")
        if args.resume:
            src, digest = args.resume
            src = Path(src)
            if sha(src) != digest:
                raise ValueError("parent segment SHA-256 mismatch")
            old = unpack(src)
            if (old["origins"] != chain or old["engine"] != engine or old["stop"] != "segment end"
                    or tuple(old["contact_branch"]) != fixed):
                raise ValueError("incompatible or terminal parent segment")
            chain = digest
            last = old["points"][-1]
            parent, z = checkpoint(last["state"], engine), np.asarray(last["z"])
            restart, used = old["restart"], old["steps_used"]
        steps = LIMIT - used if args.steps is None else args.steps
        if not 1 <= steps <= LIMIT - used:
            raise ValueError("segment exceeds remaining original open-step budget")

        def evaluate(state, coords):
            ev = op(state, coords)
            prior = state.get("progress") if engine == "independent" else (
                state["stage_progress"] if state["stage_index"] == op.stage else 0.)
            # A zero-increment replay can relabel a friction point as stick.
            # Check the arrival side at equilibrium, allowing trial iterates to cross.
            if (coords[-1] != prior and np.linalg.norm(ev.residual) <= tol * ev.denominator
                    and tuple(ev.event) != expected):
                raise TrialRejected("equilibrated candidate differs from accepted contact branch")
            return ev

        points, trials = [], []
        stop = "segment end"
        for _ in range(steps):
            plane = None
            ds = 0.02
            cuts = 4
            forward = np.zeros_like(z)
            forward[-1] = 1.
            expected = fixed
            if used < len(short):
                sample = short[used]
                expected = tuple(sample["event"])
                # Follow the original physical loads at the eight known short-path
                # coordinates. No displacement or force is imported as the solution.
                plane = forward.copy()
                incoming = op(parent, z)
                if restart:
                    incoming.k = np.asarray(restart["payload"]["k"])
                    incoming.q = np.asarray(restart["payload"]["q"])
                    forward = np.asarray(restart["payload"]["previous_direction"])
                tangent = direction(incoming, forward)
                if tangent[-1] <= 0:
                    stop = "no forward tangent to the next accepted load coordinate"
                    break
                ds = (sample["progress"] - z[-1]) / tangent[-1]
                cuts = 0
            step = trace(evaluate, parent, z, forward, Controls(ds, 1, tol, 1e-10, 18,
                         event_bisections=cuts), restart=restart, local_control=plane)
            trials.extend(step["trials"])
            if len(step["accepted"]) == 1:
                stop = step["termination_reason"]
                break
            pt = step["accepted"][-1]
            if tuple(pt["event"]) != expected:
                raise ValueError("unmatched arrival side")
            arrival = evaluate(parent, pt["z"])
            pt["arc_increment"] = pt.pop("s")
            pt["work_kn_m"] = arrival.details["work_kn_m"]
            pt["factor"] = origin["path"]["origin"] + pt["z"][-1] * origin["path"]["rate"]
            pt["relative_force"] = np.linalg.norm(arrival.residual) / arrival.denominator
            points.append(pt)
            parent, z, restart = pt["state"], pt["z"], step["restart"]
            used += 1
            if z[-1] >= 0.125:
                stop = "T7 observation reached"
                break
        if used == LIMIT and stop == "segment end":
            stop = "step budget"
        write(args.output, dict(origins=sha(d / "origins.gz"), engine=engine, parent_sha256=chain,
              stage="open", contact_branch=fixed, points=points, trials=trials, restart=restart,
              last_accepted=identity(parent), last_progress=z[-1], steps_used=used, stop=stop))
        print(engine, used, z[-1], stop)
    finally:
        if engine == "native":
            op.close()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    commands = p.add_subparsers(dest="command", required=True)
    start = commands.add_parser("prepare")
    start.add_argument("--elements", type=int, choices=(18, 144), required=True)
    start.add_argument("cli", type=Path)
    follow = commands.add_parser("follow-same-contact")
    follow.add_argument("engine", choices=("native", "independent"))
    follow.add_argument("output", type=Path)
    follow.add_argument("--steps", type=int)
    follow.add_argument("--resume", nargs=2, metavar=("SEGMENT", "SHA256"))
    for command in (start, follow):
        command.add_argument("folder", type=Path)
        command.add_argument("bridge", type=Path)
    args = p.parse_args()
    (prepare if args.command == "prepare" else advance)(args)


if __name__ == "__main__":
    main()
