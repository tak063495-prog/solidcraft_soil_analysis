"""Observe same-parent, fixed-partition equilibria along the existing independent route."""

from pathlib import Path
import inspect
import json
import sys
import numpy as np

root = Path.cwd()
sys.path.insert(0, str(root / "interim/2026-11-23-c01-method-improvements"))
from run_case import context
from methods import Methods, active_tangent

study, mesh, ground, case, control = context(1024)
ref = study.reference
control.update(
    localize_friction_events=True,
    resolve_slip_path=True,
    semismooth_contact=True,
    contact_history_bisection=True,
)
case["settings"]["attempt_budget_policy"] = "mesh_proportional"
seen = set()
aux = {}


def observe(v):
    name = v["name"]
    x = v["free"]
    sys_ = v["system"]
    if v["probe"] or name in seen:
        return
    if name == "open":
        if not (
            8 / 256 <= v["previous"] < 9 / 256
            and v["status"][5] == 2
            and v["new_status"][5] == 0
        ):
            return
    elif name == "reverse":
        if not np.any((v["slip"] * v["plastic"]).sum(1) < -1e-16):
            return
    elif name == "reload":
        if not 926 / 1024 <= v["previous"] < 927 / 1024:
            return
    else:
        return
    seen.add(name)
    evaluate = v["evaluate"]

    def jac(target):
        tangent = active_tangent(
            sys_, v["jump"], v["traction"], target, v["frozen_counts"]
        )
        return v["kr"] + np.einsum(
            "qai,qab,qbj,q->ij", v["cr"], tangent, v["cr"], sys_.weight, optimize=True
        )

    k = jac(v["target"])
    _, s, vt = np.linalg.svd(k)
    eig = np.linalg.eigvalsh((k + k.T) / 2)
    rows = []
    for axis in (vt[-1], vt[0]):
        for sign in (-1.0, 1.0):
            y = x + sign * 1e-5 * axis
            for it in range(20):
                _, j, q, flow, _, r = evaluate(y)
                if np.linalg.norm(r) < 1e-10 * v["scale"]:
                    break
                dx = np.linalg.solve(jac(j), -r)
                for cut in range(16):
                    trial = y + 2.0 ** (-cut) * dx
                    if np.linalg.norm(evaluate(trial)[-1]) < np.linalg.norm(r):
                        y = trial
                        break
                else:
                    break
            _, j, q, flow, _, r = evaluate(y)
            rows.append(
                dict(
                    res=np.linalg.norm(r) / v["scale"],
                    iterations=it,
                    du=np.linalg.norm(y - x),
                    flow_diff=np.max(abs(flow - v["slip"])),
                    partition_sufficient=bool(
                        np.all(v["partition"](j) <= v["frozen_counts"])
                    ),
                )
            )
    aux[name] = len(rows)
    print(
        json.dumps(
            dict(
                stage=name,
                parent=v["previous"],
                target=v["progress"],
                smin=s[-1],
                smax=s[0],
                sym_min=eig[0],
                runs=rows,
            )
        ),
        flush=True,
    )


source = inspect.getsource(ref._reference)
anchor = "                pd, nd, work = work_parts(target)"
assert source.count(anchor) == 1
# Evaluate frozen residual closures only. The route receives no replay results.
ref.__dict__["observe"] = observe
exec(
    compile(
        source.replace(anchor, "                observe(locals())\n" + anchor),
        ref.__file__,
        "exec",
    ),
    ref.__dict__,
)
counts = {}


def record(attempt):
    name = attempt["stage"]
    counts[name] = counts.get(name, 0) + 1
    assert counts[name] + aux.get(name, 0) <= 2065


ref._reference(
    mesh,
    ground,
    case,
    1024,
    1,
    7,
    False,
    True,
    maximum_m=5e-5,
    path_error=control,
    attempt_sink=record,
    methods=Methods(None, None, True, None),
)
print(json.dumps(dict(route_attempts=counts, diagnostic_solves=aux)), flush=True)
