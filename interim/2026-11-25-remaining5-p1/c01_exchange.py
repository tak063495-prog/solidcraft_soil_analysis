"""Cross the native and independent equilibrated jumps with both actual contact maps."""

import gzip
import inspect
import json
import subprocess
from pathlib import Path
import sys
from types import SimpleNamespace
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "2026-11-23-c01-method-improvements"))
from run_case import context
from methods import Methods, active_tangent

work, binary, output = map(Path, sys.argv[1:])
study, mesh, ground, case, control = context(1024)
ref = study.reference
control.update(
    localize_friction_events=True,
    resolve_slip_path=True,
    semismooth_contact=True,
    contact_history_bisection=True,
)
replay_route = work / "c01-route.json"
with gzip.open(HERE.parent / "2026-11-25-remaining5-diagnosis/c01-native.json.gz") as f:
    route = json.load(f)
if replay_route.exists():
    assert json.loads(replay_route.read_bytes()) == route
packed = work / "c01-paths.jsonl.gz"
if packed.exists():
    raw = gzip.open(packed, "rt")
else:
    raw = (work / "c01-paths.jsonl").open()
trials = []
for attempt in route["records"]:
    triple = [json.loads(next(raw)) for _ in range(3)]
    c, m, f = triple
    assert (
        c["from"] == m["from"] == attempt["from"]
        and c["to"] == f["to"] == attempt["to"]
    )
    assert (
        m["to"] == f["from"] and c["parent"] == m["parent"] and m["end"] == f["parent"]
    )
    if not attempt["accepted"] and attempt["event"] is None:
        trials.append((attempt, triple))
assert next(raw, None) is None and len(trials) == 70
# Retain the original independent volume assembly, constraints, contact model,
# seed, linear solve and equilibrium tolerance. Supply fixed parents/counts
# and visit only the three requested targets, bypassing the route controller.
source = inspect.getsource(ref._reference)
i = source.index("    phase('geostatic'")
source = (
    source[:i]
    + """
    def solve(parent,row):
        nonlocal u,jump,traction,status,forced_counts
        u,jump,traction,status = parent
        forced_counts = np.array(row["counts"])
        return increment(base+row["factor"]*boundary,row["factor"],"reload",row["to"],row["from"],probe=True)
    forced_counts = None
    def residual_correction(parent,row):
        nonlocal jump,traction,frozen_counts
        _,jump,traction,_ = parent
        frozen_counts = np.array(row["counts"])
        candidate = np.array(row["u"]).ravel()
        target = np.einsum('qag,g->qa',system.c,candidate)
        q = project(target)[0][:,0]
        load = base+row["factor"]*boundary
        residual = t.T@(k@candidate+seed_force+np.einsum('qag,qa,q->g',system.c,q,system.weight)-load)
        tangent = active_tangent(system,jump,traction,target,frozen_counts)
        jac = kr + np.einsum("qai,qab,qbj,q->ij",cr,tangent,cr,system.weight,optimize=True)
        dx = np.linalg.solve(jac,-residual)
        new = candidate + t@dx
        j1 = np.einsum('qag,g->qa',system.c,new)
        tr = project(j1)[0][:,0]
        return (new,j1,tr,None),np.linalg.norm(residual)/max(np.linalg.norm(load),1.)
    return solve,system,residual_correction
"""
)
marker = "        fixed = offset(factor);"
assert source.count(marker) == 1
source = source.replace(
    marker, "        frozen_counts = forced_counts.copy()\n" + marker
)
ns = dict(ref.__dict__)
exec(compile(source, ref.__file__, "exec"), ns)
ns["active_tangent"] = active_tangent
solve, system, residual_correction = ns["_reference"](
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
    methods=Methods(None, None, True, None),
)
one = SimpleNamespace(
    **{k: getattr(system, k)[:1] for k in ("kn", "kt", "mu", "normal")}
)
mat = case["joint_system"]["joints"][0]["material"]
rows, requests = [], []
for attempt, triple in trials:
    c, m, f = triple
    states = c["parent"]
    j, q = [np.array([s[key] for s in states]) for key in ("jump_m", "traction_kpa")]
    labels = np.array([("open", "stick", "slip").index(s["status"]) for s in states])
    parent = np.array(c["parent_u"]).ravel(), j, q, labels
    np.testing.assert_allclose(
        np.einsum("qag,g->qa", system.c, parent[0]), j, rtol=0, atol=5e-16
    )
    np.testing.assert_allclose(
        system.normal, [s["normal"] for s in states], rtol=0, atol=5e-15
    )
    ic, im = solve(parent, c), solve(parent, m)
    iff = solve(tuple(im[k] for k in ("u", "jump", "traction", "status")), f)
    ind = [ic, im, iff]
    for r, step in zip(ind, triple):
        np.testing.assert_array_equal(r["counts"], step["counts"])
    counts = [r["counts"][0] for r in triple]
    paths = [
        [r["end"][0]["jump_m"] for r in triple],
        [r["jump"][0].tolist() for r in ind],
    ]
    cc, rc = residual_correction(parent, c)
    cm, rm = residual_correction(parent, m)
    cf, rf = residual_correction(cm, f)
    paths.append([v[1][0].tolist() for v in (cc, cm, cf)])
    native_actual = [r["end"][0]["plastic_slip_m"][0] for r in triple]
    p = states[0]
    assert p.get("birth") is None
    maps = []
    for jumps in paths:
        requests.append(dict(material=mat, parent=p, jumps=jumps, counts=counts))

        def update(start, tr, end, n):
            tr, flow, _ = ref.contact(
                one, np.array([start]), np.array([tr]), np.array([end]), n
            )
            return tr[0, 0], flow[0, 0, 0]

        _, sc = update(j[0], q[0], jumps[0], counts[0])
        mt, sm = update(j[0], q[0], jumps[1], counts[1])
        _, sf = update(jumps[1], mt, jumps[2], counts[2])
        old = p["plastic_slip_m"][0]
        maps.append([old + sc, old + sm, (old + sm) + sf])
    assert (
        abs(native_actual[0] - native_actual[-1] - attempt["slip_errors_m"][0][0])
        < 1e-20
    )
    rows.append(
        dict(
            parent_progress=c["from"],
            target=c["to"],
            parent=p,
            counts=counts,
            jumps_m=paths,
            native_actual=native_actual,
            native_residual_in_independent_operator=[rc, rm, rf],
            independent_update_m=maps,
            residuals_native=[r["residual"] for r in triple],
            residuals_independent=[r["residual"] for r in ind],
        )
    )
res = subprocess.run(
    [binary],
    input="".join(json.dumps(r) + "\n" for r in requests),
    text=True,
    stdout=subprocess.PIPE,
    check=True,
)
values = iter(map(json.loads, res.stdout.splitlines()))
for row in rows:
    native = [next(values) for _ in range(3)]
    row["native_update_m"] = native
    np.testing.assert_allclose(native[0], row["native_actual"], rtol=0, atol=3e-18)
    # Rows: native/independent/one-Newton-corrected native jump path. Columns: native/independent law.
    e = np.array(
        [
            [x[0] - x[2] for x in pair]
            for pair in zip(native, row["independent_update_m"])
        ]
    )
    row["coarse_minus_fine_m"] = e.tolist()
    row["equilibrium_contribution_m"] = float(e[0, 0] - e[1, 0])
    row["update_contribution_m"] = float(e[1, 0] - e[1, 1])
    row["interaction_m"] = float((e[0, 0] - e[0, 1]) - (e[1, 0] - e[1, 1]))
assert next(values, None) is None
errors = np.array([r["coarse_minus_fine_m"] for r in rows])
budget = 1e-4 / 1024**2
summary = dict(
    paths=["native", "independent", "corrected_native"],
    updates=["native", "independent"],
    abs_error_min_median_max_m=np.quantile(abs(errors), [0, 0.5, 1], axis=0).tolist(),
    over_budget=np.sum(abs(errors) > budget, axis=0).tolist(),
    max_fixed_jump_update_delta_m=float(abs(errors[:, :, 0] - errors[:, :, 1]).max()),
    equilibrium_abs_min_median_max_m=np.quantile(
        abs(errors[:, 0, 0] - errors[:, 1, 0]), [0, 0.5, 1]
    ).tolist(),
)
output.with_name("c01-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
output.write_bytes(
    gzip.compress(
        (
            json.dumps(
                dict(
                    budget_m=1e-4 / 1024**2,
                    route_attempts=route["attempts"],
                    independent_fixed_solves=len(trials) * 3,
                    records=rows,
                ),
                allow_nan=False,
            )
            + "\n"
        ).encode(),
        mtime=0,
    )
)
