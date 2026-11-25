"""Prepare high-order bulk operators and audit every accepted native increment."""
import json
from pathlib import Path
import subprocess
import sys
import time
from types import FunctionType

import numpy as np
from scipy.linalg import cho_factor, cho_solve

HERE = Path(__file__).resolve().parent
P1 = HERE.parent / "2026-11-25-remaining5-p1"
sys.path.insert(0, str(P1))
from c02_sensitivity import ref
import c02_geometry


def geometry(inputs, degree, q):
    # Duffy-Legendre volume integration: degree 2p-2 plus the Duffy Jacobian.
    original = c02_geometry.simplex_rule
    def rule(dim, order):
        return original(dim, max(degree+1, order) if dim == 3 else order)
    build = FunctionType(ref.geometry.__code__, dict(ref.geometry.__globals__, simplex_rule=rule))
    return build(*inputs, degree, q)


def run(work, out, binary, degree, n, sub):
    start = time.monotonic()
    work.mkdir(parents=True)
    out.mkdir(parents=True)
    inputs = [json.loads((ref.ROOT / "examples" / (s + ".json")).read_bytes()) for s in
              ("joint-water-column-tet10-mesh", "ground-joint-water-mc", "stages-joint-water-cycle")]
    mesh, _, case = inputs
    g = geometry(inputs, degree, 10)
    common = geometry(inputs, degree, 6)
    free = np.isnan(g.constraint)
    fixed = np.nan_to_num(g.constraint)
    face = np.any(g.scalar_c != 0, axis=0)
    a = np.flatnonzero(np.repeat(face, 3))
    i = np.flatnonzero(free & ~np.repeat(face, 3))
    assert free[a].all() and not np.any(common.scalar_c[:, ~face])
    kai = g.k[np.ix_(a,i)]
    inside = cho_factor(g.k[np.ix_(i,i)])
    lift = cho_solve(inside, kai.T)
    schur = g.k[np.ix_(a,a)] - kai @ lift
    drive, boundary = g.load - g.initial, -g.k @ fixed
    def reduce(f):
        return (f[a] - kai @ cho_solve(inside, f[i])).tolist()
    switches = json.loads((P1 / "c02/summary.json").read_bytes())["refined_base_intervals"]
    stages = ["geostatic"] + [s["name"] for s in case["stages"]]
    targets = [0.] + [s["target_factor"] for s in case["stages"]]
    grids = []
    for name in stages:
        grid = []
        for k in range(1,n+1):
            base = (256*k+n-1)//n
            split = sub if name == "reverse" or base in switches[name] else 1
            grid.extend((k-1+j/split)/n for j in range(1,split+1))
        assert len(grid) <= (2065*len(mesh["tetrahedra"])+11)//12
        grids.append(grid)
    data = dict(elements=len(mesh["tetrahedra"]), stiffness=schur.tolist(), rhs_load=reduce(drive),
                rhs_boundary=reduce(boundary), load=g.load.tolist(), initial=g.initial.tolist(),
                trace=g.scalar_c[:,face].tolist(), observation=common.scalar_c[:,face].tolist(),
                normal=g.system.normal.tolist(), observation_normal=common.system.normal.tolist(),
                material=case["joint_system"]["joints"][0]["material"], area=g.area.tolist(),
                targets=targets, grids=grids)
    model, raw = work / "operators.json", work / "steps.jsonl"
    model.write_text(json.dumps(data))
    subprocess.run([binary, "study-joint-p-cycle", model, raw], check=True)
    samples, total, steps = [], [], []
    inner_load, inner_bc = (cho_solve(inside,f[i]) for f in (drive,boundary))
    c = np.kron(g.scalar_c[:,face], np.eye(3))
    with raw.open() as f:
        for phase, grid in enumerate(grids):
            origin = targets[phase-1] if phase else 0.
            for t in grid:
                row = json.loads(next(f))
                assert (row["phase"], row["t"]) == (phase,t)
                x = np.asarray(row["x"])
                frac = t if phase == 0 else 1.
                factor = origin+(targets[phase]-origin)*t
                u = factor*fixed
                u[a] = x
                u[i] = frac*inner_load + factor*inner_bc - lift@x
                force = g.k@u - frac*drive
                force[a] += c.T@(np.asarray(row["traction"])*g.area[:,None]).ravel()
                rel = np.linalg.norm(force[free])/max(1.,np.linalg.norm(g.initial + frac*drive))
                y = max(float(np.max(ref.normalized_yield(mc,s0-(b@u[dofs])@d.T)))
                        for b,dofs,d,s0,mc in g.operators)
                assert rel < 1e-10 and y <= 1e-8, (phase,t,rel,y)
                steps.append([phase,t,rel,row["linear"],y,row["iterations"]])
                if t*64 == int(t*64):
                    samples.append(row["slip"])
                    total.append(row["integral"])
        assert not f.readline()
    assert len(samples) == len(total) == 448
    np.savez_compressed(out/"principal.npz", slip=samples, signed_area_integral=total, weight=common.area,
                        position=common.positions,parent=common.parents)
    np.savez_compressed(out/"increments.npz",phase_t_fullforce_linear_yield_iterations=steps)
    report = dict(degree=degree,base_steps=n,switch_reverse_split=sub,active_face_order=10,
                  volume_order=max(5,degree+1), elements=len(mesh["tetrahedra"]),dofs=g.nodes.size,
                  attempts_per_stage=list(map(len,grids)), accepted=len(steps),anchors=len(samples),
                  max_force=max(s[2] for s in steps),max_linear=max(s[3] for s in steps),
                  max_bulk_yield=max(s[4] for s in steps),seconds=time.monotonic()-start)
    (out/"cost.json").write_text(json.dumps(report,indent=2)+"\n")


if __name__ == "__main__":
    run(*map(Path,sys.argv[1:4]),*map(int,sys.argv[4:]))
