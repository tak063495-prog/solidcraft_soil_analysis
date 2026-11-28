"""Matched shifted surfaces and an embedded post-solve quadrature on every leaf."""
import inspect
import json
from math import fsum
from pathlib import Path
import sys

import numpy as np

from freeze import HERE, P7, save
sys.path.insert(0, str(P7))
import face_shift as p7

p6 = p7.p6


def geometry(faces):
    # Install child rules only in the passive geometry constructor.
    rules = {}
    bary, w = p6.face.triangle_rule(1)
    for face in faces:
        leaves = [child for tile in face["leaves"] for child in p6.split(np.array(tile))]
        tiles = np.asarray(leaves)
        xy = (bary @ tiles).reshape(-1, 2)
        weights = np.outer(abs(np.linalg.det(tiles[:,1:]-tiles[:,:1])), w).ravel()
        rules[face["parent"]] = np.column_stack((1-xy.sum(axis=1),xy)), weights
    source = inspect.getsource(p6.face.surface)
    source = source.replace("    bary, weights = triangle_rule(divisions)\n","")
    source = source.replace("        corners = xyz[left]", "        corners = xyz[left]\n        bary, weights = rules[i]")
    source = source.replace("assert_allclose(sum(areas)", "assert_allclose(fsum(areas)")
    namespace = dict(vars(p6.face), rules=rules, fsum=fsum)
    exec(compile(source,p6.face.__file__,"exec"),namespace)
    inputs = [json.loads((p6.face.ref.ROOT/"examples"/(name+".json")).read_text()) for name in
              ("joint-water-column-tet10-mesh","ground-joint-water-mc","stages-joint-water-cycle")]
    return namespace["surface"](*inputs,4,8)


def solve(out, sign):
    p7.SHIFT = sign*np.array([1/16,1/32])
    factory = p6.face.coupled_solver
    displacements, clock = [], []
    def instrument():
        original = factory()
        def solve_observed(mesh, ground, case, p, n, div, folder, grids, observe):
            def observation(phase, t, g, u, status, width):
                displacements.append(u.copy())
                clock.append((phase,t,width))
                observe(phase,t,g,u,status,width)
            return original(mesh,ground,case,p,n,div,folder,grids,observation)
        return solve_observed
    p6.face.coupled_solver = instrument
    p7.run(out)
    np.savez_compressed(out/"kinematics.npz",u=displacements,clock=clock)


def residual(out):
    faces = json.loads((out/"tiles.json").read_text())["faces"]
    g = geometry(faces)
    kinematics = np.load(out/"kinematics.npz")
    active = np.load(out/"slip.npz")
    weights = g.area.reshape(-1,28)
    parent_w = active["weight"].reshape(-1,7)
    np.testing.assert_allclose(weights.sum(axis=1), parent_w.sum(axis=1), atol=1e-16)
    cols = np.any(g.scalar_c != 0,axis=0)
    op = np.ascontiguousarray(g.scalar_c[:,cols])
    old = np.zeros((len(g.area),3))
    traction, flow = old.copy(), old.copy()
    estimates, mixed, changed, child_areas = [], [], [], []
    previous = None
    transitions = np.zeros(len(weights), dtype=bool)
    active_slip, active_status = active["slip"], active["status"]
    anchor = 0
    for (phase,t,width), u in zip(kinematics["clock"],kinematics["u"],strict=True):
        target = op @ u.reshape(-1,3)[cols]
        counts = np.maximum(1,np.ceil(np.linalg.norm(target-old,axis=1)/width)).astype(int)
        q, inc, labels = p6.face.ref.project(g.system,old,traction,target[:,None,:],counts)
        old, traction = target,q[:,0]
        flow += inc[:,0]
        st = labels[:,0].reshape(-1,28)
        transitions |= (st != st[:,:1]).any(axis=1)
        if previous is not None:
            transitions |= (st != previous).any(axis=1)
        previous = st.copy()
        if t*64 != int(t*64):
            continue
        a = active_status[anchor].reshape(-1,7)
        crossing = (st != st[:,:1]).any(axis=1) | (a != st[:,:1]).any(axis=1)
        coarse = np.einsum("lqc,lq->lc",active_slip[anchor].reshape(-1,7,3),parent_w)
        fine = np.einsum("lqc,lq->lc",flow.reshape(-1,28,3),weights)
        estimates.append(fine-coarse)
        child_areas.append(fine)
        mixed.append(crossing)
        changed.append(transitions.copy())
        transitions[:] = False
        anchor += 1
    assert anchor == 448
    np.savez_compressed(out/"residual.npz", signed=estimates, child=child_areas,
                        mixed=mixed, interval_switch=changed, area=weights.sum(axis=1),
                        centroid=np.average(g.positions.reshape(-1,28,3),axis=1,weights=weights[0]),
                        parent=g.parents.reshape(-1,28)[:,0])


if __name__ == "__main__":
    out = Path(sys.argv[2])
    if sys.argv[1] == "solve":
        solve(out, int(sys.argv[3]))
    elif sys.argv[1] == "residual":
        residual(out)
