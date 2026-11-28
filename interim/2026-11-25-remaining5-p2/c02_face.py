"""Composite symmetric triangle integration for the coupled wet JOINT path."""
import argparse
from itertools import permutations, product
import inspect
import json
from math import factorial, prod
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
P1 = HERE.parent / "2026-11-25-remaining5-p1"
sys.path.insert(0, str(P1))
from c02_sensitivity import ref

# Proposed before the runs. Applied separately to each complete signed history.
PROPOSED_RELATIVE_WIDTH = 0.001


def triangle_rule(n):
    """Equal-area subtriangles, each with a positive symmetric degree-five rule."""
    points = [[1 / 3] * 3]
    weights = [0.225]
    for a, w in [(0.470142064105115, 0.132394152788506),
                 (0.101286507323456, 0.125939180544827)]:
        orbit = sorted(set(permutations((a, a, 1 - 2*a))))
        points.extend(orbit)
        weights.extend([w] * 3)
    tiles = []
    for i in range(n):
        for j in range(n-i):
            tiles.append([(i,j), (i+1,j), (i,j+1)])
            if i+j+1 < n:
                tiles.append([(i+1,j), (i+1,j+1), (i,j+1)])
    xy = (np.array(points) @ (np.array(tiles)/n)).reshape(-1, 2)
    bary = np.column_stack((1-xy.sum(axis=1), xy))
    w = np.tile(weights, len(tiles))/n**2
    for i in range(6):
        for j in range(6-i):
            exact = 2*factorial(i)*factorial(j)/factorial(i+j+2)
            np.testing.assert_allclose(w @ (xy[:,0]**i * xy[:,1]**j), exact, atol=3e-15, rtol=0)
    return bary, w


def surface(mesh, ground, case, degree, divisions):
    # Hold volume operators and external loads fixed. The wet seam is planar,
    # horizontal, with spatially constant water pressure on each original face.
    g = ref.geometry(mesh, ground, case, degree, 6)
    joint, = case["joint_system"]["joints"]
    xyz = np.array(mesh["nodes"])
    keys = {}
    for cell in mesh["tetrahedra"]:
        for a in product(range(degree+1), repeat=4):
            if sum(a) == degree:
                key = tuple(sorted((v,k) for v,k in zip(cell,a) if k))
                idx = keys.setdefault(key, len(keys))
                np.testing.assert_allclose(g.nodes[idx], np.array(a) @ xyz[cell]/degree, atol=2e-15)
    bary, weights = triangle_rule(divisions)
    rows, areas, coords, normals, parents = [], [], [], [], []
    lefts, rights = [mesh["face_sets"][joint[name]] for name in ("side_a", "side_b")]
    for i, left in enumerate(lefts):
        corners = xyz[left]
        assert np.ptp(corners[:,2]) == 0
        right, = [r for r in rights if set(map(tuple,xyz[r])) == set(map(tuple,corners))]
        right = [next(v for v in right if np.array_equal(xyz[v], p)) for p in corners]
        cross = np.cross(corners[1]-corners[0], corners[2]-corners[0])
        cell, = [c for c in mesh["tetrahedra"] if set(left)<=set(c)]
        fourth, = set(cell)-set(left)
        if cross @ (xyz[fourth]-corners[0]) > 0:
            cross = -cross
        area = np.linalg.norm(cross)/2
        block = np.zeros((len(bary),len(keys)))
        for a in product(range(degree+1), repeat=3):
            if sum(a) != degree:
                continue
            phi = prod((degree*bary[:,j]-k)/(k+1)
                       for j in range(3) for k in range(a[j]))
            for side, sign in ((left,-1), (right,1)):
                idx = keys[tuple(sorted((v,k) for v,k in zip(side,a) if k))]
                block[:,idx] = sign*phi
        np.testing.assert_allclose(block.sum(axis=1), 0, atol=3e-14)
        rows.extend(block)
        areas.extend(area*weights)
        normals.extend([cross/(2*area)]*len(bary))
        coords.extend(bary@corners)
        parents.extend([i]*len(bary))
    g.scalar_c = np.array(rows)
    g.area, g.positions, g.parents = map(np.asarray, (areas, coords, parents))
    law = joint["material"]
    g.system = SimpleNamespace(normal=np.array(normals),
        kn=np.full(len(areas), law["normal_stiffness_kpa_m"]),
        kt=np.full(len(areas), law["tangential_stiffness_kpa_m"]),
        mu=np.full(len(areas), law["friction_coefficient"]))
    np.testing.assert_allclose(sum(areas), 1., atol=2e-14, rtol=0)
    return g


def face_force(shape, traction, area):
    return np.einsum("qci,qc,q->i", shape, traction, area, optimize=True)


def face_tangent(shape, local, area):
    return np.einsum("qci,qcd,qdj,q->ij", shape, local, shape, area, optimize=True)


def coupled_solver():
    # Keep the archived Newton loop, volume condensation and acceptance checks.
    # Replace both weak-form contractions, including the final full force audit.
    source = inspect.getsource(ref.reference)
    edits = {
        "face_order,folder):": "face_order,folder,grids,observe):",
        "for step in range(1,steps+1):": "for t in grids[phase]:",
        "            t = step/steps\n": "",
        "if step%(steps//64)==0:": "observe(phase,t,g,u,status,width)\n            if t*64 == int(t*64):",
        "c.T@(q[:,0]*g.area[:,None]).ravel()": "face_force(ca,q[:,0],g.area)",
        "c.T@(q*g.area[:,None]).ravel()": "face_force(ca,q,g.area)",
        "c.T@(g.area[:,None,None]*(local@ca)).reshape(3*n,-1)": "face_tangent(ca,local,g.area)",
    }
    for before,after in edits.items():
        assert source.count(before)==1
        source=source.replace(before,after)
    namespace=dict(ref.__dict__,geometry=surface,face_force=face_force,face_tangent=face_tangent)
    exec(compile(source,ref.__file__,"exec"),namespace)
    return namespace["reference"]


def run(out, p, div, n):
    out.mkdir(parents=True)
    mesh,ground,case = [json.loads((ref.ROOT/"examples"/(s+".json")).read_bytes())
        for s in ("joint-water-column-tet10-mesh","ground-joint-water-mc","stages-joint-water-cycle")]
    solve = coupled_solver()
    common = ref.geometry(mesh,ground,case,p,6)
    old = np.zeros((len(common.area),3))
    traction, flow = old.copy(), old.copy()
    samples, states = [], []
    def observe(phase,t,g,u,status,width):
        nonlocal old, traction, flow
        target = common.scalar_c@u.reshape(-1,3)
        counts = np.maximum(1,np.ceil(np.linalg.norm(target-old,axis=1)/width)).astype(int)
        q, increment, st = ref.project(common.system,old,traction,target[:,None,:],counts)
        old, traction = target, q[:,0]
        flow += increment[:,0]
        if t*64 == int(t*64):
            samples.append(flow.copy())
            states.append(st[:,0].copy())
    grids = [np.arange(1,n+1)/n for _ in range(7)]
    cost = solve(mesh,ground,case,p,n,div,out,grids,observe)
    del cost["face_gauss_order"]
    cost.update(subdivisions_per_edge=div, rule="symmetric degree 5", proposal=PROPOSED_RELATIVE_WIDTH)
    with np.load(out/"slip.npz") as data:
        total = np.einsum("tqc,q->tc",data["slip"],data["weight"])
    assert total.shape==(448,3) and len(samples)==448
    np.savez_compressed(out/"principal.npz",slip=samples,status=states,weight=common.area,
                         parent=common.parents,position=common.positions,signed_area_integral=total)
    (out/"cost.json").write_text(json.dumps(cost,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out",type=Path)
    parser.add_argument("degree",type=int,choices=[3,4])
    parser.add_argument("divisions",type=int,choices=[2,4,8])
    parser.add_argument("steps",type=int,choices=[256,512])
    args=parser.parse_args()
    run(args.out,args.degree,args.divisions,args.steps)
