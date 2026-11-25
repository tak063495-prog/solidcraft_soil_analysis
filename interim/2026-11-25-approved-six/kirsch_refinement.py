"""Refine inner radial cells and score raw stresses and the two wall intercepts."""
import argparse
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np


def score(mesh, final, initial):
    from kirsch import A, solution, tensor
    xyz = np.array(mesh["nodes"])
    cells = np.array(mesh["tetrahedra"])
    vertices = xyz[cells]
    inverse = np.linalg.inv((vertices[:, 1:] - vertices[:, :1]).transpose(0, 2, 1))
    stress = np.empty((len(cells), 2, 2))
    for p, s in zip(final["integration_points"], final["stress_tension_positive_kpa"], strict=True):
        stress[p["element_index"]] = tensor(s)
    du = np.array(final["displacement_m"])[:,:2] - np.array(initial["displacement_m"])[:len(xyz),:2]

    def errors(position, computed, movement):
        stress_exact, u_exact = solution(position)
        er = position[:2] / np.linalg.norm(position[:2])
        frame = np.stack((er, [-er[1],er[0]]))
        d = frame @ (computed - stress_exact) @ frame.T / 30000
        return [d[0,0],d[1,1],d[0,1],np.linalg.norm(d),*(1000*frame @ (movement-u_exact))]

    field, weights = [], []
    for p in final["integration_points"]:
        x = np.array(p["position_m"])
        if A <= np.linalg.norm(x[:2]) <= 2.5*A:
            i = p["element_index"]
            field.append(errors(x,stress[i],du[cells[i]].mean(axis=0)))
            weights.append(p["weight_m3"])
    points = []
    for r in (1.1,1.5,2.25):
        for theta in np.deg2rad([7,27,52,79]):
            x = np.array([r*np.cos(theta),r*np.sin(theta),0.37])
            natural = np.einsum("eij,ej->ei",inverse,x-vertices[:,0])
            bary = np.column_stack((1-natural.sum(axis=1),natural))
            inside = np.flatnonzero((bary>0).all(axis=1))
            if len(inside)!=1:
                raise ValueError(f"sample is not inside one Tet4 at {x}")
            i = inside[0]
            points.append(errors(x,stress[i],bary[i]@du[cells[i]]))
    def summarize(rows, weights=None):
        a=np.array(rows)
        names=("rr_over_30MPa","hoop_over_30MPa","shear_over_30MPa","tensor_over_30MPa","radial_mm","tangential_mm")
        return dict(count=len(rows),rms=dict(zip(names,np.sqrt(np.average(a*a,axis=0,weights=weights)),strict=True)),
                    max_abs=dict(zip(names,np.max(np.abs(a),axis=0),strict=True)))
    return dict(points=summarize(points),band=summarize(field,weights))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("root",type=Path)
    p.add_argument("output",type=Path)
    p.add_argument("--theta",type=int)
    grid = p.add_mutually_exclusive_group(required=True)
    grid.add_argument("--inner-splits",type=int)
    grid.add_argument("--band-cells",type=int,help="cells between a and 2.5a out of 24 radial cells")
    grid.add_argument("--registered-grid",choices=["inverse-radius-22x21"])
    args=p.parse_args()
    if args.registered_grid:
        if args.theta not in (None,22):
            p.error("inverse-radius-22x21 fixes theta at 22")
        args.theta=22
    elif args.theta is None:
        p.error("--theta is required with --inner-splits or --band-cells")
    root=args.root.resolve()
    sys.path[:0]=[str(root/"interim/2026-11-24-approved-four"),str(root/"interim/2026-11-24-t9d-benchmarks")]
    from kirsch import A,R,study
    from kirsch_wall import extrapolate
    start=time.monotonic()
    folder=args.output
    folder.mkdir()
    rings=np.geomspace(A,R,17)
    if args.registered_grid:
        refined=A/np.linspace(1,A/R,22)
        refined[[0,-1]]=[A,R]
    elif args.band_cells is None:
        splits=args.inner_splits
        if not 0<=splits<=16:
            p.error("inner splits outside the initial 16 radial cells")
        refined=np.sort(np.r_[rings,np.sqrt(rings[:splits]*rings[1:splits+1])])
    else:
        if not 0<args.band_cells<24:
            p.error("band cells must leave at least one cell outside the band")
        edge=2.5*A
        refined=np.r_[np.geomspace(A,edge,args.band_cells+1),
                      np.geomspace(edge,R,24-args.band_cells+1)[1:]]
    nt=args.theta
    nr=len(refined)-1
    if (6*nt*nr in (768,6144) or 4*nt*(nr+1)>2000):
        raise ValueError("mesh exceeds the study bounds")
    run=study(root/"target/release/soil-fem-3d",folder,None,divisions=(nt,nr,1),radii=None if args.inner_splits == 0 else refined)
    read=lambda path: json.loads((folder/path).read_bytes())
    final=read("release-2/result.json")
    field=score(read("release-2/mesh.json"),final,read("intact/result.json"))
    wall=[x["wall_kpa"] for x in extrapolate(final,nt,nr)]
    ref=run["stages"][-1]["wall_hoop_exact_kpa"]
    wall_error=np.array(wall)/ref-1
    result=dict(theta=nt,radial=nr,radii_m=refined.tolist(),elements=run["retained_tets"],solves=4,
                wall_kpa=wall,wall_relative_errors=wall_error.tolist(),wall_passed=bool(max(abs(wall_error))<=0.03),
                **field,max_free_relative=max(s["free_relative"] for s in [run["initial"],*run["stages"]]),
                max_linear_relative=max(s["linear_relative"] for s in [run["initial"],*run["stages"]]),
                seconds=time.monotonic()-start,python_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                solver_peak_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    print(json.dumps(result,indent=2,allow_nan=False))

if __name__=="__main__":
    main()
