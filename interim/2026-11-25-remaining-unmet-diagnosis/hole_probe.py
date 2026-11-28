"""Compare the same twelve interior locations across angular and radial refinements."""
import json
import sys
from pathlib import Path
import numpy as np

root=Path(sys.argv[1]).resolve()
sys.path.insert(0,str(root/"interim/2026-11-24-t9d-benchmarks"))
from kirsch import solution, study, tensor

out=Path(sys.argv[2])
out.mkdir()
solver=root/"target/release/soil-fem-3d"
rows=[]
for nt,nr in [(4,2),(8,2),(4,4),(8,4)]:
    folder=out/f"theta{nt}-radial{nr}"
    folder.mkdir()
    history=study(solver,folder,None,divisions=(nt,nr,1))
    read=lambda path:json.loads((folder/path).read_bytes())
    mesh=read("release-2/mesh.json")
    final=read("release-2/result.json")
    nodes=np.array(mesh["nodes"])
    cells=np.array(mesh["tetrahedra"])
    du=np.array(final["displacement_m"])[:,:2]-np.array(read("intact/result.json")["displacement_m"])[:len(nodes),:2]
    stress={p["element_index"]:tensor(s) for p,s in zip(final["integration_points"],final["stress_tension_positive_kpa"])}
    vertices=nodes[cells]
    transforms=np.linalg.inv((vertices[:,1:]-vertices[:,:1]).transpose(0,2,1))
    points=[]
    for radius in (1.1,1.5,2.25):
        for angle in (7,27,52,79):
            theta=np.deg2rad(angle)
            xyz=np.array([radius*np.cos(theta),radius*np.sin(theta),0.37])
            b=np.einsum("ijk,ik->ij",transforms,xyz-vertices[:,0])
            weights=np.column_stack((1-b.sum(axis=1),b))
            inside=np.flatnonzero(np.all(weights>0.,axis=1))
            assert len(inside)==1
            i=inside[0]
            sigma,disp=solution(xyz)
            interpolated=weights[i]@du[cells[i]]
            err_sigma=np.linalg.norm(stress[i]-sigma)/30000
            err_u_mm=1000*np.linalg.norm(interpolated-disp)
            points.append(dict(r=radius,angle=angle,position_m=xyz.tolist(),element=int(i),
                               barycentric=weights[i].tolist(),stress_kPa=stress[i].tolist(),
                               reference_stress_kPa=sigma.tolist(),displacement_m=interpolated.tolist(),
                               reference_displacement_m=disp.tolist(),
                               stress_norm_over_30MPa=err_sigma,u_error_mm=err_u_mm))
    rows.append(dict(theta_cells=nt,radial_cells=nr,tets=len(cells),points=points,
                     stress_rms_over_30MPa=float(np.sqrt(np.mean([p["stress_norm_over_30MPa"]**2 for p in points]))),
                     displacement_rms_mm=float(np.sqrt(np.mean([p["u_error_mm"]**2 for p in points]))),
                     max_free_force=max(r["free_relative"] for r in [history["initial"],*history["stages"]])))
(out/"summary.json").write_text(json.dumps(rows,indent=2,allow_nan=False)+"\n",encoding="utf-8")
