"""Excavate a core in a plane-strain sector using native Tet4 mechanics."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import subprocess

import numpy as np

HERE = Path(__file__).resolve().parent
SPEC = json.loads((HERE/"criteria.json").read_text())
CFG = SPEC["kirsch"]
A = CFG["radius_m"]
R = CFG["outer_radius_m"]
K, G = CFG["bulk_kpa"], CFG["shear_kpa"]
NU = (3*K-2*G)/(2*(3*K+G))
S0 = -np.diag(CFG["horizontal_vertical_compression_kpa"])
FACES = ((1,2,3), (0,3,2), (0,1,3), (0,2,1))
# Six-point (degree-four) triangle quadrature for the exterior traction.
QUAD = []
for repeated, w in [(0.445948490915965, 0.223381589678011),
                    (0.091576213509771, 0.109951743655322)]:
    for i in range(3):
        l = np.full(3, repeated)
        l[i] = 1-2*repeated
        QUAD.append((l,w))


def solution(x):
    """Tension-positive Cartesian stress and displacement change from the intact state."""
    r = np.linalg.norm(x[:2])
    er = x[:2]/r
    et = np.array([-er[1],er[0]])
    c, s = er[0]**2-er[1]**2, 2*er[0]*er[1]
    m = (S0[0,0]+S0[1,1])/2
    d = (S0[0,0]-S0[1,1])/2
    t = (A/r)**2
    rr = m*(1-t)+d*(1-4*t+3*t*t)*c
    tt = m*(1+t)-d*(1+3*t*t)*c
    rt = -d*(1+2*t-3*t*t)*s
    basis = np.column_stack((er,et))
    stress = basis @ np.array([[rr,rt],[rt,tt]]) @ basis.T
    ur = (m+d*(4*(1-NU)-t)*c)*A*A/(2*G*r)
    ut = -d*(2*(1-2*NU)+t)*s*A*A/(2*G*r)
    return stress, ur*er+ut*et


def geometry(n, *, divisions=None):
    nt, nr, nz = divisions if divisions is not None else (2*n, n, n)
    rings = np.geomspace(A,R,nr+1)
    # Ring vertices precede core vertices, keeping retained node IDs across stages.
    coords = []
    node = {}
    for k in range(nz+1):
        for j,r in enumerate(rings):
            for i,theta in enumerate(np.linspace(0,np.pi/2,nt+1)):
                node[k,j,i] = len(coords)
                coords.append([r*np.cos(theta) if i<nt else 0.,
                               r*np.sin(theta) if i else 0.,k/nz])
    retained = len(coords)
    for k in range(nz+1):
        node[k,-1,0] = len(coords)
        coords.append([0.,0.,k/nz])
    xyz = np.array(coords)
    ring_cells, core_cells = [], []
    for k in range(nz):
        for i in range(nt):
            triangles = [([(-1,0),(0,i),(0,i+1)],core_cells)]
            for j in range(nr):
                a,b,c,d = (j,i),(j+1,i),(j+1,i+1),(j,i+1)
                triangles.extend([([a,b,c],ring_cells),([a,c,d],ring_cells)])
            for points, cells in triangles:
                p = sorted(points)
                a,b,c = [node[k,*v] for v in p]
                x,y,z = [node[k+1,*v] for v in p]
                for tet in ([a,b,c,z],[a,b,y,z],[a,x,y,z]):
                    if np.linalg.det(xyz[tet[1:]]-xyz[tet[0]]) < 0:
                        tet[1],tet[2] = tet[2],tet[1]
                    cells.append(tet)
    return xyz, retained, ring_cells, core_cells


def boundary(tets):
    faces = {}
    for e,tet in enumerate(tets):
        for ids in FACES:
            f = tuple(tet[i] for i in ids)
            key = tuple(sorted(f))
            if key in faces:
                del faces[key]
            else:
                faces[key] = (e,list(f))
    return faces


def traction(xyz, faces, radius, evaluate):
    """Consistent face nodal forces, using the actual outward facet normal."""
    load = np.zeros_like(xyz)
    for key,(_,f) in faces.items():
        points = xyz[f]
        if not np.all(np.isclose(np.linalg.norm(points[:,:2],axis=1),radius)):
            continue
        area = np.cross(points[1]-points[0],points[2]-points[0])/2
        if np.linalg.norm(area[:2])==0:
            continue
        for l,w in QUAD:
            force = evaluate(key,l@points) @ area[:2]
            load[f,:2] += np.outer(l,w*force)
    return load


def solve(binary, folder, xyz, cells, load):
    folder.mkdir()
    constraints, sets = [], {}
    for axis, name in enumerate(("sym_x","sym_y","plane_strain")):
        sets[name] = (list(range(len(xyz))) if axis==2
                      else np.where(xyz[:,axis]==0)[0].tolist())
        constraints.append(dict(node_set=name, direction=np.eye(3)[axis].tolist(),value_m=0))
    nodal = []
    for i,force in enumerate(load):
        if np.any(force):
            name = f"traction_{i}"
            sets[name] = [i]
            nodal.append(dict(node_set=name,total_force_kn=force.tolist()))
    mesh = dict(schema="soil-fem-mesh/1",units="m-kN-kPa",axis="z_up",crs="LOCAL",
                origin_m=[0,0,0],nodes=xyz.tolist(),tetrahedra=cells,
                material_ids=[1]*len(cells),node_sets=sets)
    case = dict(schema="soil-fem-mechanics/1",units="m-kN-kPa",
                materials=[dict(id=1,young_kpa=2*G*(1+NU),poisson=NU,body_force_kn_m3=[0,0,0])],
                boundary=dict(constraints=constraints,nodal_loads=nodal),
                relative_linear_tolerance=SPEC["linear_relative_max"],maximum_free_dofs=2000,
                linear_solver=dict(method="dense_lu"))
    paths = [folder/name for name in ("mesh.json","case.json","result.json")]
    for p,data in zip(paths,(mesh,case)):
        p.write_text(json.dumps(data),encoding="utf-8")
    subprocess.run([str(binary),"solve-mechanics",*map(str,paths)],check=True)
    return json.loads(paths[-1].read_bytes())


def tensor(s):
    return np.array([[s[0],s[3]],[s[3],s[1]]])


def balances(result):
    return dict(free_relative=result["true_relative_free_force_residual"],
                linear_relative=result["linear"]["true_relative_residual"],
                applied_kn=np.sum(result["applied_nodal_load_kn"],axis=0).tolist(),
                reaction_kn=np.sum(result["reaction_kn"],axis=0).tolist(),
                force_error_kn=result["force_balance_error_kn"],
                moment_error_kn_m=result["moment_balance_error_kn_m"])


def study(binary, folder, n, *, divisions=None):
    xyz,count,ring,core = geometry(n, divisions=divisions)
    full = ring+core
    initial = solve(binary,folder/"intact",xyz,full,traction(xyz,boundary(full),R,lambda k,x:S0))
    cells = {p["element_index"]:tensor(s) for p,s in zip(initial["integration_points"],
                                                         initial["stress_tension_positive_kpa"])}
    core_side = boundary(core)
    outer = xyz[:count]
    exposed = boundary(ring)
    # On shared facets n_retained=-n_core. Handoff uses computed core stress.
    held = traction(outer,exposed,A,lambda key,x:cells[len(ring)+core_side[key][0]])
    u0 = np.array(initial["displacement_m"])[:count,:2]
    report = {"retained_tets":len(ring),"removed_tets":len(core),"initial":balances(initial),
              "held_interface_kn":held.sum(axis=0).tolist(),"stages":[]}
    radius = np.linalg.norm(outer[:,:2],axis=1)
    wall = np.isclose(radius,A)
    du_exact = np.array([solution(x)[1] for x in outer])
    for j,factor in enumerate(CFG["release"]):
        outer_load = traction(outer,exposed,R,lambda key,x:(1-factor)*S0+factor*solution(x)[0])
        result = solve(binary,folder/f"release-{j}",outer,ring,outer_load+(1-factor)*held)
        du = np.array(result["displacement_m"])[:,:2]-u0
        profile = []
        wall_values = defaultdict(list)
        for p,s in zip(result["integration_points"],result["stress_tension_positive_kpa"]):
            x = np.array(p["position_m"])
            r = np.linalg.norm(x[:2])
            actual = tensor(s)
            if A <= r <= 2.5*A:
                exact = (1-factor)*S0+factor*solution(x)[0]
                profile.append(dict(position_m=x.tolist(),weight=p["weight_m3"],
                                    stress_kpa=actual.tolist(), reference_kpa=exact.tolist(),
                                    relative_to_30MPa=float(np.linalg.norm(actual-exact)/30000)))
            # Volume-weighted recovery at the two wall-axis nodes on each layer.
            indices = ring[p["element_index"]]
            for axis in (0,1):
                if any(wall[i] and outer[i,1-axis]==0 for i in indices):
                    wall_values[axis].append((p["weight_m3"],-actual[1-axis,1-axis]))
        wall_hoop = [float(np.average(np.array(wall_values[k])[:,1],
                                      weights=np.array(wall_values[k])[:,0])) for k in (0,1)]
        row = dict(release=factor,**balances(result),wall_hoop_recovered_kpa=wall_hoop,
                   wall_hoop_exact_kpa=((1-factor)*(-np.diag(S0)[::-1])
                                        +factor*np.array(CFG["wall_hoop_compression_kpa_horizontal_vertical"])).tolist(),
                   wall_axis_displacement_m=[float(du[wall & (outer[:,1-k]==0),k].mean()) for k in (0,1)],
                   wall_axis_reference_m=[float(factor*du_exact[wall & (outer[:,1-k]==0),k].mean()) for k in (0,1)],
                   outer_load_kn=outer_load.sum(axis=0).tolist(),
                   interface_kn=((1-factor)*held.sum(axis=0)).tolist(),
                   max_displacement_change_m=float(np.max(np.linalg.norm(du,axis=1))),
                   wall_displacement_error_max_m=float(np.max(np.linalg.norm(
                       du[wall]-factor*du_exact[wall],axis=1))), near_field=profile)
        report["stages"].append(row)
    return report


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("solver",type=Path)
    parser.add_argument("output",type=Path)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False,parents=True)
    rows = []
    for n in (1,2):
        folder = args.output/str(12*n**3)
        folder.mkdir()
        rows.append(study(args.solver.resolve(),folder,n))
    (args.output/"kirsch.json").write_text(json.dumps(rows,indent=2,allow_nan=False)+"\n")
