"""Dry two-sided patch from the 12 Tet10 geometry; independent full-node audit."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import fem_reference as ref
from plate6_offset_reference import Element


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def inputs():
    mesh = json.loads((ROOT / "examples/plate6-offset-mesh.json").read_bytes())
    xyz = np.array(mesh["nodes"])
    cells = np.array(mesh["tetrahedra10"])
    upper = np.array(mesh["region_ids"]) == 2
    shared = np.intersect1d(cells[upper], cells[~upper])
    duplicates = {int(n): len(xyz) + i for i, n in enumerate(shared)}
    coords = np.vstack((xyz, xyz[shared]))
    for i in np.flatnonzero(upper):
        cells[i] = [duplicates.get(int(n), int(n)) for n in cells[i]]
    coords[np.unique(cells[upper]), 2] += 0.1
    mesh.update(nodes=coords.tolist(), tetrahedra10=cells.tolist(),
                tetrahedra=cells[:, :4].tolist(), node_sets={}, face_sets={})
    for bound in mesh["regions"][1]["bounds_m"]:
        bound[2] += 0.1
    owners = ref.topology(mesh)
    key = next(key for key in owners if set(key) <= set(shared))
    cell_a = owners[key][0]
    nodes_a = ref.face_nodes(mesh, cell_a, key).tolist()
    if np.cross(coords[nodes_a[1]] - coords[nodes_a[0]], coords[nodes_a[2]] - coords[nodes_a[0]])[2] < 0:
        nodes_a = [nodes_a[j] for j in (0, 2, 1, 5, 4, 3)]
    nodes_b = [duplicates[n] for n in nodes_a]
    cell_b = owners[tuple(sorted(nodes_b[:3]))][0]
    section = dict(id="liner",young_kpa=200000.,poisson=0.25,thickness_m=0.1,
                   unit_weight_kn_m3=2.,drilling_factor=0.001)
    contact = dict(id="B",soil_b=dict(cell=cell_b,nodes=nodes_b),
                   material=dict(normal_stiffness_kpa_m=1000.,tangential_stiffness_kpa_m=500.,friction_coefficient=0.3))
    plate = dict(id="P",section=section,soil_a=dict(cell=cell_a,nodes=nodes_a),offset_m=-0.05,
                 exposed_side=1,tie_at_physical_surface=True,contacts=[contact])
    case = dict(plate=plate,soil={str(i):dict(young_kpa=12000.*i,poisson=0.25) for i in (1,2)},steps=[])
    n = coords.size + 18
    u0 = np.zeros(n)
    u0[:coords.size:3] = 0.0003 * coords[:, 0]
    u0[coords.size+1::3] = 0.0007 + 0.0002 * coords[nodes_a, 0]
    def load(compression=0., shear=0., flex=0., free=False):
        d = u0.copy()
        d[:coords.size:3] += flex*coords[:, 0]
        d[coords.size+1::3] += flex*coords[nodes_a, 0]
        upper_ids = np.unique(cells[upper])
        d[3*upper_ids] += shear + 0.1*flex*coords[upper_ids,0]
        d[3*upper_ids+2] -= compression
        force = np.zeros(n)
        released = [3*nodes_a[3],3*nodes_b[3]] if free else []
        if released:
            force[released] = 0.008
        return dict(event="load",prescribed_configuration=d.tolist(),free_dofs=released,
                    nodal_load=force.tolist(),plate_traction_kpa=[flex*10,0.,0.])
    install = dict(event="install", contacts={"B":dict(policy="stress_free")}, holder="installation supports")
    # Independent constant-force area projection. At a midside integral N_i=A/3.
    traction = np.array([1.,0.,-10.])
    actions = np.zeros(54)
    for i in range(3,6):
        f = traction/6
        actions[6*i:6*i+3] = -f
        actions[6*i+3:6*i+6] = -np.cross([0,0,0.1],f)
        actions[36+3*i:39+3*i] = f
    preload = dict(event="install", contacts={"B":dict(policy="transfer_force",supplier="jack",
                   traction_kpa=[traction.tolist() for _ in range(16)],requested_nodal_actions=actions.tolist())},holder="jack supports")
    case["steps"] = [load(),install,load(0.01),load(0.01,0.0005,0.0002),
                     load(0.01,0.0005,0.0002,True),load(0.01),load(-0.002),
                     load(0.01),dict(event="remove",holder="removal supports"),preload,load(0.012)]
    return mesh,case


def independent(mesh,case):
    """Separate GL125 soil and GL25 contact integration; no production fields."""
    xyz = np.array(mesh["nodes"])
    plate = case["plate"]
    a = np.array(plate["soil_a"]["nodes"])
    b = np.array(plate["contacts"][0]["soil_b"]["nodes"])
    nd = xyz.size + 18
    kernel = Element(xyz[a],plate["section"],plate["offset_m"])
    rotations = np.arange(xyz.size,nd).reshape(6,3)
    axes = np.arange(3)
    plate_rows = np.column_stack([3*a[:,None]+axes,rotations]).ravel()
    ksoil = np.zeros((nd,nd))
    stress_ops = []
    for ei, mid in enumerate(mesh["material_ids"]):
        soil = case["soil"][str(mid)]
        c = ref.elastic_tensor(soil["young_kpa"],soil["poisson"])
        gp = ref.volume(mesh,ei,independent=True)
        rows = (3*gp["nodes"][:,None]+axes).ravel()
        ksoil[np.ix_(rows,rows)] += np.einsum("qai,ab,qbj,q->ij",gp["b"],c,gp["b"],gp["weight"],optimize=True)
        stress_ops.append((rows,ref.volume(mesh,ei)["b"],c))
    # Explicit plate-side translations are retained. Six independent node ties
    # impose equality to A via multipliers for the free equilibrium solve below.
    plate_k = np.zeros_like(ksoil)
    plate_k[np.ix_(plate_rows,plate_rows)] = kernel.stiffness
    gp_jump = []
    weights = []
    for l,w in kernel.neutral.quadrature:
        shapes = kernel.neutral.shapes(*(l @ kernel.neutral.xy[:3]))[0]
        j = np.zeros((3,nd))
        for i,s in enumerate(shapes):
            j[:,3*b[i]:3*b[i]+3] += np.eye(3)*s
            j[:,3*a[i]:3*a[i]+3] -= np.eye(3)*s
            for axis in axes:
                j[:,rotations[i,axis]] -= s*np.cross(np.eye(3)[axis],[0.,0.,0.1])
        gp_jump.append(j)
        weights.append(w)
    j = np.array(gp_jump)
    w = np.array(weights)
    kc = np.einsum("gai,a,gaj,g->ij",j,[500.,500.,1000.],j,w,optimize=True)
    frames = []
    d = np.zeros(nd)
    previous_load = previous_reaction = np.zeros(nd)
    work = np.zeros(2)
    installation = None
    q0 = np.zeros(3)
    for i,step in enumerate(case["steps"],1):
        start = d.copy()
        if step["event"] == "install":
            installation = d.copy()
            q0 = np.array(step["contacts"]["B"].get("traction_kpa",[[0.,0.,0.]])[0])
        elif step["event"] == "remove":
            installation = None
        else:
            d = np.array(step["prescribed_configuration"])
        active = installation is not None
        traction = np.array(step.get("plate_traction_kpa",[0.,0.,0.])) if active else np.zeros(3)
        if active:
            traction[2] -= 0.2
        applied = np.array(step["nodal_load"]) if "nodal_load" in step else previous_nodal.copy()
        previous_nodal = applied.copy()
        applied[plate_rows] += kernel.load(traction)
        def actions(x):
            force = ksoil @ x
            q = np.zeros((len(w),3))
            if active:
                delta = x-installation
                force += plate_k @ delta
                q = q0 + (j @ delta) * [500.,500.,1000.]
                # The sequence unloads shear before opening and recloses in pure compression.
                q[q[:,2]>0] = 0.
                assert np.all(np.linalg.norm(q[:,:2],axis=1) <= -0.3*q[:,2]+1e-12)
                force += np.einsum("gai,ga,g->i",j,q,w)
            return force,q
        free = np.array(step.get("free_dofs",[]),dtype=int)
        if len(free):
            # Mixed system with all 18 plate-side translations independent;
            # prescribed soil rows remain known. Multiplier rows enforce ties.
            exposed = np.concatenate([np.arange(3*n,3*n+3) for n in a])
            split = np.zeros((nd,len(free)+18))
            tie = np.column_stack([np.zeros((18,len(free))),np.eye(18)])
            for k,row in enumerate(free):
                match = np.flatnonzero(exposed == row)
                if len(match):
                    tie[match[0],k] = -1.
                else:
                    split[row,k] = 1.
            split[exposed,len(free)+np.arange(18)] = 1.
            matrix = split.T @ (plate_k+kc) @ split
            matrix[:len(free),:len(free)] += ksoil[np.ix_(free,free)]
            saddle = np.block([[matrix,tie.T],[tie,np.zeros((18,18))]])
            force,_ = actions(d)
            soil_internal = ksoil @ d
            member_rhs = (applied-previous_nodal) - (force-soil_internal)
            rhs = np.r_[split.T @ member_rhs,np.zeros(18)]
            rhs[:len(free)] += (previous_nodal-soil_internal)[free]
            correction = np.linalg.solve(saddle,rhs)
            d[free] += correction[:len(free)]
        internal,q = actions(d)
        reaction = internal-applied
        if len(free):
            np.testing.assert_allclose(reaction[free],0,atol=1e-10)
            reaction[free]=0.
        work += [0.5*(previous_load+applied)@(d-start),0.5*(previous_reaction+reaction)@(d-start)]
        previous_load,previous_reaction=applied,reaction
        frames.append(dict(d=d.copy(),reaction=reaction,applied=applied,work=work.copy(),
                           instantaneous=(reaction+applied)@d,stress=[np.einsum("ij,qjk,k->qi",c,op,d[rows]) for rows,op,c in stress_ops],
                           q=q,sections=kernel.recover((d-installation)[plate_rows]) if active else None))
    return frames


def run(output,cli):
    output.mkdir(parents=True,exist_ok=False)
    mesh,case = inputs()
    save(output/"mesh.json",mesh)
    save(output/"case.json",case)
    expected = independent(mesh,case)
    command = [str(cli.resolve()),"offset-composite",str(output/"mesh.json"),str(output/"case.json")]
    dest = output/"run.json"
    subprocess.run([*command,str(dest)],check=True)
    native = json.loads(dest.read_bytes())
    attachment = case["plate"]["soil_a"]["nodes"]
    element = Element(np.array(mesh["nodes"])[attachment],case["plate"]["section"],case["plate"]["offset_m"])
    xy = [(l[1],l[2]) for l,_ in element.neutral.quadrature]
    rows = []
    for frame,target in zip(native["frames"][1:],expected,strict=True):
        def error(field,data):
            return float(np.max(np.abs(np.array(frame[field])-data)))
        row = dict(index=frame["accepted_index"],event=frame["event"],
                   displacement_max=error("configuration",target["d"]),reaction_max=error("reactions",target["reaction"]),
                   applied_max=error("applied",target["applied"]),soil_stress_max=error("soil_stress_kpa",target["stress"]),
                   instantaneous_error=error("instantaneous_external_kn_m",target["instantaneous"]),
                   cumulative_work_error=float(np.max(np.abs(np.array([frame["cumulative_applied_work_kn_m"],frame["cumulative_reaction_work_kn_m"]])-target["work"]))))
        if frame["active"]:
            contact = frame["active"]["contacts"][0]
            row["contact_generation"]=contact["generation"]
            row["local_nodal_work_difference"]=abs(contact["local_instantaneous_work_kn_m"]-contact["nodal_instantaneous_work_kn_m"])
            # The independent GL25 and native Jacobi16 are distinct physical points.
            # Compare traction fields by evaluating their quadratic traces at native points.
            xi = ref.jacobi(4,1)[0]
            eta = ref.jacobi(4,0)[0]
            pred=[]
            def polys(r,s):
                return [1,r,s,r*r,r*s,s*s]
            coeff=np.linalg.lstsq([polys(r,s) for r,s in xy],target["q"],rcond=None)[0]
            for r in xi:
                for s in eta:
                    pred.append(polys(r,(1-r)*s)@coeff)
            row["traction_max"]=float(np.max(np.abs(np.array([p["traction_kpa"] for p in contact["points"]])-pred)))
            samples=frame["active"]["last_sections"]["samples"]
            for key in ("membrane_force_kn_m","bending_moment_kn","transverse_shear_force_kn_m"):
                row[key+"_max"]=float(np.max(np.abs(np.array([s[key] for s in samples])-[s[key] for s in target["sections"]["samples"]])))
        rows.append(row)
    restarts=[]
    for stop in (8,9,10):
        cp=output/f"checkpoint-{stop}.json"
        resumed=output/f"resumed-{stop}.json"
        subprocess.run([*command,str(cp),"--until",str(stop)],check=True)
        subprocess.run([*command,str(resumed),"--resume",str(cp)],check=True)
        restarts.append(dict(index=stop,identical=resumed.read_bytes()==dest.read_bytes()))
    assert all(item["identical"] for item in restarts)
    refusals={}
    for field in ("arm","owner","reference_jump","initial_traction","recounted_work"):
        damaged=json.loads((output/"checkpoint-9.json").read_bytes())
        frozen=damaged["frames"][-1]["archived"][0]["contacts"][0]
        if field=="arm":
            frozen["arm_m"][2]+=0.01
        elif field=="owner":
            frozen["soil_b"]["cell"]=frozen["soil_a"]["cell"]
        elif field=="reference_jump":
            frozen["points"][0]["birth"]["jump_m"][0]+=0.001
        elif field=="initial_traction":
            frozen["points"][0]["birth"]["initial_traction_kpa"][2]-=1.
        else:
            damaged["frames"][-1]["cumulative_reaction_work_kn_m"]+=0.1
        changed=output/f"altered-{field}.json"
        save(changed,damaged)
        refused=output/f"invalid-run-{field}.json"
        result=subprocess.run([*command,str(refused),"--resume",str(changed)],capture_output=True,text=True)
        assert result.returncode==1 and "accepted path replay mismatch" in result.stderr and not refused.exists()
        refusals[field]=result.stderr.strip()
    summary=dict(tet10=len(mesh["tetrahedra"]),soil_nodes=len(mesh["nodes"]),frames=rows,
                 restarts=restarts,refusals=refusals,altered_archive_rejected=True,passed=False,b02_complete=False,b03_complete=False)
    save(output/"patch-results.json",summary)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("output",type=Path)
    parser.add_argument("--cli",type=Path,default=ROOT/"target/release"/("soil-fem-3d.exe" if os.name=="nt" else "soil-fem-3d"))
    args=parser.parse_args()
    run(args.output,args.cli)
