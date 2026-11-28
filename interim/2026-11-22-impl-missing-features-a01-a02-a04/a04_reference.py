"""Simultaneous saturated RT0 / elastic stick solve for the horizontal Tet4 patch."""
import numpy as np
from joint_reference import System
from joint_ext_reference import elastic_operators
from one_way_reference import constraint_space, geometry, fluid_load
from verify_high_order import physical

def reference(mesh, ground, case, dt, steps):
    xyz = np.array(mesh["nodes"])
    cells = np.array(mesh["tetrahedra"])
    mechanics = case["mechanics"]
    joint = System(mesh, mechanics["joint_system"])
    basis, offset = constraint_space(mesh, mechanics, ground)
    flow = {"head_datum_m": 0.0, **case["flow"]}
    case = {**case, "flow": flow}
    if (len(cells) != 12 or mesh.get("tetrahedra10")
            or any(mat["volume_source_per_s"] != 0. or
                   mat["retention"]["saturated_specific_storage_per_m"] != 0.
                   for mat in flow["materials"])
            or len(flow["analysis"]["periods"]) != 1
            or case["bishop"]["kind"] != "effective_saturation"):
        raise ValueError("reference requires the saturated, zero-source, zero bulk storage Tet4 patch")
    gamma = ground["initialization"]["water_unit_weight_kn_m3"]
    materials = {m["id"]:m for m in flow["materials"]}
    # Each coincident triangle of this continuous, matching planar column
    # has a single head. Mechanical nodes on opposite sides stay distinct.
    def key(nodes):
        return tuple(sorted(tuple(p) for p in xyz[nodes]))
    face_ids = {}
    element_faces = []
    inverse_mass = []
    for nodes, mid in zip(cells,mesh["material_ids"],strict=True):
        x=xyz[nodes]
        v=abs(np.linalg.det(x[1:]-x[0]))/6
        d=x-x.mean(axis=0)
        invk=np.linalg.inv(materials[mid]["saturated_conductivity_m_s"])
        gram=d @ invk @ d.T
        mass=(gram+np.trace(gram)/20)/(9*v)
        inverse_mass.append(np.linalg.inv(mass))
        element_faces.append([face_ids.setdefault(key(np.delete(nodes,j)),len(face_ids)) for j in range(4)])
    nf,nc,nd=len(face_ids),len(cells),xyz.size
    zface=np.array([np.mean(np.array(k)[:,2]) for k in face_ids])
    G=np.zeros((nf,nd))
    for p,c,n,w in zip(joint.points,joint.c,joint.normal,joint.weight,strict=True):
        G[face_ids[key(p["nodes_a"])]]+=w*(n @ c)
    C=joint.kt[:,None,None]*np.eye(3)+(joint.kn-joint.kt)[:,None,None]*joint.normal[:,:,None]*joint.normal[:,None,:]
    K=np.einsum("qai,qab,qbj,q->ij",joint.c,C,joint.c,joint.weight)
    seed=np.zeros(nd)
    geo=ground["initialization"]["geostatic"]
    bodies={r["id"]:r["body"] for r in mesh["regions"]}
    for e,(qp,_,ids,ke) in enumerate(elastic_operators(mesh,ground)):
        stress=geo["seed"]["states"][bodies[mesh["region_ids"][e]]]["effective_stress_compression_positive_kpa"]
        seed[ids]-=np.einsum("qai,a,q->i",qp["b"],stress,qp["weight"])
        K[np.ix_(ids,ids)]+=ke
    grid=[geometry(mesh,e,independent=True) for e in range(nc)]
    base,_=physical(mesh,ground,np.ones(nc,bool),ground["initialization"]["surface_surcharge_kpa"],
                    mechanics["pressure_on_original_surface_only"],joint_system=joint)
    initial=flow["initial_head"]["head_m"]
    h0=np.full(nc+nf,initial)
    load=fluid_load(mesh,ground,case,{"relative_cell_head_m":h0[:nc]-flow["head_datum_m"]},[],grid,base)["equivalent_mechanical_rhs_kn"]
    load+=gamma*(G.T @ (h0[nc:]-zface))
    B=np.zeros((nd,nc+nf))
    for e,op in enumerate(grid):
        ids=(3*op["nodes"][:,None]+np.arange(3)).ravel()
        B[ids,e]=gamma*np.einsum("qai,q->ai",op["gradient"],op["weight"]).ravel()
    B[:,nc:]=gamma*G.T
    kr=basis.T @ K @ basis
    rhs=basis.T @ (load-seed-K @ offset)
    x0=np.linalg.solve(kr,rhs)
    # Unknowns are free displacements, cell heads and all unique face heads.
    nh=nc+nf
    water=np.zeros((nh,nh))
    for e,(faces,A) in enumerate(zip(element_faces,inverse_mass,strict=True)):
        ids=nc+np.asarray(faces)
        water[e,e]=A.sum()
        water[e,ids]=-A.sum(axis=0)
        water[ids,e]+=A.sum(axis=1)
        water[np.ix_(ids,ids)]-=A
    storage=np.zeros((nh,basis.shape[1]))
    if case.get("joint_aperture_storage_m") is not None:
        storage[nc:]=G @ basis/dt
    system=np.block([[kr,-basis.T @ B],[-storage,water]])
    load_heads=np.zeros(nh)
    for bc in flow["analysis"]["periods"][0]["boundaries"]:
        value=bc["condition"]["value"]["head_m"]
        for f in mesh["face_sets"][bc["face_set"]]:
            fi=nc+face_ids[key(f)]
            row=len(x0)+fi
            system[row]=0.
            system[row,row]=1.
            load_heads[fi]=value-initial
    previous=x0
    result=[]
    for step in range(steps):
        full=np.linalg.solve(system,np.r_[rhs,load_heads-storage @ previous])
        x=full[:len(x0)]
        heads=h0+full[len(x0):]
        u=offset+basis @ x
        q=np.array([A @ (heads[e]-heads[nc+np.asarray(f)]) for e,(A,f) in enumerate(zip(inverse_mass,element_faces,strict=True))])
        traction=np.einsum("qab,qb->qa",C,joint.c @ u)
        normal=np.einsum("qa,qa->q",traction,joint.normal)
        shear=traction-normal[:,None]*joint.normal
        assert np.all(normal<0.) and np.all(np.linalg.norm(shear,axis=1)<-joint.mu*normal)
        assert np.min(heads)>xyz[:,2].max()
        result.append(dict(time=(step+1)*dt,u=u.reshape(-1,3),heads=heads[:nc],flux=q,
                           face_heads={k:heads[nc+i] for k,i in face_ids.items()},
                           volume_change={k:float(change) for k,change in
                                         zip(face_ids,G @ (basis @ (x-previous)),strict=True)}))
        previous=x
    return result
