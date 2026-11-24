"""Affine tetrahedra with complete degree-p Lagrange displacement spaces."""
from itertools import product
from types import SimpleNamespace
from pathlib import Path
import sys
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"tools"))
import fem_reference as fem
from hydrostatic_cut_volume_reference import wet_tetrahedra


def simplex_rule(dimension, order):
    x,w = np.polynomial.legendre.leggauss(order)
    x,w = (x+1)/2,w/2
    rows, weights = [], []
    for indices in product(range(order),repeat=dimension):
        rem, weight, bary = 1., 1., []
        for i in indices:
            weight *= w[i]*rem
            bary.append(rem*x[i])
            rem *= 1-x[i]
        rows.append([rem,*bary])
        weights.append(weight)
    return np.array(rows), np.array(weights)


def lagrange(bary, alpha, degree):
    """Value and barycentric gradient, also defined when a factor vanishes."""
    q, dim = bary.shape
    value = np.ones((q,len(alpha)))
    grad = np.zeros((q,len(alpha),dim))
    for j in range(dim):
        for k in range(degree):
            take = alpha[:,j] > k
            factor = (degree*bary[:,j]-k)/(k+1)
            grad[:,take,:] *= factor[:,None,None]
            grad[:,take,j] += value[:,take]*degree/(k+1)
            value[:,take] *= factor[:,None]
    return value,grad


def strain_matrix(gradient):
    b = np.zeros((len(gradient),6,3*gradient.shape[1]))
    gx,gy,gz = gradient.transpose(2,0,1)
    b[:,0,0::3], b[:,1,1::3], b[:,2,2::3] = gx,gy,gz
    b[:,3,0::3], b[:,3,1::3] = gy,gx
    b[:,4,1::3], b[:,4,2::3] = gz,gy
    b[:,5,0::3], b[:,5,2::3] = gz,gx
    return b


def geometry(mesh,ground,case,degree,face_order=6):
    vertices = np.array(mesh["nodes"])
    cells = np.array(mesh["tetrahedra"])
    alpha = np.array([a for a in product(range(degree+1),repeat=4) if sum(a)==degree])
    numbering, nodes, elements = {}, [], []
    for cell in cells:
        ids = []
        for a in alpha:
            key = tuple(sorted((int(v),int(c)) for v,c in zip(cell,a) if c))
            if key not in numbering:
                numbering[key] = len(nodes)
                nodes.append((a@vertices[cell])/degree)
            ids.append(numbering[key])
        elements.append(np.array(ids))
    nodes = np.array(nodes)
    nd = nodes.size
    stiff = np.zeros((nd,nd))
    initial = np.zeros(nd)
    load = np.zeros_like(nodes)
    constrained = np.full_like(nodes,np.nan)
    assert not case["boundary"].get("surface_loads") and not case["boundary"].get("nodal_loads")
    for eq in case["boundary"]["constraints"]:
        selection = set(mesh["node_sets"][eq["node_set"]])
        axis, = np.flatnonzero(eq["direction"])
        assert eq["direction"][axis]==1
        for key,idx in numbering.items():
            if all(v in selection for v,w in key):
                constrained[idx,axis] = eq["value_m"] + np.asarray(eq.get("gradient",[0,0,0]))@(nodes[idx]-eq.get("reference_point_m",[0,0,0]))
    iv = ground["initialization"]
    geo = iv["geostatic"]
    level, gamma = iv["water_level_m"],iv["water_unit_weight_kn_m3"]
    assert not geo.get("hydrostatic_face_sets") and iv["surface_surcharge_kpa"]==0
    mat_by_id = {m["id"]:m for m in ground["materials"]}
    region = {r["id"]:r["body"] for r in mesh["regions"]}
    volume_bary,weights = simplex_rule(3,5)
    shape,gr = lagrange(volume_bary,alpha,degree)
    assert np.max(abs(shape.sum(axis=1)-1))<1e-13
    grad_ref = gr[:,:,1:]-gr[:,:,:1]
    operators,faces = [],{}
    for ei,(cell,ids) in enumerate(zip(cells,elements,strict=True)):
        xyz = vertices[cell]
        jac = (xyz[1:]-xyz[0]).T
        inv = np.linalg.inv(jac)
        w = abs(np.linalg.det(jac))*weights
        b = strain_matrix(grad_ref@inv)
        dofs = (ids[:,None]*3+np.arange(3)).ravel()
        mat = mat_by_id[mesh["material_ids"][ei]]
        mc = mat["mohr_coulomb"]
        d = fem.elastic_tensor(mc["young_kpa"],mc["poisson"])
        s0 = np.asarray(geo["seed"]["states"][region[mesh["region_ids"][ei]]]["effective_stress_compression_positive_kpa"])
        stiff[np.ix_(dofs,dofs)] += np.einsum("qai,ab,qbj,q->ij",b,d,b,w,optimize=True)
        initial[dofs] -= np.einsum("qai,a,q->i",b,s0,w)
        load[ids,2] -= (w@shape)*mat["unit_weight_kn_m3"]
        for tet in wet_tetrahedra(xyz,level):
            pos = volume_bary@tet
            local = (pos-xyz[0])@inv.T
            _,dg = lagrange(np.column_stack((1-local.sum(axis=1),local)),alpha,degree)
            gradient = (dg[:,:,1:]-dg[:,:,:1])@inv
            pressure = np.maximum(level-pos[:,2],0)*gamma
            dw = abs(np.linalg.det((tet[1:]-tet[0]).T))*weights
            initial[dofs] -= np.einsum("qai,q,q->ai",gradient,pressure,dw).ravel()
        operators.append((b,dofs,d,s0,mc))
        for face in ((1,2,3),(0,3,2),(0,1,3),(0,2,1)):
            key = tuple(sorted(cell[list(face)]))
            faces[key] = ei,face
    fb,fw = simplex_rule(2,face_order)
    def trace(face_ids, bary):
        ei,loc = faces[tuple(sorted(face_ids))]
        # Orient the integration points by the supplied triangle, independently
        # of the parent tetrahedron's local face numbering.
        l = np.zeros((len(bary),4))
        for j,v in enumerate(face_ids):
            l[:,list(cells[ei]).index(v)] = bary[:,j]
        val = lagrange(l,alpha,degree)[0]
        xyz = vertices[cells[ei][list(loc)]]
        normal = np.cross(xyz[1]-xyz[0],xyz[2]-xyz[0])
        return elements[ei],val,normal
    for surf in geo["base_loads"]["surface_loads"]:
        assert surf["kind"]=="pressure"
        for face in mesh["face_sets"][surf["face_set"]]:
            ids,phi,n = trace(face,fb)
            load[ids] -= np.outer(fw@phi,n*surf["pressure_kpa"])
    joint, = case["joint_system"]["joints"]
    assert joint["water"]["kind"]=="continuous"
    jump,area,normals,positions,parents = [],[],[],[],[]
    for parent,left in enumerate(mesh["face_sets"][joint["side_a"]]):
        xyz = vertices[left]
        right, = [r for r in mesh["face_sets"][joint["side_b"]] if set(map(tuple,vertices[r]))==set(map(tuple,xyz))]
        right = [next(i for i in right if np.array_equal(vertices[i],x)) for x in xyz]
        a,pa,n = trace(left,fb)
        b,pb,nb = trace(right,fb)
        assert np.allclose(n,-nb)
        length = np.linalg.norm(n)
        coords = fb@xyz
        water = gamma*np.maximum(level-coords[:,2],0)
        load[a] -= np.outer((fw*water)@pa,n)
        load[b] += np.outer((fw*water)@pb,n)
        row = np.zeros((len(fb),len(nodes)))
        row[:,a] -= pa
        row[:,b] += pb
        jump.extend(row)
        area.extend(fw*length)
        positions.extend(coords)
        normals.extend(np.tile(n/length,(len(fw),1)))
        parents.extend([parent]*len(fw))
    area = np.asarray(area)
    assert np.isclose(sum(area),1.) and set(parents)=={0,1}
    law = joint["material"]
    system = SimpleNamespace(normal=np.asarray(normals),
        kn=np.full(len(area),law["normal_stiffness_kpa_m"]),kt=np.full(len(area),law["tangential_stiffness_kpa_m"]),
        mu=np.full(len(area),law["friction_coefficient"]))
    return SimpleNamespace(nodes=nodes,k=stiff,initial=initial,load=load.ravel(),constraint=constrained.ravel(),
        scalar_c=np.asarray(jump),system=system,area=area,positions=np.asarray(positions),parents=np.asarray(parents),operators=operators)
