"""Input-only p/time comparison for the wet opening and reversal cycle."""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from c02_geometry import ROOT, geometry
from mc_reference import normalized_yield


def project(law, old, traction, targets, counts):
    """Coulomb update on the portion of each straight subpath in closed contact."""
    values = np.zeros_like(targets)
    slip = np.zeros_like(targets)
    status = np.zeros(targets.shape[:2], dtype=np.uint8)
    for count in np.unique(counts):
        where = counts == count
        n = law.normal[where, None, :]
        kt = law.kt[where, None, None]
        kn = law.kn[where, None]
        mu = law.mu[where, None]
        start = old[where, None, :]
        delta = (targets[where] - start) / count
        dg = (delta*n).sum(axis=-1)
        tangential = delta - dg[...,None]*n
        gap0 = (start*n).sum(axis=-1)
        tau0 = traction[where, None, :]
        tau = tau0 - (tau0*n).sum(axis=-1,keepdims=True)*n
        flow = np.zeros_like(delta)
        for step in range(count):
            ga = gap0 + step*dg
            gb = gap0 + (step+1)*dg
            opening = (ga <= 0) & (gb > 0)
            closing = (ga > 0) & (gb <= 0)
            crossing = np.zeros_like(dg)
            np.divide(-ga, dg, out=crossing, where=opening|closing)
            closed_length = np.where(opening, crossing, np.where(closing, 1-crossing, gb <= 0))
            trial = tau + kt*tangential*closed_length[...,None]
            length = np.linalg.norm(trial, axis=-1)
            pressure = -kn*np.minimum(gb,0)
            radius = mu*pressure
            sliding = length > radius
            direction = np.divide(trial, length[...,None], out=np.zeros_like(trial),
                                  where=length[...,None]>0)
            correction = np.maximum(length-radius,0)[...,None]*direction
            flow += correction/kt
            tau = np.where(sliding[...,None],radius[...,None]*direction,trial)
        values[where] = tau - pressure[...,None]*n
        slip[where] = flow
        status[where] = np.where(gb>0,0,np.where(sliding,2,1))
    return values, slip, status


def reference(mesh,ground,case,degree,steps,face_order,folder):
    import resource  # Linux route profiling

    start = time.monotonic()
    g = geometry(mesh,ground,case,degree,face_order)
    n = len(g.area)
    free = np.isnan(g.constraint)
    fixed = np.nan_to_num(g.constraint)
    seam = np.repeat(np.any(g.scalar_c!=0,axis=0),3)
    a = np.flatnonzero(seam & free)
    i = np.flatnonzero(~seam & free)
    assert not (seam & ~free).any()
    c = np.kron(g.scalar_c,np.eye(3))[:,a]
    ca = c.reshape(n,3,-1)
    # Condense only the free interior coordinates. Interface displacements
    # remain Newton unknowns and receive the independently integrated tractions.
    kai = g.k[np.ix_(a,i)]
    inside = cho_factor(g.k[np.ix_(i,i)])
    lift = cho_solve(inside,kai.T)
    condensed = g.k[np.ix_(a,a)] - kai@lift
    gfixed = g.k@fixed
    material_steps = 0
    worst_residual = worst_linear = 0.
    worst_yield = -np.inf
    attempts = []
    old_jump = np.zeros((n,3))
    traction = old_jump.copy()
    plastic = old_jump.copy()
    u = np.zeros_like(fixed)
    anchor_slip,anchor_status = [],[]
    local_work_min = 0.
    width = 5e-5
    holding = g.initial-g.load
    path = [dict(name="geostatic",target_factor=0.),*case["stages"]]
    origin = 0.
    for phase,stage in enumerate(path):
        stage_tries = 0
        for step in range(1,steps+1):
            stage_tries += 1
            assert stage_tries <= (2065*len(mesh["tetrahedra"])+11)//12
            t = step/steps
            factor = origin + (stage["target_factor"]-origin)*t
            load = g.load + ((1-t)*holding if phase==0 else 0.)
            rhs = load-g.initial-factor*gfixed
            inner = cho_solve(inside,rhs[i])
            reduced_rhs = rhs[a]-kai@inner
            scale = max(1.,np.linalg.norm(load))
            counts = np.ones(n,dtype=int)
            x = u[a].copy()
            def evaluate(x):
                jump = (c@x).reshape(n,3)
                q,flow,status = project(g.system,old_jump,traction,jump[:,None,:],counts)
                res = condensed@x + c.T@(q[:,0]*g.area[:,None]).ravel()-reduced_rhs
                return jump,q[:,0],flow[:,0],status[:,0],res
            for iteration in range(60):
                jump,q,flow,status,res = evaluate(x)
                if np.linalg.norm(res)<1e-10*scale:
                    needed = np.maximum(1,np.ceil(np.linalg.norm(jump-old_jump,axis=1)/width)).astype(int)
                    if np.any(needed>counts):
                        counts = np.maximum(counts,needed)
                        continue
                    break
                h = 2e-9
                offsets = np.concatenate((np.eye(3),-np.eye(3)))*h
                probes = jump[:,None,:]+offsets
                pq,_,_ = project(g.system,old_jump,traction,probes,counts)
                local = ((pq[:,:3]-pq[:,3:])/(2*h)).transpose(0,2,1)
                jac = condensed + c.T@(g.area[:,None,None]*(local@ca)).reshape(3*n,-1)
                dx = np.linalg.solve(jac,-res)
                linear = np.linalg.norm(jac@dx+res)/max(np.linalg.norm(res),1e-30)
                assert linear <= 1e-10
                worst_linear = max(worst_linear,linear)
                for power in range(16):
                    new = x+dx/2**power
                    if np.linalg.norm(evaluate(new)[-1])<np.linalg.norm(res):
                        x = new
                        break
                else:
                    raise ArithmeticError(f"contact line search {stage['name']} {t} {iteration}")
            else:
                raise ArithmeticError(f"Newton iteration limit {stage['name']} {t}")
            u = fixed*factor
            u[a] = x
            u[i] = inner-lift@x
            # Reconstruct every interior displacement before checking the full
            # uncondensed force residual and the bulk elastic yield condition.
            full = g.k@u + g.initial-load
            full[a] += c.T@(q*g.area[:,None]).ravel()
            balance = np.linalg.norm(full[free])/scale
            assert balance < 1e-10
            worst_residual = max(worst_residual,balance)
            for b,dofs,d,s0,mc in g.operators:
                sigma = s0-(b@u[dofs])@d.T
                worst_yield = max(worst_yield,float(np.max(normalized_yield(mc,sigma))))
            assert worst_yield <= 1e-8
            local_work_min = min(local_work_min,float(np.min((flow*q).sum(axis=1))))
            material_steps += int(counts.sum())
            old_jump,traction = jump.copy(),q.copy()
            plastic += flow
            if step%(steps//64)==0:
                anchor_slip.append(plastic.copy())
                anchor_status.append(status.copy())
        origin = stage["target_factor"]
        attempts.append(stage_tries)
    assert len(anchor_slip)==448
    np.savez_compressed(folder/"slip.npz",slip=anchor_slip,status=anchor_status,weight=g.area,
                        position=g.positions,parent=g.parents)
    return dict(degree=degree,elements=len(mesh["tetrahedra"]),steps=steps,face_gauss_order=face_order,
                nodes=len(g.nodes),dofs=g.nodes.size,free_dofs=int(free.sum()),area_points=n,accepted=sum(attempts),
                attempts_per_stage=attempts,accepted_contact_partitions=material_steps,
                anchors=448,contact_area_m2=float(g.area.sum()),maximum_relative_force=float(worst_residual),
                maximum_relative_linear=float(worst_linear),maximum_normalized_bulk_yield=worst_yield,
                min_endpoint_flow_work=local_work_min,seconds=time.monotonic()-start,
                peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("output",type=Path)
    p.add_argument("--degree",type=int,choices=(2,3),required=True)
    p.add_argument("--steps",type=int,choices=(64,128),required=True)
    p.add_argument("--face-order",type=int,default=6)
    args=p.parse_args()
    args.output.mkdir()
    mesh,ground,case = [json.loads((ROOT/"examples"/(name+".json")).read_bytes()) for name in (
        "joint-water-column-tet10-mesh","ground-joint-water-mc","stages-joint-water-cycle")]
    result=reference(mesh,ground,case,args.degree,args.steps,args.face_order,args.output)
    (args.output/"cost.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
