"""Local static continuation with the original single holding-load coefficient."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import sys
import time
import zipfile
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/"tools"))
from verify_high_order import prepare, physical, internal


def digest(data):
    return hashlib.sha256(data).hexdigest()


class Section:
    def __init__(self, inputs, binary):
        self.binary = binary
        source = HERE.parent/"2026-11-23-nc01-coordinate-polish"
        metadata = json.loads((source/"measurements.json").read_bytes())
        for name, sha in metadata["input_sha256"].items():
            assert digest((inputs/name).read_bytes()) == sha
        with zipfile.ZipFile(source/"evidence.zip") as archive:
            raw = archive.read("run.json")
        self.run_sha = digest(raw)
        assert self.run_sha == metadata["run_sha256"]
        run = json.loads(raw)
        self.cp = cp = run["last_committed"]
        self.start = cp["stage_progress"]
        assert self.start == 0.5647809505462646
        mesh = json.loads((inputs/"mesh.json").read_bytes())
        self.settings = settings = run["case"]["settings"]
        live = np.array(cp["active_elements"])
        assert len(live) == 96 and live.sum() == 48
        assert not any(live[p["owner_a"]] and live[p["owner_b"]] for p in run["joint_integration_points"])
        self.hold = np.array(run["excavation_audits"][0]["initial_holding_load_kn"])
        q, offsets = prepare(mesh,run["integration_points"])
        ground = run["ground"]
        self.base = physical(mesh,ground,live,cp["surface_pressure_kpa"],run["case"]["pressure_on_original_surface_only"])[0].ravel()
        assert ground["initialization"]["water_level_m"] is None
        np.testing.assert_allclose(self.load(self.start),cp["applied_nodal_load_kn"],atol=3e-12,rtol=0)
        p1 = HERE.parent/"2026-11-25-remaining5-p1/nc/increment-4.npz"
        with np.load(p1) as f:
            self.free = f["free"]
            last, = [a for a in run["attempts"] if a["accepted"] and a["stage"]=="excavate_upper" and a["progress_to"]==self.start]
            secant=np.r_[f["du"], self.start-last["progress_from"]]
            self.seed=secant/np.linalg.norm(secant)
        nodes=np.unique(np.array(mesh["tetrahedra10"])[live])
        expected=(3*np.setdiff1d(nodes,mesh["node_sets"]["bottom"])[:,None]+np.arange(3)).ravel()
        np.testing.assert_array_equal(np.sort(self.free), expected)
        self.h = self.hold[self.free]
        materials=ground["materials"]
        ids=[m["id"] for m in materials]
        self.laws=[dict(m["mcc"],integration=settings["mcc_integration"]) for m in materials]
        self.b,self.w,self.dofs,self.ids,self.point_ids=[],[],[],[],[]
        for e in np.flatnonzero(live):
            self.b.append(q[e]["b"])
            self.w.append(q[e]["weight"])
            self.dofs.append((3*q[e]["nodes"][:,None]+np.arange(3)).ravel())
            points=list(range(offsets[e],offsets[e+1]))
            self.ids.extend([ids.index(cp["element_material_ids"][e])]*len(points))
            self.point_ids.extend(points)
        self.b,self.w,self.dofs=map(np.array,(self.b,self.w,self.dofs))
        self.initial=[cp["material_states"][j] for j in self.point_ids]
        force,_=internal(mesh,ground,cp["material_states"],live,q,offsets,material_ids=cp["element_material_ids"])
        self.seed_force=np.linalg.norm((force.ravel()-self.load(self.start))[self.free])
        assert self.seed_force < 2e-10
        self.width=settings["material_maximum_substep"]

    def load(self, lam):
        return self.base+(1-lam)*self.hold

    def evaluate(self, old, du, counts):
        flat=np.zeros_like(self.base)
        flat[self.free]=du
        de=-np.einsum("eqai,ei->eqa",self.b,flat[self.dofs]).reshape(-1,6)
        needed=np.maximum(1,np.ceil(np.abs(de).sum(axis=1)/self.width)).astype(int)
        if np.any(needed>counts):
            if np.max(needed)>50_000:
                raise ArithmeticError("material partition limit")
            counts[:] = np.where(needed>counts,2*needed,counts)
            return None
        request=[self.laws,self.ids,old,de.tolist(),counts.tolist()]
        call=subprocess.run([self.binary],input=json.dumps(request).encode(),capture_output=True)
        if call.returncode:
            raise ArithmeticError(call.stderr.decode())
        states,c,active=zip(*json.loads(call.stdout))
        c=np.array(c).reshape(48,27,6,6).swapaxes(-1,-2)
        stress=np.array([s["effective_stress_compression_positive_kpa"] for s in states]).reshape(48,27,6)
        f=-np.einsum("eqai,eqa,eq->ei",self.b,stress,self.w,optimize=True)
        ke=np.einsum("eqai,eqab,eqbj,eq->eij",self.b,c,self.b,self.w,optimize=True)
        force=np.zeros_like(flat)
        k=np.zeros((len(flat),len(flat)))
        for ids,fi,ki in zip(self.dofs,f,ke):
            force[ids]+=fi
            k[np.ix_(ids,ids)]+=ki
        return force[self.free],k[np.ix_(self.free,self.free)],states,int(np.count_nonzero(active))


def run(inputs,binary,out,ds,steps):
    start=time.monotonic()
    out.mkdir(parents=True)
    m=Section(inputs,binary)
    old,lam,t=m.initial,m.start,m.seed.copy()
    parent_u=np.array(m.cp["displacement_m"]).ravel()
    force_tol,lin_tol=(m.settings[k] for k in ("relative_force_tolerance","relative_linear_tolerance"))
    rows,attempts=[],[]
    report=dict(run_sha256=m.run_sha,start=lam,seed_free_force_kn=m.seed_force,
                step_m=ds,requested_steps=steps,original_force_gate=force_tol,original_linear_gate=lin_tol,
                fixed_budget=2000,mesh_budget=16520,material_point_indices=m.point_ids)
    stop=None
    with gzip.open(out/"accepted.jsonl.gz","wt",compresslevel=3) as history:
        for step in range(steps):
            assert step < 2000
            counts=np.ones(len(old),dtype=int)
            dz=t*ds
            record=dict(step=step+1,parent_lambda=lam,iterations=[])
            attempts.append(record)
            try:
                for it in range(m.settings["maximum_newton_iterations"]+1):
                    load=lam+dz[-1]
                    evaluated=m.evaluate(old,dz[:-1],counts)
                    if evaluated is None:
                        record["iterations"].append(dict(iteration=it,repartition=int(counts.max())))
                        dz=t*ds
                        continue
                    force,k,new,plastic=evaluated
                    external=m.load(load)
                    rhs=external[m.free]-force
                    rel=np.linalg.norm(rhs)/max(1.,np.linalg.norm(external))
                    g=ds-t@dz
                    item=dict(iteration=it,trial_lambda=load,relative_force=rel,arc_residual_m=g)
                    record["iterations"].append(item)
                    # One corrector is mandatory to resolve the shallow load extremum.
                    if any("augmented_true_linear" in r for r in record["iterations"]) and rel<=min(1e-12,force_tol) and abs(g)<=1e-12:
                        break
                    a=np.vstack((np.column_stack((k,m.h)),t))
                    b=np.r_[rhs,g]
                    dx=np.linalg.solve(a,b)
                    # Audit the physical force equations on their own unscaled RHS.
                    lr=np.linalg.norm(a@dx-b)/max(np.linalg.norm(b),1e-30)
                    fr=np.linalg.norm(k@dx[:-1]+m.h*dx[-1]-rhs)/max(np.linalg.norm(rhs),1e-30)
                    item.update(augmented_true_linear=lr,force_true_linear=fr)
                    if not max(lr,fr)<=lin_tol:
                        raise ArithmeticError("true linear residual")
                    dz+=dx
                else:
                    raise ArithmeticError("Newton iteration limit")
                left,s,vh=np.linalg.svd(k)
                sign=np.sign(left[:,-1]@m.h) or 1
                oriented=sign*left[:,-1]
                mode=vh[-1]*(np.sign(vh[-1]@dz[:-1]) or 1)
                row=dict(step=step+1,arc_m=(step+1)*ds,load=load,delta_lambda=dz[-1],
                         relative_force=rel,arc_residual_m=g,plastic_points=plastic,
                         partitions_max=int(counts.max()),sigma_min=s[-1],sigma_next=s[-2],load_projection=oriented@m.h,
                         relative_projection=(oriented@m.h)/np.linalg.norm(m.h),
                         signed_soft_stiffness=oriented@k@mode,
                         increment_work=dz[:-1]@k@dz[:-1])
                parent_u[m.free]+=dz[:-1]
                history.write(json.dumps(dict(step=step+1,load=load,displacement_m=parent_u.reshape(-1,3).tolist(),
                                              active_material_states=new),allow_nan=False)+"\n")
                rows.append(row)
                old,lam,t=new,load,dz/np.linalg.norm(dz)
                record["accepted"]=True
            except (ArithmeticError,np.linalg.LinAlgError) as e:
                stop=str(e)
                record.update(accepted=False,reason=stop)
                break
    report.update(rows=rows,attempts=attempts,stop=stop,seconds=time.monotonic()-start,
                  python_peak_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  material_worker_peak_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    (out/"report.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs",type=Path)
    parser.add_argument("binary",type=Path)
    parser.add_argument("out",type=Path)
    parser.add_argument("step",type=float,choices=[0.0008,0.0004,0.0002])
    args=parser.parse_args()
    run(args.inputs,args.binary.resolve(),args.out,args.step,round(0.0128/args.step))
