"""Independent differential (explicit) MCC integration vs Rust implicit return.

The reference advances a 6-component elastoplastic rate equation, not the Rust
four-unknown Newton equations. It intentionally uses much smaller strain steps.
"""
import argparse
import json
from pathlib import Path
import subprocess
import numpy as np


def invariants(stress):
    p = stress[:3].mean()
    dev = stress.copy()
    dev[:3] -= p
    q = np.sqrt(1.5*(np.square(dev[:3]).sum()+2*np.square(dev[3:]).sum()))
    return p,q,dev


def explicit_reference(params, initial, strain, count=40000):
    s = np.asarray(initial["effective_stress_compression_positive_kpa"],dtype=float)
    pc = initial["preconsolidation_pressure_kpa"]
    v = 1+initial["void_ratio"]
    de = np.asarray(strain)/count
    m2 = params["critical_state_m"]**2
    nu = params["poisson"]
    for _ in range(count):
        p,q,dev = invariants(s)
        bulk = v*p/params["kappa"]
        shear = 3*bulk*(1-2*nu)/(2*(1+nu))
        d = np.zeros((6,6))
        d[:3,:3] = bulk-2*shear/3
        d[np.arange(3),np.arange(3)] += 2*shear
        d[np.arange(3,6),np.arange(3,6)] = shear
        trial = s+d@de
        pt,qt,_ = invariants(trial)
        if qt**2+m2*pt*(pt-pc) > 0:
            gradient = 3*dev
            gradient[:3] += m2*(2*p-pc)/3
            gradient[3:] *= 2
            stiffness_gradient = d@gradient
            loading = gradient@(d@de)
            hardening_modulus = m2*p*pc*v/(params["lambda"]-params["kappa"])*gradient[:3].sum()
            multiplier = max(loading/(gradient@stiffness_gradient+hardening_modulus),0.)
            s += d@(de-multiplier*gradient)
            pc += pc*v/(params["lambda"]-params["kappa"])*multiplier*gradient[:3].sum()
        else:
            s = trial
        v *= np.exp(-de[:3].sum())
    return {"stress":s,"pc":pc,"e":v-1}


def verify(solver, output):
    output.mkdir(parents=True,exist_ok=True)
    params = {"critical_state_m":1.2,"lambda":0.2,"kappa":0.04,"poisson":0.3,
              "reference_pressure_kpa":1.,"specific_volume_ncl_at_reference":3.5}
    reports = []
    for name,pc,strain in [
        ("nc_isotropic",100.,[.01,.01,.01,0,0,0]),
        ("oc_isotropic",400.,[.01,.01,.01,0,0,0]),
        ("nc_undrained_shear",100.,[-.015,-.015,.03,0,0,0]),
        ("oc_undrained_shear",400.,[-.015,-.015,.03,0,0,0]),
    ]:
        initial = {"effective_stress_compression_positive_kpa":[100.,100.,100.,0.,0.,0.],
                   "void_ratio":3.5-.2*np.log(pc)+.04*np.log(pc/100)-1,
                   "preconsolidation_pressure_kpa":pc,"plastic_strain_compression_positive":[0.]*6}
        reference = explicit_reference(params,initial,strain)
        errors = []
        for suffix,h in (("coarse",1e-4),("fine",1e-5)):
            path = output/f"{name}-{suffix}.json"
            result_path = output/f"{name}-{suffix}-result.json"
            path.write_text(json.dumps({"schema":"soil-fem-mcc-path/1","parameters":params,"initial_state":initial,
                                       "strain_increments_compression_positive":[strain],"maximum_substep":h},indent=2),encoding="utf-8")
            subprocess.run([str(solver),"mcc-path",str(path),str(result_path)],check=True,capture_output=True,text=True)
            result = json.loads(result_path.read_text())["states"][-1]
            stress = np.asarray(result["effective_stress_compression_positive_kpa"])
            error = np.linalg.norm(stress-reference["stress"])/np.linalg.norm(reference["stress"])
            errors.append(error)
            if suffix == "fine":
                assert error<1e-3,(name,error)
                assert abs(result["preconsolidation_pressure_kpa"]/reference["pc"]-1)<1e-3
                p,q,_ = invariants(stress)
                reports.append({"case":name,"p_kpa":p,"q_kpa":q,
                                "pc_kpa":result["preconsolidation_pressure_kpa"],
                                "relative_stress_error_vs_explicit_rate_reference":error})
        assert errors[1]<errors[0],(name,errors)
        reports[-1]["coarse_error"] = errors[0]
    nc = next(r for r in reports if r["case"]=="nc_undrained_shear")
    oc = next(r for r in reports if r["case"]=="oc_undrained_shear")
    assert nc["p_kpa"]<100 and oc["p_kpa"]>100
    (output/"mcc-verification.json").write_text(json.dumps(reports,indent=2),encoding="utf-8")
    print(json.dumps(reports,indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--solver",type=Path,required=True)
    parser.add_argument("--output-dir",type=Path,required=True)
    a = parser.parse_args()
    verify(a.solver.resolve(),a.output_dir.resolve())
