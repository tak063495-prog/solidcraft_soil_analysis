"""Pair only accepted positive-direction observations on the same dry contact branch."""
import gzip
import hashlib
import json

import numpy as np

from sealed import HERE, archive

STUDY = HERE.parent / "2026-11-22-impl-t7-physical-instability"


def paired(a, b, target):
    for r in (a,b):
        if r.get("matched_side") is False:
            return "accepted points on different sides of the switch"
        if "fields" not in r:
            return "accepted field unavailable"
        if target:
            if r.get("matched_side") is not True or r.get("actual") != target:
                return "same side or progress unavailable"
            left, right = r["bracket"]
            if not left <= target <= right or right <= left:
                return "accepted direction differs"
    for key in ("stage", "generations", "event"):
        if a["fields"][key] != b["fields"][key]:
            return f"different {key}"
    return None


def support_reason(pair, case, mesh, progress):
    """Check every prescribed support direction on both accepted displacement fields."""
    stages = case["stages"]
    stage = pair[0]["fields"]["stage"]
    index = next(i for i, entry in enumerate(stages) if entry["name"] == stage)
    if index == 0 or stages[index]["kind"] != "boundary" or stages[index - 1]["kind"] != "boundary":
        return "proportional boundary history unavailable"
    origin = stages[index - 1]["target_factor"]
    factor = origin + progress * (stages[index]["target_factor"] - origin)
    for record in pair:
        u = np.asarray(record["fields"]["u"]).reshape(len(mesh["nodes"]), 3)
        for constraint in case["boundary"]["constraints"]:
            nodes = mesh["node_sets"][constraint["node_set"]]
            actual = u[nodes] @ constraint["direction"]
            if not len(nodes) or not np.all(actual == factor * constraint["value_m"]):
                return "prescribed support displacement differs"
    return None


def evaluate():
    spectra = json.loads((STUDY / "accepted-spectra.json").read_bytes())
    results = []
    with archive(STUDY, "archive-index.json") as z:
        comparison = json.loads(gzip.decompress(z.read("science/comparison.json.gz")))
        for n, data in comparison["grids"].items():
            raw = z.read(f"science/mesh-{n}/case.json")
            source = hashlib.sha256(raw).hexdigest()
            case = json.loads(raw)
            mesh = json.loads(z.read(f"science/mesh-{n}/mesh.json"))
            for engine in ("native", "independent"):
                trace = spectra[n][engine+"-0.02-1"]["accepted_points"]
                works = [p["work_kn_m"] for p in trace if p.get("work_kn_m") is not None]
                results.append(dict(elements=int(n), engine=engine, case_sha256=source,
                                    stage="birth" if (engine,n)==("independent","18") else "open",
                                    accepted_points=len(trace), accepted_increments=len(trace)-1,
                                    maximum_progress=max(r["progress"] for r in trace),
                                    weak_residual_projection_max=max(abs(r["weak_residual_projection"]) for r in trace),
                                    max_abs_interval_work_kn_m={k:max(abs(w[k]) for w in works) for k in works[0]}
                                    if works else {}))
            for progress in (0., .02, .125):
                key=str(progress)
                a,b=[data["engines"][e].get(key,{}) for e in ("native","independent")]
                reason=paired(a,b,progress) or support_reason((a,b),case,mesh,progress)
                row=dict(elements=int(n), stage="open", target=progress, scope="same support input",
                         case_sha256=source, hold=reason, fields=None)
                if reason is None:
                    fields={}
                    for field in ("u","stress","traction","slip","jump"):
                        x,y=[np.array(r["fields"][field]) for r in (a,b)]
                        if x.shape != y.shape: raise ValueError("different field layout")
                        delta=float(np.linalg.norm(x-y)); den=float(np.linalg.norm(y))
                        fields[field]=dict(absolute_l2=delta,reference_l2=den,relative_l2=delta/den if den else None)
                    row.update(fields=fields, generations=a["fields"]["generations"],
                               surface={k:{s:r["fields"]["surfaces"][s] for s in ("UpperSeam","LowerSeam")}
                                        for k,r in (("native",a),("independent",b))},
                               weak_space_angles_rad=data["common_progress"][key]["left_principal_angles_rad"])
                if reason is None and progress:
                    origins = [data["engines"][e]["0.0"] for e in ("native", "independent")]
                    if paired(*origins, 0.) is None and support_reason(origins, case, mesh, 0.) is None:
                        row["friction_work_from_zero_kn_m"] = {
                            seam: [
                                r["fields"]["surfaces"][seam]["friction_work_kn_m"] -
                                origin["fields"]["surfaces"][seam]["friction_work_kn_m"]
                                for r, origin in zip((a, b), origins, strict=True)
                            ] for seam in ("UpperSeam", "LowerSeam")
                        }
                results.append(row)
    return results
