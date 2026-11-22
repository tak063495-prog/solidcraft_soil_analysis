"""Fixed physical point and full plate area comparisons for B-02, 12/96 Tet10."""
import argparse
import json
from pathlib import Path
import shutil
import sys

import numpy as np

from sealed import HERE, archive

STUDY = HERE.parent / "2026-12-01-tri6-offset-composite-b02"
sys.path.insert(0, str(STUDY))
import convergence as old

# Six-point degree-four Dunavant rule, normalized to the physical triangle area.
WEIGHTS = (0.109951743655322,) * 3 + (0.223381589678011,) * 3
PATHS = ("beam-positive", "mesh-1")
STREAMS = {"native": "native/accepted.jsonl", "independent": "independent.jsonl"}
FIELDS = ("plate_membrane_global_kn_m", "plate_bending_global_kn",
          "plate_attachment_moment_global_kn", "plate_shear_global_kn_m")


def section_integrals(folder, engine):
    cfg = old.read(folder / "stages.json")["beam_system"]
    offsets = {p["id"]: p["attachment_offset_m"] for p in cfg["plates"]}
    result = {}
    with (folder / STREAMS[engine]).open() as stream:
        for line in stream:
            row = json.loads(line)
            key = row["stage_index"], row["progress"]
            if key[1] not in (.25, .5, .75, 1.):
                continue
            plates = row["beam_results"].get("plates", [])
            items = plates.items() if isinstance(plates, dict) else ((p["id"], p) for p in plates)
            total = {}
            for name, plate in items:
                axes = np.array(plate["local_axes"])
                for gp, weight in zip(plate["samples"], WEIGHTS, strict=True):
                    n = old.tensor(gp["membrane_force_kn_m"], axes)
                    m = old.tensor(gp["bending_moment_kn"], axes)
                    q = axes.T @ np.r_[gp["transverse_shear_force_kn_m"], 0.]
                    for field, value in zip(FIELDS, (n, m, m-offsets[name]*n, q), strict=True):
                        total[field] = total.get(field, 0.) + plate["area_m2"] * weight * value
            result[key] = total
    return result


def evaluate(root):
    with archive(STUDY, "evidence-manifest.json") as z:
        for case in PATHS:
            for name in ("mesh.json", "ground.json", "stages.json", "refinement.json", *STREAMS.values()):
                dest = root / case / name
                dest.parent.mkdir(parents=True, exist_ok=True)
                with z.open(f"runs/b02-formal-01/{case}/{name}") as src, dest.open("xb") as output:
                    shutil.copyfileobj(src, output)
    mesh = old.read(root / PATHS[0] / "mesh.json")
    cfg = old.read(root / PATHS[0] / "stages.json")["beam_system"]
    model = old.System(mesh, cfg)
    origin = model.positions[model.node_ids[cfg["members"][0]["nodes"][0]]]
    xyz = [np.array(l) @ model.positions[nodes][:3]
           for nodes, _, _ in model.plate_data for l in old.SAMPLES]
    result = dict(probes_m=np.asarray(xyz).tolist(), connection_m=np.asarray(origin).tolist(),
                  probe_distance_m=[float(np.linalg.norm(x-origin)) for x in xyz], rows=[])
    for engine in STREAMS:
        series = [old.sample(root / case, mesh, engine)[0] for case in PATHS]
        areas = [section_integrals(root / case, engine) for case in PATHS]
        if not all(x.keys() == series[0].keys() for x in (*series, *areas)):
            raise ValueError("accepted coordinates differ")
        for point in series[0]:
            a, b = [s[point] for s in series]
            if a.keys() != b.keys():
                raise ValueError("active field inventory differs")
            row = dict(engine=engine, stage=point[0], progress=point[1], fields={})
            for key in (*FIELDS, "load_work_kn_m", "support_resultant_kn"):
                if key not in a:
                    continue
                data = dict(coarse=a[key].tolist(), fine=b[key].tolist(),
                            absolute_l2=float(np.linalg.norm(a[key]-b[key])))
                if key in FIELDS:
                    ca, fa = [area[point][key] for area in areas]
                    data.update(area_integrals=[ca.tolist(), fa.tolist()],
                                integral_absolute_l2=float(np.linalg.norm(ca-fa)))
                row["fields"][key] = data
            result["rows"].append(row)
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    old.write(args.output / "comparison.json", evaluate(args.output))
