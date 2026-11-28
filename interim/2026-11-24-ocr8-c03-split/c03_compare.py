"""Replay the accepted FEM strains continuously through the native MC return."""
import argparse
import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time

import numpy as np

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "2026-11-24-independent-local-references"
sys.path.insert(0, str(SOURCE))
from c03_point import LABELS, POINTS, WORK_LIMIT, REFERENCE_LIMIT, discrepancies, tensor, work

TARGETS = {108: (5, 31/32), 185: (2, 31/32)}


def read(path):
    return json.loads(gzip.decompress(path.read_bytes()) if path.suffix == ".gz" else path.read_bytes())


def replay(data, probe, mode):
    seeds = read(HERE.parents[1] / "examples/ground-joint-water-mc.json")["initialization"]["geostatic"]["seed"]["states"]
    seed = seeds["LowerClay"]
    assert all(s == seed for s in seeds.values())
    requests = []
    for route, history in data.items():
        assert history["points"] == list(POINTS)
        for i, frame in enumerate(history["frames"]):
            coord = (frame["stage"], frame["progress"])
            assert coord == (i//32, (i%32+1)/32)
            for j, gp in enumerate(history["points"]):
                if mode == "fine_parent" and (route != "m0c1" or gp != 108 or coord > TARGETS[108]):
                    continue
                count = frame["fem_substeps"][j]
                assert seed["effective_stress_compression_positive_kpa"] == history["initial_stress_kpa"][j]
                req = dict(route=route, gp=gp, increment=i, parameters=history["material"],
                           old=seed, strain=frame["strain"][j], count=count, save_substeps=True)
                if mode != "default" and route == "m0c1":
                    if ((gp == 185 and coord <= TARGETS[185]) or mode == "fine_parent"
                            or (gp == 108 and coord == TARGETS[108])):
                        req["count"] = max(128, count)
                    if gp == 108 and coord == TARGETS[108]:
                        req["locate_on"] = count
                requests.append(req)
    with tempfile.TemporaryDirectory() as tmp:
        inp, out = (Path(tmp)/n for n in ("request", "response"))
        inp.write_text(json.dumps(requests), encoding="utf8")
        start = time.perf_counter()
        subprocess.run([probe, inp, out, "--history"], check=True)
        elapsed = time.perf_counter() - start
        records = read(out)
    assert len(records) == len(requests)
    return records, elapsed


def assess(data, refs, checks, records, elapsed, mode):
    series, counts, parents, events = {}, {}, {}, []
    for rec in records:
        route, gp, i = rec["route"], rec["gp"], rec["increment"]
        key = route, gp
        j = list(POINTS).index(gp)
        row = data[route]["frames"][i]
        steps = rec["substeps"]
        assert len(steps) == rec["count"]
        stress = tensor([s["stress_kpa"] for s in steps])
        plastic = tensor([s["plastic_increment"] for s in steps], engineering=True)
        inc = work(stress, plastic).sum(axis=0)
        end = rec["state"]
        if mode == "default":
            assert end["last_return"] == row["fem_label"][j]
            np.testing.assert_allclose(inc, row["fem_work_increment_kpa"][j], atol=2e-15, rtol=0)
        labels, history = series.setdefault(key, ([], []))
        assert len(labels) == i
        labels.append(LABELS.index(end["last_return"]))
        history.append(inc)
        counts[key] = counts.get(key, 0) + len(steps)
        if route == "m0c1" and (row["stage"], row["progress"]) == TARGETS.get(gp):
            eig = np.linalg.eigvalsh(stress)
            events.append(dict(gp=gp, preceding_frames=i, parent=parents[key],
                               parent_substeps=counts[key]-len(steps),
                               accepted_substeps=rec["count"], searches=rec.get("search_steps",0),
                               boundaries=rec.get("boundaries",[]), labels=rec["returns"],
                               labels_in_order=[s["label"] for s in steps],
                               corrected_upper_gap_kpa=(eig[:,2]-eig[:,1]).tolist(),
                               end_label=end["last_return"], end_stress=end["effective_stress_compression_positive_kpa"]))
        parents[key] = end
    points = []
    for (route,gp), (labels, increments) in series.items():
        j = list(POINTS).index(gp)
        check = checks["routes"][route][j]
        ref = refs[route]
        n = len(labels)
        assert n == (191 if mode == "fine_parent" else 224)
        err = discrepancies(np.cumsum(increments, axis=0), np.asarray(ref["reference_running_work_kpa"])[:n,j])
        diff = np.flatnonzero(np.array(labels) != np.asarray(ref["reference_label"])[:n,j])
        points.append(dict(route=route,gp=gp,substeps=counts[route,gp],matched=n-len(diff),
                           mismatched_frames=(diff+1).tolist(), relative_max=err.tolist(),
                           classification_passed=bool(check["reference_labels_stable"] and len(diff)==0),
                           work_passed=bool(all(err<=WORK_LIMIT) and max(check["refinement"][-1])<=REFERENCE_LIMIT)))
    return dict(seconds=elapsed, accepted_substeps=sum(counts.values()),
                search_substeps=sum(e["searches"] for e in events),
                matched=sum(p["matched"] for p in points), work_passed=sum(p["work_passed"] for p in points),
                work_max=np.max([p["relative_max"] for p in points],axis=0).tolist(),
                points=points, events=events)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("probe",type=Path)
    p.add_argument("output",type=Path)
    a = p.parse_args()
    data, refs, checks = [read(SOURCE / name) for name in ("c03-inputs.json.gz","c03-history.json.gz","c03.json")]
    assert set(data) == {"m0c1","m1c0","m1c1"}
    result = {}
    for mode in ("default", "refined", "fine_parent"):
        records, elapsed = replay(data, a.probe.resolve(), mode)
        result[mode] = assess(data, refs, checks, records, elapsed, mode)
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf8")
