"""Compare the same accepted coordinates and the same signed plastic-work components."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
from zipfile import ZipFile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "2026-11-23-impl-root-causes-c"))
from c03 import STAGES, measure, read
from c03diag_study import write

ROUTES = ("m0c1", "m1c0", "m1c1")
LABELS = ("material_return", "contact_status")


def sha(p):
    with p.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def stream(p):
    with p.open(encoding="utf-8") as f:
        for line in f:
            yield json.loads(line)


def first_differences(folders, differences):
    # Report the first nonidentical stored stress/jump separately from the first
    # categorical divergence. A bit-level difference carries no acceptance width.
    points = {(r["field"], r["gp"]) for r in differences}
    result = {f"{field}/{gp}": {} for field, gp in points}
    if not points:
        return result
    histories = [stream(path / "fields.jsonl") for path in folders]
    for a, b in zip(*histories, strict=True):
        coord = {k: a["row"][k] for k in ("stage", "progress")}
        assert all(b["row"][k] == v for k, v in coord.items())
        for field, gp in sorted(points):
            record = result[f"{field}/{gp}"]
            quantity = "stress" if field == "material_return" else "jump"
            left, right = (sample["fields"][quantity][gp] for sample in (a, b))
            if left != right and "first_field_difference" not in record:
                record["first_field_difference"] = dict(
                    **coord,
                    quantity=quantity,
                    max_component_difference=max(abs(x - y) for x, y in zip(left, right, strict=True)),
                )
            labels = [sample["fields"][field][gp] for sample in (a, b)]
            if labels[0] != labels[1] and "first_classification_difference" not in record:
                record["first_classification_difference"] = coord
    return result


def compare(matrix, after, work, mcc, regressions, out, baseline):
    historical = baseline / "interim/2026-11-23-impl-root-causes-c"
    native_map = matrix / "m1c1-native/raw/run.json"
    raw = read(native_map)
    counts = dict(
        material_return=len(raw["integration_points"]),
        contact_status=len(raw["joint_integration_points"]),
    )
    anchors = [(s, i / 32) for s in STAGES for i in range(1, 33)]
    assert counts == dict(material_return=324, contact_status=32)
    comparisons = {}
    for file in sorted((historical / "results").glob("c03-*.json")):
        saved = read(file)
        entry = {}
        for mode in ("before", "after"):
            folders = [
                (matrix if (mode == "before" or name.endswith("-native")) else after) / name
                for name in (saved["left"], saved["right"])
            ]
            result = measure(matrix, *folders, native_map, saved["axis"])
            assert result["accepted"] == len(anchors) and len(result["full_fields"]) == 23
            if mode == "before":
                assert (
                    result["full_fields"] == saved["full_fields"]
                    and result["mismatches"] == saved["mismatches"]
                )
            entry[mode] = {
                f: dict(
                    matches=len(anchors) * counts[f] - result["full_fields"][f]["mismatches"],
                    total=len(anchors) * counts[f],
                )
                for f in LABELS
            }
            entry[mode]["differences"] = [
                {k: r[k] for k in ("stage", "progress", "field", "gp", "left", "right")}
                for r in result["mismatches"]
            ]
            entry[mode]["onset"] = first_differences(folders, result["mismatches"])
            if mode == "after":
                entry["after_full_fields"] = result["full_fields"]
        comparisons[file.stem] = entry
    assert len(comparisons) == 7
    work_summaries = {}
    digests = {}
    csv_path = out / "mc-work.csv"
    with csv_path.open("x", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(
            (
                "route",
                "solver",
                "stage",
                "progress",
                "before_substeps",
                "after_substeps",
                "before_direct",
                "after_direct",
                "before_volume",
                "after_volume",
                "before_shear",
                "after_shear",
            )
        )
        for route in ROUTES:
            for engine in ("native", "independent"):
                paths = (
                    [work / (route + "-native.jsonl")] * 2
                    if engine == "native"
                    else [work / (route + "-" + mode + ".jsonl") for mode in ("before", "after")]
                )
                rows = [list(stream(p)) for p in paths]
                for side in rows:
                    assert [(r["stage"], r["progress"]) for r in side] == anchors
                sums = [[math.fsum(r["work_kn_m"][i] for r in side) for i in range(3)] for side in rows]
                name = route + "-" + engine
                b, a = sums
                work_summaries[name] = dict(
                    before=b,
                    after=a,
                    delta=[y - x for x, y in zip(b, a)],
                    max_increment_direct_delta=max(
                        abs(r["work_kn_m"][0] - s["work_kn_m"][0]) for r, s in zip(*rows, strict=True)
                    ),
                    before_split_residual=b[0] - (b[1] + b[2]),
                    after_split_residual=a[0] - (a[1] + a[2]),
                )
                for r, s in zip(*rows, strict=True):
                    writer.writerow(
                        (
                            route,
                            engine,
                            r["stage"],
                            r["progress"],
                            r["substeps"],
                            s["substeps"],
                            *(x for pair in zip(r["work_kn_m"], s["work_kn_m"]) for x in pair),
                        )
                    )
                before = matrix / name
                revised = after / name if engine == "independent" else before
                assert (before / "stages.json").read_bytes() == (revised / "stages.json").read_bytes()
                digests[name] = dict(
                    before_fields=sha(before / "fields.jsonl"),
                    after_fields=sha(revised / "fields.jsonl"),
                    case=sha(before / "stages.json"),
                )
    fixed = baseline / "interim/2026-11-23-c03-mcc-plastic-work/results"
    mc_old = read(fixed / "summary.json")
    mc_new = read(mcc / "results/summary.json")
    assert mc_old["totals"] == mc_new["totals"] and mc_old["pairs"] == mc_new["pairs"]
    assert sha(fixed / "increments.csv") == sha(mcc / "results/increments.csv")
    previous = baseline / "interim/2026-11-23-impl-root-causes-a-nc"
    pinned = baseline / "interim/2026-11-23-ac-scope-decisions/results"
    checks = {}
    pairs = [
        ("A-01", pinned / "a01.json", regressions / "a01.json"),
        ("A-02 shear", pinned / "a02.json", regressions / "shear/comparison.json"),
        ("C-02", historical / "c02-results.json", regressions / "c02.json"),
    ]
    for name, original, repeat in pairs:
        assert sha(original) == sha(repeat), name
        checks[name] = sha(repeat)
    with ZipFile(previous / "observations.zip") as z:
        for name, key, folder in (
            ("A-02 compression", "a02", "compression"),
            ("A-04", "a04", "saturated-longdouble"),
        ):
            comparison = regressions / folder / "results.json"
            old = json.loads(z.read(key + "/results.json"))
            new = read(comparison)
            if name == "A-04":
                assert new.pop("reference") == "longdouble"
            assert old == new, name
            checks[name] = dict(identical_keys=list(new), observed=sha(comparison))
    return dict(
        classifications=comparisons,
        mc_work_kn_m=work_summaries,
        mc_work_csv_sha256=sha(csv_path),
        mcc_work_kn_m={
            k: dict(before=v["work_kn_m"], after=mc_new["totals"][k]["work_kn_m"])
            for k, v in mc_old["totals"].items()
        },
        mcc_increments_sha256=sha(fixed / "increments.csv"),
        regressions=checks,
        accepted_fields_sha256=digests,
        matrix_inputs_sha256={
            name: sha(matrix / name) for name in ("mesh.json", "ground.json", "stages.json")
        },
        native_trace_sha256={
            route: {
                name: sha(matrix / (route + "-trace") / name)
                for name in ("accepted.jsonl", "constitutive.jsonl")
            }
            for route in ROUTES
        },
        quantitative_status="参照なし（全文待ち）・未達",
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("matrix", "after", "work", "mcc", "regressions", "output"):
        p.add_argument(name, type=Path)
    p.add_argument("--baseline-checkout", type=Path, default=HERE.parents[1])
    a = p.parse_args()
    result = compare(a.matrix, a.after, a.work, a.mcc, a.regressions, a.output, a.baseline_checkout)
    write(a.output / "comparison.json", result)
