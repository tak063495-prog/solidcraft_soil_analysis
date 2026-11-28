"""Derive the prospective width table and stratified parent list solely from P7."""
from collections import Counter, defaultdict
import gzip
import hashlib
import json
from pathlib import Path
from statistics import median

HERE = Path(__file__).resolve().parent
P7 = HERE.parent / "2026-11-26-remaining-p7"


def save(path, value):
    path.parent.mkdir(exist_ok=True, parents=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def causes(row):
    a = row["audit"]
    return [name for name, active in [
        ("event", a["event"] is not None),
        ("slip", a["max_slip"] > a["slip_budget"]),
        ("work", a["max_work"] > a["work_budget"]),
    ] if active]


def freeze_contact():
    bins = 256
    labels = [set() for _ in range(bins)]
    ratios = defaultdict(list)
    rejected = []
    hashes = {}
    for n in (256, 512, 1024):
        path = P7 / f"c01/original{n}/trace.jsonl.gz"
        hashes[str(n)] = digest(path)
        with gzip.open(path, "rt") as stream:
            rows = [r for r in map(json.loads, stream) if r.get("stage") == "reload"]
        for i, row in enumerate(rows):
            if row["accepted"]:
                continue
            why = causes(row)
            assert why
            if n != 1024:
                # Cover the actual coarse/middle rejected interval, including its parent.
                for k in range(int(row["from"]*bins), min(bins, int(row["to"]*bins)+1)):
                    labels[k].update(why)
                continue
            after = next((r for r in rows[i+1:] if r["accepted"]), None)
            ratio = None  # Last rejected series is right-censored at the P7 cap.
            if after is not None:
                assert after["from"] == row["from"]
                ratio = (after["to"]-after["from"])/(row["to"]-row["from"])
                for cause in why:
                    ratios[cause].append(ratio)
            rejected.append(dict(attempt=i+1, parent=row["from"], end=row["to"],
                                 causes=why, first_accepted_width_ratio=ratio))
    assert len(rejected) == 163
    # Use cause-conditioned ratios to estimate an absolute cap, avoiding
    # repeated fractional shrinkage of a proposal that is already small.
    factors = {k: median(v) for k,v in ratios.items()}
    caps = [min([1.] + [factors[k] for k in label])/1024 for label in labels]
    save(HERE/"c01/predictor.json", dict(
        source_sha256=hashes, bins=bins, causes_per_bin=[sorted(s) for s in labels],
        median_ratio=factors, caps=caps, training_rejections=rejected,
        rule="Take min(ordinary proposed width, cap at floor(256*parent)); preserve native target grid and all rejection conditions.",
        original_trials=1787, mesh_limit=2065, default_limit=2000))


def freeze_nc():
    path = P7/"nc/continuation.jsonl"
    rows = [json.loads(line) for line in path.read_text().splitlines()]
    groups = []
    for row in rows:
        if not groups or groups[-1][-1]["accepted"]:
            groups.append([])
        groups[-1].append(row)
    assert groups[-1][-1]["accepted"]
    # Keep a complete retry series in the stratum containing its first trial.
    # Sample first, middle and last parent per reported local failure point
    # in each third of the archived chronological attempt interval.
    intervals = []
    for lo, hi in ((106,282),(283,459),(460,637)):
        sample = [g for g in groups if lo <= g[0]["attempt"] <= hi]
        by_point = defaultdict(list)
        for group in sample:
            tags = sorted({r.get("reason","").split("active point ")[-1].split(",")[0]
                           for r in group if not r["accepted"]})
            for tag in tags:
                by_point[tag].append(group)
        selected = set()
        for collection in by_point.values():
            selected.update(collection[i][0]["attempt"] for i in (0,len(collection)//2,-1))
        intervals.append(dict(first=lo,last=hi,parents=len(sample),
            lambda_range=[sample[-1][0]["parent_lambda"],sample[0][0]["parent_lambda"]],
            trials=sum(map(len,sample)), rejected=sum(len(g)-1 for g in sample),
            rejected_parents=sum(len(g)>1 for g in sample),
            parents_by_failure_tag={k:len(v) for k,v in by_point.items()},
            selected=[dict(first=g[0]["attempt"],last=g[-1]["attempt"],
                           parent_lambda=g[0]["parent_lambda"],
                           widths=[r["ds_m"] for r in g]) for g in sample if g[0]["attempt"] in selected]))
    save(HERE/"nc/selection.json", dict(source_sha256=digest(path), points=[351,272],
        unit="series of retries from one committed parent, through first acceptance",
        selection="first, upper middle and last rejected parent per failure point and interval, deduplicated",
        intervals=intervals))


if __name__ == "__main__":
    freeze_contact()
    freeze_nc()
