"""Apply the p11 symmetric endpoint decomposition to p11's accepted half steps."""
import gzip
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "2026-11-26-remaining-p11/c01/alternate"
LEVELS = (256, 512, 1024)


def halves(path):
    candidates = {}
    attempts = []
    with gzip.open(path, "rt") as stream:
        for line in stream:
            if not line.startswith("{"):
                continue
            row = json.loads(line)
            if "p11_span" in row:
                l, m, r = row["p11_span"]
                candidates[l, r] = (m, row["all"])
            if row.get("stage") != "recontact":
                continue
            attempts.append(row)
    for attempt, row in enumerate(attempts, 1):
        if row["accepted"]:
            l, r = row["from"], row["to"]
            m, (x, y, z) = candidates[l, r]
            yield attempt, 1, l, m, x, y
            yield attempt, 2, m, r, y, z


def run(out):
    base_knots = json.loads((BASE / "plan.json").read_text())["knots"]
    knots = set(base_knots)
    paths = [list(halves(BASE / str(n) / "trace.jsonl.gz")) for n in LEVELS]
    clocks = []
    for path in paths:
        clock = {0.: path[0][-2]}
        for _, _, l, r, x, y in path:
            assert clock[l] == x
            clock[r] = y
        clocks.append(clock)
    shared = sorted(set.intersection(*(set(c) for c in clocks)))
    assert shared[0] == 0 and shared[-1] == 1
    points = {}
    reported = json.loads((BASE / "contributions.json").read_text())["recontact"]["p11"]["endpoint"]
    with gzip.open(out / "signed-halves.jsonl.gz", "xt") as stream:
        for point in (14, 26):
            levels = np.array([c[1][point][1] for c in clocks])
            A, B = -np.diff(levels, axis=0)
            np.testing.assert_allclose(A @ B, reported[point], rtol=1e-14, atol=0)
            weights = [B/2, (A-B)/2, -A/2]
            bins = []
            for t in shared:
                a, b = -np.diff([c[t][point][1] for c in clocks], axis=0)
                bins.append(float((a @ B + b @ A)/2))
            j = int(np.argmin(np.diff(bins)))
            lo, hi = shared[j:j+2]
            local, changes, sums = [], [], []
            for n, path, weight in zip(LEVELS, paths, weights, strict=True):
                total = 0.
                for attempt, half, l, r, x, y in path:
                    delta = np.subtract(y[point][1], x[point][1])
                    term = float(weight @ delta)
                    total += term
                    row = dict(point=point, level=n, attempt=attempt, half=half, span=[l,r],
                               labels=[x[point][0],y[point][0]], slip_increment_m=delta.tolist(), term_m2=term)
                    print(json.dumps(row), file=stream)
                    if lo <= l < r <= hi:
                        local.append(row)
                    if row["labels"] == ["open","slip"]:
                        changes.append(row)
                        if n == 512:
                            knots.update((l,r))
                sums.append(total)
            assert len(changes) == 3 and {r["level"] for r in changes} == set(LEVELS)
            knots.update((min(s["span"][0] for s in changes), max(s["span"][1] for s in changes)))
            np.testing.assert_allclose(bins[0]+sum(sums), A @ B, rtol=1e-8, atol=1e-28)
            points[point] = dict(endpoint_m2=float(A@B), initial_m2=bins[0], half_sums_m2=sums,
                                 worst_common_interval=[lo,hi], interval_term_m2=bins[j+1]-bins[j],
                                 accepted_halves=local, open_to_slip=changes,
                                 shared_interval_terms=[[l,r,d] for l,r,d in zip(shared[:-1],shared[1:],np.diff(bins))])
    plan = dict(p11_knots=base_knots, knots=sorted(knots), added=sorted(knots-set(base_knots)), points=points)
    (out / "plan.json").write_text(json.dumps(plan, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    run(Path(sys.argv[1]))
