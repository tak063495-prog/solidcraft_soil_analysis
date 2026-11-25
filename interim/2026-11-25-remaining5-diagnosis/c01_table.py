import gzip
import json
import sys
from pathlib import Path
from collections import Counter
import numpy as np

root = Path.cwd()
pack = root / "interim/2026-11-25-remaining5-diagnosis"
old = json.load(
    gzip.open(root / "interim/2026-11-25-approved-six/c01-reload-audit.json.gz")
)
new = json.loads(Path(sys.argv[1]).read_bytes())
assert (new["last"], new["attempts"], new["complete"]) == (
    old["last"],
    old["attempts"],
    False,
)
assert len(new["records"]) == len(old["records"]) == 234
rank = Counter()
labels = Counter()
ranges = {}
for row, prior in zip(new["records"], old["records"], strict=True):
    for k, x in prior.items():
        assert row[k] == x, (k, row[k], x)
    if row["accepted"] or row["event"] is not None:
        continue
    a = np.asarray(row["slip_errors_m"])
    i, j = np.unravel_index(np.abs(a).argmax(), a.shape)
    rank[(int(i), int(j))] += 1
    labels[(row["old"][i], row["mid"][i], row["coarse"][i], row["fine"][i])] += 1
    for name in ("mid_force", "fine_force"):
        ranges.setdefault(name, []).append(row[name])
report = dict(
    identical_archived_trials=len(new["records"]),
    peak_point_component_counts=[
        dict(point=i, component=j, count=v) for (i, j), v in rank.most_common()
    ],
    peak_status_sequences=[dict(statuses=k, count=v) for k, v in labels.most_common()],
    force_min_median_max={
        k: np.quantile(v, [0, 0.5, 1]).tolist() for k, v in ranges.items()
    },
)
print(json.dumps(report, indent=2))
