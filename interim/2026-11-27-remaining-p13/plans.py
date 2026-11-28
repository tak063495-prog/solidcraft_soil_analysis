"""Separate only the p12 knots attributed to each measured switch."""
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
P12 = HERE.with_name("2026-11-27-remaining-p12")
P11 = HERE.with_name("2026-11-26-remaining-p11")


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def freeze(out):
    plan = json.loads((P12 / "c01/alternate/plan.json").read_text())
    base = set(plan["p11_knots"])
    additions = []
    for point in (14, 26):
        changes = plan["points"][str(point)]["open_to_slip"]
        middle, = [r["span"] for r in changes if r["level"] == 512]
        selected = set(middle) | {min(r["span"][0] for r in changes),
                                  max(r["span"][1] for r in changes)}
        extra = selected - base
        additions.append(extra)
        dest = out / f"point{point}"
        dest.mkdir(parents=True, exist_ok=True)
        save(dest / "plan.json", dict(point=point, added=sorted(extra), knots=sorted(base | extra)))
    assert additions[0].isdisjoint(additions[1])
    assert set.union(*additions) == set(plan["added"])


if __name__ == "__main__":
    freeze(Path(sys.argv[1]))
