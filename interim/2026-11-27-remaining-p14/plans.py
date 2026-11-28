"""Freeze the two disjoint point-26 insertions before solving any path."""
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
P13 = HERE.with_name("2026-11-27-remaining-p13")
P12 = HERE.with_name("2026-11-27-remaining-p12")
P11 = HERE.with_name("2026-11-26-remaining-p11")
ROUTES = ("middle", "outer")


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def freeze(out):
    plan = json.loads((P12 / "c01/alternate/plan.json").read_text())
    base = set(plan["p11_knots"])
    changes = plan["points"]["26"]["open_to_slip"]
    mid, = [set(r["span"]) for r in changes if r["level"] == 512]
    outer = {min(r["span"][0] for r in changes), max(r["span"][1] for r in changes)}
    full = json.loads((P13 / "c01/alternate/point26/plan.json").read_text())
    assert len(mid | outer) == 4 and mid | outer == set(full["added"])
    for name, extra in zip(ROUTES, (mid, outer), strict=True):
        dest = out / name
        dest.mkdir(parents=True, exist_ok=True)
        save(dest / "plan.json", dict(point=26, route=name, added=sorted(extra), knots=sorted(base | extra)))


if __name__ == "__main__":
    freeze(Path(sys.argv[1]))
