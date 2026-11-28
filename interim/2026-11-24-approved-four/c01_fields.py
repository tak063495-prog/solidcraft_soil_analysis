"""Check six completed or budget-limited histories at their shared load coordinates."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "2026-11-24-c01-event4-retry"))
from collect import History, check_dissipation_anchors, compare, inspect


def measure(raw, mode):
    histories, routes, streams = {}, {}, {}
    for side in ("native", "independent"):
        for n in (256, 512, 1024):
            key = f"{side}-{n}-{mode}"
            folder = raw / key
            r = inspect(folder, side, n, mode, None, None)
            h = History(raw / (key + "-fields"), side, n, r["stages"])
            r["confirmed_dissipation_anchors"] = check_dissipation_anchors(h, folder / "dissipation.jsonl")
            routes[key] = r
            histories[side, n, mode] = h
            streams[key] = h.meta
    return dict(routes=routes, streams=streams, **compare(histories, modes=(mode,)))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("raw", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--mode", choices=("friction-events", "local-active","resolved-slip"), default="friction-events")
    args = p.parse_args()
    result = measure(args.raw, args.mode)
    with args.output.open("x", encoding="utf-8") as f:
        json.dump(result, f, allow_nan=False, indent=2)
        f.write("\n")
