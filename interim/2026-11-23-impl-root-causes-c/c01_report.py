"""Read the independent and native compact trial traces without running a solver."""
import argparse
import gzip
import json
from pathlib import Path

from c01 import PHASES, ROUTES, compare, summarize


def report(source):
    result = {}
    for steps in (512, 1024):
        for solver, prefix in (("native", "native-"), ("independent", "")):
            runs = []
            terminals = {}
            for route in ROUTES:
                with gzip.open(source/f"{prefix}{steps}-{route}.json.gz", "rt", encoding="utf-8") as stream:
                    run = json.load(stream)
                if (run["steps"],run["route"]) != (steps,route):
                    raise ValueError("unexpected trace identity")
                runs.append(run)
                stages = {}
                for stage in PHASES:
                    rows = run["trials"][stage]
                    final = run["last_accepted"].get(stage)
                    target = next((r["end"] for r in rows if r["start"]==final and not r["accepted"]), None)
                    point = None
                    if stage == "open" and steps == 1024:
                        point = 44/1024
                    if stage == "reload" and steps == 512:
                        point = 472/512
                    seen = next((i+1 for i,r in enumerate(rows) if r["accepted"] and r["end"]==point), None)
                    stages[stage] = dict(last_accepted=final, first_rejected_after_final=target,
                                         last_rejected_trial=run["last_rejected"].get(stage),
                                         observation=point, attempts_to_observation=seen,
                                         full_observed_prefix=summarize(rows))
                terminals[route] = dict(stop=run["stop"], attempts=run["counts"], stages=stages)
            result[f"{solver}/{steps}"] = dict(terminals=terminals, common_prefix=compare(runs))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = report(args.source)
    with args.output.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
