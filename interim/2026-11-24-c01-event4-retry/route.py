"""Run a bounded C-01 comparison with the approved mesh budget."""

import argparse
from pathlib import Path
from types import SimpleNamespace

from c01_retry_fields import Fields
from run_case import execute


def run(folder, n, side, mode, seeds, binary, full_fields=True):
    path = folder / f"{side}-{n}-{mode}"
    warm = mode in ("ws", "local-ws")
    options = SimpleNamespace(
        folder=path,
        steps=n,
        solver=side,
        active=warm or mode in ("local-active", "friction-events"),
        arc=False,
        seed=seeds / f"{side}-{n // 2}-base" if warm else None,
        save_seed=False,
        attempt_budget_policy="mesh_proportional",
        history=True,
        localize_friction_events=mode == "friction-events",
        error_scaled_retry=mode == "retry",
        contact_history_bisection=mode in ("local", "local-ws", "local-active", "friction-events"),
        binary=binary,
    )
    if full_fields:
        with Fields(folder / (path.name + "-fields"), n, side) as fields:
            execute(options, fields)
    else:
        execute(options)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder", type=Path)
    p.add_argument("--steps", type=int, choices=(256, 512, 1024), required=True)
    p.add_argument("--solver", choices=("native", "independent"), required=True)
    p.add_argument(
        "--mode", choices=("base", "ws", "retry", "local", "local-ws", "local-active", "friction-events"), required=True
    )
    p.add_argument("--seeds", type=Path, default=Path("seeds"))
    p.add_argument("--binary", type=Path, default=Path("target/release/examples/c01_methods"))
    p.add_argument("--counts-only", action="store_true")
    a = p.parse_args()
    run(a.folder, a.steps, a.solver, a.mode, a.seeds, a.binary, not a.counts_only)
