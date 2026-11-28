"""Replay the archived independent surface solver on the native accepted time grid."""
import inspect
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
P3 = HERE.parent / "2026-11-25-remaining-p3/c02/p4"
sys.path.insert(0, str(HERE.parent / "2026-11-25-remaining5-p2"))
import c02_face


def run(output, divisions):
    table = np.load(P3 / "increments.npz")["phase_t_fullforce_linear_yield_iterations"]
    grids = [table[table[:, 0] == phase, 1] for phase in range(7)]
    cost = json.loads((P3 / "cost.json").read_text())
    assert cost["degree"] == 4 and cost["base_steps"] == 512
    assert [len(grid) for grid in grids] == cost["attempts_per_stage"]
    for grid in grids:
        assert (np.diff(grid) > 0).all() and grid[-1] == 1 and len(grid) <= 2065
    source = inspect.getsource(c02_face.run)
    uniform = "grids = [np.arange(1,n+1)/n for _ in range(7)]"
    assert source.count(uniform) == 1 and divisions in (4, 8)
    namespace = dict(vars(c02_face), matched=grids)
    exec(compile(source.replace(uniform, "grids = matched"), c02_face.__file__, "exec"), namespace)
    namespace["run"](output, 4, divisions, 512)


if __name__ == "__main__":
    run(Path(sys.argv[1]), int(sys.argv[2]))
