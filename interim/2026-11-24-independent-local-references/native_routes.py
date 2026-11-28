"""Prepare the three existing C-03 cases and run the production observation probe."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "2026-12-18-c03-plastic-bulk-cause-diagnosis"))
from c03diag_study import configuration, initialize


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path)
    p.add_argument("probe", type=Path)
    a = p.parse_args()
    initialize(a.output, a.probe.resolve(), "tet10", 30.0, open_factor=0.12)
    for name, m, c in (
        ("m0c1", 1e-4, 5e-5),
        ("m1c0", 5e-5, 1e-4),
        ("m1c1", 5e-5, 5e-5),
    ):
        path = a.output / (name + ".json")
        path.write_text(json.dumps(configuration(a.output, 32, m, c)), encoding="utf8")
        subprocess.run(
            [
                a.probe.resolve(),
                a.output / "mesh.json",
                a.output / "ground.json",
                path,
                a.output / name,
                "--diagnose",
            ],
            check=True,
        )
