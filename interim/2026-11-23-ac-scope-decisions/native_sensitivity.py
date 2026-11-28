"""Repeat the three C-03 native paths with the original partitions and with a stricter force gate."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

STUDY = Path(__file__).resolve().parents[1] / "2026-12-18-c03-plastic-bulk-cause-diagnosis"
sys.path.insert(0, str(STUDY))
from c03diag_study import binary


def run(root):
    probe = binary(root)
    for route in ("m1c1", "m0c1", "m1c0"):
        case = root / (route + "-native") / "stages.json"
        inputs = [root / "mesh.json", root / "ground.json"]
        # The ordinary c03_cloud entry records states; --diagnose also records
        # the accepted quadrature partitions and the material strain increments.
        subprocess.run(
            [probe, *inputs, case, root / (route + "-trace"), "--diagnose"], check=True
        )
        original = json.loads(case.read_bytes())
        assert original["settings"]["relative_force_tolerance"] == 1e-8
        original["settings"]["relative_force_tolerance"] = 1e-10
        dest = root / (route + "-tight.json")
        with dest.open("x") as f:
            json.dump(original, f)
        subprocess.run([probe, *inputs, dest, root / (route + "-tight")], check=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    run(parser.parse_args().root.resolve())
