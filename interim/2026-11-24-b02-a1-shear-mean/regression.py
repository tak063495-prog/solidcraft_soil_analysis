"""Same-host parent/updated binaries, identical B-01 and A3 inputs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np

from study import relative
from cylinders import cylinder_input, standard

SUFFIX = ".exe" if os.name == "nt" else ""


def shell(cli, case, dof, mode):
    with tempfile.TemporaryDirectory() as temp:
        folder = Path(temp)
        source, output = folder/"input.json", folder/"output.json"
        source.write_text(json.dumps(case, allow_nan=False), encoding="utf-8")
        args = [cli, "curved-shell", source, output]
        if mode == "element_mean":
            args.append("--shear-mean")
        run = subprocess.run(args, capture_output=True, text=True)
        if run.returncode:
            return dict(diagnostic=run.stderr.strip(), displacement_m=None)
        data = json.loads(output.read_bytes())
        # The binary hash necessarily differs. Hash every frame including
        # configurations, reactions, energies, soil fields and shell resultants.
        frames = data["frames"]
        digest = hashlib.sha256()
        for frame in frames:
            digest.update(json.dumps(frame, separators=(",", ":"), allow_nan=False).encode())
        return dict(diagnostic=None, displacement_m=frames[-1]["configuration"][dof],
                    input_sha256=data["input_sha256"], frames_sha256=digest.hexdigest())


def benchmark_cases():
    for card in json.loads(standard.CARDS.read_bytes())["benchmarks"]:
        probe = card["probe"]
        for distorted in (False, True):
            for n in (2, 4, 8):
                case = standard.build(card, n, distorted)
                xyz = np.array(case["nodes"])
                (node,) = np.flatnonzero(np.linalg.norm(xyz-probe["xyz"], axis=1) < probe["coordinate_tolerance"])
                dof = 6*node + ("ux", "uy", "uz").index(probe["component"])
                yield (card["id"], n, distorted), case, dof
    for name in ("b01-regular16", "b01-edf364"):
        case, node = cylinder_input(name)
        yield (name,), case, 6*node+2


def regression(parent, current):
    rows = []
    for key, case, dof in benchmark_cases():
        old, new = (shell(release/("soil-fem-3d"+SUFFIX), case, dof, "centroid")
                    for release in (parent, current))
        row = dict(case=key, parent_default=old, current_default=new, identical=(old == new))
        if len(key) == 1:
            row["opt_in"] = shell(current/("soil-fem-3d"+SUFFIX), case, dof, "element_mean")
            value = row["opt_in"]["displacement_m"]
            row["opt_in_cylinder_passed"] = value is not None and relative(value, -1.8248e-5, .03)["passed"]
        rows.append(row)
        print(key, row["identical"], flush=True)
    plates = []
    for grid in ([], ["12x6"], ["18x9"]):
        outputs = [subprocess.check_output([release/"examples"/("t10_layered_cantilever"+SUFFIX), *grid])
                   for release in (parent, current)]
        tips = [json.loads(output)["tip_dz_m"] for output in outputs]
        plates.append(dict(grid=next(iter(grid), "10x5"), parent_a3=tips[0], current_a3=tips[1],
                           identical=(outputs[0] == outputs[1]),
                           stdout_sha256=[hashlib.sha256(v).hexdigest() for v in outputs]))
    return dict(parent_commit="7a47af617041bf8a409baf77e0602af6cf1fbac0", b01=rows, b02=plates)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parent_release", type=Path)
    parser.add_argument("release", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = regression(args.parent_release.resolve(), args.release.resolve())
    with args.output.open("x", encoding="utf-8", newline="\n") as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write("\n")
