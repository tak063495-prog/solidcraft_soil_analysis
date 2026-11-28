"""Build a disposable crate with two read-only probes in the native trial loop."""

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DIAG = HERE.parent / "2026-11-25-remaining5-diagnosis"
NC = HERE.parent / "2026-11-23-nc01-coordinate-polish"
SCRATCH = Path(sys.argv[1]).resolve()
MODE = sys.argv[2]
CRATE = SCRATCH / "crate"
ENV = dict(os.environ, CARGO_TARGET_DIR=str(SCRATCH / "target"))
if MODE == "prepare":
    SCRATCH.mkdir()
    CRATE.mkdir()
    for name in ("Cargo.toml", "Cargo.lock", "rust-toolchain.toml"):
        shutil.copy(ROOT / name, CRATE / name)
    shutil.copytree(ROOT / "src", CRATE / "src")
    (CRATE / "examples").mkdir()
    shutil.copy(ROOT / "examples/nc_history_trace.rs", CRATE / "examples")
    shutil.copy(DIAG / "c01_reload.rs", CRATE / "examples")
    source = CRATE / "src/stages.rs"
    code = source.read_text()
    inputs = SCRATCH / "inputs"
    subprocess.run(
        [
            sys.executable,
            HERE.parent / "2026-11-22-impl-linear-solver/prepare_inputs.py",
            inputs,
            "--cloud-original",
        ],
        check=True,
    )
    metadata = json.loads((NC / "measurements.json").read_bytes())
    for name, expected in metadata["input_sha256"].items():
        assert hashlib.sha256((inputs / name).read_bytes()).hexdigest() == expected
    with zipfile.ZipFile(NC / "evidence.zip") as z:
        raw = z.read("run.json")
        assert hashlib.sha256(raw).hexdigest() == metadata["run_sha256"]
        run = json.loads(raw)
        selected = z.read("stages-coordinate_polish.json")
        assert hashlib.sha256(selected).hexdigest() == metadata["selected_case_sha256"]
        (inputs / "coordinate.json").write_bytes(selected)
    targets = [
        a for a in run["attempts"] if a["stage"] == "excavate_upper" and a["accepted"]
    ][-5:]
    targets.append(run["attempts"][-1])
    targets = [(a["progress_from"], a["progress_to"]) for a in targets]
    code_targets = ", ".join(f"({a!r}_f64,{b!r}_f64)" for a, b in targets)
    trigger = "            let relative = reduced.norm() / scale;"
    assert code.count(trigger) == 1
    probe = (HERE / "capture_hook.rs.inc").read_text().replace("NC_PAIRS", code_targets)
    source.write_text(code.replace(trigger, trigger + "\n" + probe))
    subprocess.run(
        [
            "cargo",
            "build",
            "--offline",
            "--release",
            "--locked",
            "--example",
            "c01_reload",
            "--example",
            "nc_history_trace",
        ],
        env=ENV,
        cwd=CRATE,
        check=True,
    )
elif MODE == "c01":
    with (SCRATCH / "c01-route.json").open("x") as out:
        subprocess.run(
            [
                SCRATCH / "target/release/examples/c01_reload",
                ROOT / "examples/joint-water-column-tet10-mesh.json",
                ROOT / "examples/ground-joint-water-mc.json",
                DIAG / "c01-case.json",
            ],
            cwd=SCRATCH,
            stdout=out,
            check=True,
            timeout=23 * 3600,
        )
elif MODE == "nc":
    inputs = SCRATCH / "inputs"
    subprocess.run(
        [
            SCRATCH / "target/release/examples/nc_history_trace",
            "direct",
            inputs / "mesh.json",
            inputs / "ground.json",
            inputs / "coordinate.json",
            SCRATCH / "nc-run",
        ],
        cwd=SCRATCH,
        check=True,
        timeout=23 * 3600,
    )
else:
    raise ValueError(MODE)
