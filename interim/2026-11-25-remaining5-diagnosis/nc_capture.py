"""Build an isolated copy of the existing NC driver with rejected-system logging."""

from pathlib import Path
import json
import os
import shutil
import subprocess
import sys
import zipfile
import hashlib

root = Path.cwd().resolve()
scratch = Path(sys.argv[1]).resolve()
scratch.mkdir()
crate = scratch / "observer"
crate.mkdir()
for name in ("Cargo.toml", "Cargo.lock", "rust-toolchain.toml"):
    shutil.copyfile(root / name, crate / name)
shutil.copytree(root / "src", crate / "src")
(crate / "examples").mkdir()
shutil.copyfile(
    root / "examples/nc_history_trace.rs", crate / "examples/nc_history_trace.rs"
)
source = crate / "src/stages/dense_residual.rs"
code = source.read_text()
marker = "        if let Some(absolute) = absolute {"
assert code.count(marker) == 1
source.write_text(
    code.replace(
        marker,
        marker
        + """
            if let Ok(folder) = std::env::var("NC_OBSERVER_DIR") {
                let file = std::path::Path::new(&folder).join(format!("reject-{:016x}.json", residual.to_bits()));
                let out = std::fs::File::create_new(file).expect("unique rejection record");
                serde_json::to_writer(std::io::BufWriter::new(out), &serde_json::json!({
                    "a_column_major": a.as_slice(), "b": rhs.as_slice(), "x": x.as_slice(),
                    "absolute": absolute, "relative": residual, "initial": initial
                })).expect("write rejected system");
            }""",
    )
)
inputs = scratch / "inputs"
subprocess.run(
    [
        sys.executable,
        root / "interim/2026-11-22-impl-linear-solver/prepare_inputs.py",
        inputs,
        "--cloud-original",
    ],
    check=True,
)
evidence = root / "interim/2026-11-23-nc01-coordinate-polish"
metadata = json.loads((evidence / "measurements.json").read_bytes())
for name, sha in metadata["input_sha256"].items():
    assert hashlib.sha256((inputs / name).read_bytes()).hexdigest() == sha
with zipfile.ZipFile(evidence / "evidence.zip") as archive:
    settings = archive.read("stages-coordinate_polish.json")
assert hashlib.sha256(settings).hexdigest() == metadata["selected_case_sha256"]
(inputs / "coordinate.json").write_bytes(settings)
dumps = scratch / "rejects"
dumps.mkdir()
env = dict(
    os.environ, CARGO_TARGET_DIR=str(scratch / "target"), NC_OBSERVER_DIR=str(dumps)
)
subprocess.run(
    [
        "cargo",
        "build",
        "--release",
        "--locked",
        "--offline",
        "--example",
        "nc_history_trace",
    ],
    cwd=crate,
    env=env,
    check=True,
)
subprocess.run(
    [
        scratch / "target/release/examples/nc_history_trace",
        "direct",
        inputs / "mesh.json",
        inputs / "ground.json",
        inputs / "coordinate.json",
        scratch / "run",
    ],
    env=env,
    check=True,
    timeout=23 * 3600,
)
