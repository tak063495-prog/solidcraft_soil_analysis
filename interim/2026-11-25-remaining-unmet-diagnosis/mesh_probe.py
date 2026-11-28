"""Call the archived specimen function with extra material widths in a temporary test binary."""
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

root = Path(sys.argv[1]).resolve()
output = Path(sys.argv[2]).resolve()
output.mkdir()
source = root / "examples/t9d_mcc/main.rs"
reference_dir = root / "interim/2026-11-24-independent-local-references"
sys.path.insert(0, str(reference_dir))
from mcc_history import FIELDS, differences

driver = """
#[test]
fn measure_material_widths() -> Trial<()> {
    let cfg: Value = serde_json::from_str(SPEC)?;
    let mut rows = Vec::new();
    for n in [1, 2] {
        for h in [2.5e-5, 1.25e-5] {
            rows.push(specimen(40., n, h, &cfg["mcc"], true, true, Some(ElasticPath::Analytic))?);
        }
    }
    std::fs::write(std::env::var("MCC_DIAGNOSTIC_OUTPUT")?, serde_json::to_vec(&rows)?)?;
    Ok(())
}
"""
with tempfile.TemporaryDirectory() as tmp:
    temp = Path(tmp)
    entry = temp / "examples/t9d_mcc/main.rs"
    entry.parent.mkdir(parents=True)
    (temp/"interim").symlink_to(root/"interim", target_is_directory=True)
    entry.write_bytes(source.read_bytes() + driver.encode())
    manifest = subprocess.check_output(
        ["cargo", "build", "--offline", "--locked", "--release", "--lib", "--message-format=json"],
        cwd=root, text=True)
    libs = {}
    for line in manifest.splitlines():
        item = json.loads(line)
        if item["reason"] == "compiler-artifact":
            for file in item["filenames"]:
                if file.endswith(".rlib"):
                    libs[item["target"]["name"]] = file
    target = Path(libs["soil_fem_3d"]).parent
    cmd = ["rustc","--test","--edition=2024","-O","-Dwarnings",str(entry),
           "-L","dependency="+str(target/"deps")]
    for name in ("soil_fem_3d","serde_json","nalgebra"):
        cmd.extend(["--extern", name+"="+libs[name]])
    binary = temp/"mesh-probe"
    subprocess.run([*cmd,"-o",str(binary)],check=True)
    full = output/"mesh-histories.json"
    subprocess.run([str(binary),"--exact","measure_material_widths"],check=True,
                   env=dict(os.environ,MCC_DIAGNOSTIC_OUTPUT=str(full)))

reference = np.array(json.loads(gzip.decompress((reference_dir/"mcc-history.json.gz").read_bytes()))["fine_reference"]["40"])
summary = []
for row in json.loads(full.read_bytes()):
    history = row["history"]
    assert len(history)==512 and all(
        r["stage"]=="shear" and r["fraction"]==i/512 for i,r in enumerate(history,1))
    values = np.array([[r[k] for k in FIELDS] for r in history])
    metrics = {k:row[k] for k in ("elements","material_substep","cell_substeps",
                                   "control_substeps","seconds","max_free_relative","max_lateral_error_kpa")}
    metrics["maxima_fraction"] = differences(values,reference)
    summary.append(metrics)
(output/"mesh.json").write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n")
