"""Compare final CLI and archival B-03 physical states with verified producer IDs."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
B03 = ROOT / "interim/2026-12-02-tri6-offset-kink-b03"
sys.path.insert(0, str(ROOT / "interim/2026-12-01-tri6-offset-composite-b02"))
from study import exact


def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def check(cli, formal, out):
    out.mkdir(parents=True, exist_ok=False)
    read = lambda p: json.loads(p.read_bytes())
    manifests = {
        r["path"]: r["sha256"] for r in read(B03 / "evidence-manifest.json")["members"]
    }
    original = {r["case"]: r for r in read(B03 / "measurements.json")["cases"]}
    runs = read(formal / "measurements.json")["cases"]
    assert set(original) == {r["case"] for r in runs} and len(runs) == 21
    cli_sha = digest(cli)
    observer_sha = digest(next((formal / "solver").iterdir()))
    rows = []
    for r in runs:
        name = r["case"]
        folder = formal / name
        inputs = [folder / (k + ".json") for k in ("mesh", "ground", "stages")]
        fingerprints = {p.name: digest(p) for p in inputs}
        assert fingerprints == original[name]["input_sha256"] == r["input_sha256"]
        file = out / (name + ".json")
        subprocess.run(
            [cli, "solve-stages", *inputs, file], check=True, stdout=subprocess.DEVNULL
        )
        a = read(file)
        b = read(folder / "plain/run.json")
        # Unlike observer plain vs observer solve, these are different binaries.
        # Authenticate both identity fields before comparing every other field.
        assert a["provenance"].pop("solver_binary_sha256") == cli_sha
        assert b["provenance"].pop("solver_binary_sha256") == observer_sha
        assert exact(a, b), name

        def same(rel):
            return digest(folder / rel) == manifests["fresh/" + name + "/" + rel]

        rows.append(
            dict(
                case=name,
                input_sha256=fingerprints,
                accepted=r["comparison"]["accepted_points"],
                cli_and_observer_equal_after_producer_check=True,
                native_accepted_identical_to_archived=same("native/accepted.jsonl"),
                independent_identical_to_archived=same("independent.jsonl"),
            )
        )
    (out / "results.json").write_text(
        json.dumps(
            dict(cli_sha256=cli_sha, observer_sha256=observer_sha, rows=rows), indent=2
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("cli", type=Path)
    p.add_argument("formal", type=Path)
    p.add_argument("out", type=Path)
    args = p.parse_args()
    check(args.cli.resolve(), args.formal, args.out)
