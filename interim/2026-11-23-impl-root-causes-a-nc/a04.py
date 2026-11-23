"""Compare a flooded soil cut's boundary transfer with the input-only RT0/elastic model."""
import argparse
import copy
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "2026-11-29-richards-wet-joint-a04"))
import study


def solve(cli, folder, case, inputs):
    path = folder / "case.json"
    study.write(path, case)
    return subprocess.run([cli, "solve-one-way-construction", *inputs, path, folder / "native.json"],
                          capture_output=True, text=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cli", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    folder = args.output
    folder.mkdir(parents=True, exist_ok=False)
    mesh, ground, case = study.fixture(4, True)
    inputs = [folder / "mesh.json", folder / "ground.json"]
    for path, value in zip(inputs, (mesh, ground), strict=True):
        study.write(path, value)
    cli = args.cli.resolve()
    result = solve(cli, folder, case, inputs)
    (folder / "native.log").write_text(result.stdout + result.stderr, encoding="utf-8")
    result.check_returncode()
    study.independent(folder)
    study.native_fields(folder)
    comparison = study.compare(folder / "native-fields.jsonl", folder / "independent-fields.jsonl")
    native_rows = [row for row, _ in study.fields(folder / "native-fields.jsonl")]
    ref_rows = [row for row, _ in study.fields(folder / "independent-fields.jsonl")]
    assert native_rows == ref_rows and len({study.exact(r) for r in native_rows}) == len(native_rows)
    all_fields = comparison["all_fields"]
    assert all(v.get("mismatches", 0) == 0 for v in all_fields.values())
    for k in ("u", "stress"):
        assert all_fields[k]["relative_l2"] < 1e-7
    for k in ("effective_internal", "fluid_internal", "water_body", "exterior", "physical_load",
              "equivalent_rhs", "applied_load"):
        assert all_fields[k]["maximum_absolute"] < 1e-10
    legacy = copy.deepcopy(case)
    for bc in legacy["one_way"]["flow"]["analysis"]["periods"][2]["boundaries"]:
        bc["condition"]["value"]["head_m"] = 1.8
    for p in legacy["one_way"]["exterior_periods"][2:]:
        p["reservoirs"][0]["level_m"] = 2.2
    control = folder / "legacy"
    control.mkdir()
    rejected = solve(cli, control, legacy, inputs)
    assert rejected.returncode == 1 and not (control / "native.json").exists()
    assert "exterior reservoir and prescribed hydraulic head disagree" in rejected.stderr
    study.write(folder / "results.json", dict(limited_pass=True, states=len(native_rows),
                 fields=len(all_fields), comparisons=all_fields, old_case_exit=rejected.returncode,
                 old_case_message=rejected.stderr.strip()))


if __name__ == "__main__":
    main()
