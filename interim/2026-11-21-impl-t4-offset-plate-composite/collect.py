"""Collect the rerun evidence. Preserve incomplete B-02 and B-03 judgments."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def collect(root, output, patch):
    output.mkdir(parents=True, exist_ok=True)
    records = {"passed": False, "b02_complete": False, "b03_complete": False}
    records["dry_patch"] = read(patch / "patch-results.json")
    records["solver_sha256"] = read(patch / "run.json")["solver_binary_sha256"]
    records["input_sha256"] = {name: digest(patch / name) for name in ("mesh.json", "case.json")}
    for suite in ("b02", "b03"):
        data = read(root / suite / "results.json")
        rows = data["fields"]
        expected_count = {"b02": 8, "b03": 16}[suite]
        if len({row["case"] for row in rows}) != expected_count or len(rows) != expected_count:
            raise ValueError(f"incomplete or duplicate original paths for {suite}")
        if any(row["tet10"] != 12 for row in rows):
            raise ValueError(f"unexpected mesh inventory for {suite}")
        scope = data["scope"]
        mismatched = [
            row["name"] for row in scope
            if row.get("matched_expected") is False
            or (suite == "b03" and not row["accepted"] and not row["refused_as_expected"])
        ]
        records[suite] = {
            "fields": rows,
            "compared_states": sum(row["accepted_points"] for row in rows),
            "restart_count": sum(len(row["checkpoint_indices"]) for row in rows),
            "subsequent_states_compared": sum(row["later_states_compared"] for row in rows),
            "bounded_field_paths_within_original_gates": sum(row["bounded_gates_passed"] for row in rows),
            "scope_mismatches": mismatched,
            "scope_counts": {
                "total": len(scope),
                "refused_with_expected_diagnostic": sum(
                    (row.get("refused") is True and row.get("matched_expected") is True)
                    or row.get("refused_as_expected") is True for row in scope
                ),
                "completed": sum(row.get("complete") is True for row in scope),
                "incomplete_output": sum(row.get("complete") is False for row in scope),
            },
        }
    supp = read(root / "b03/supplemental/supplemental.json")
    records["b03"]["supplemental"] = {
        "cases": len(supp["records"]),
        "refusal_contract_passed": supp["refusal_contract_passed"],
        "legacy_controls_completed": supp["legacy_controls_completed"],
    }
    archive = output / "evidence.zip"
    with zipfile.ZipFile(archive, "x", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zipped:
        for name, folder in (("patch", patch), ("b02", root / "b02"), ("b03", root / "b03")):
            for file in sorted(folder.rglob("*")):
                if file.is_file():
                    zipped.write(file, name + "/" + file.relative_to(folder).as_posix())
    records["evidence_sha256"] = digest(archive)
    with (output / "measurements.json").open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(records, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--patch", type=Path)
    args = parser.parse_args()
    collect(args.root, args.output, args.patch or args.root / "patch")
